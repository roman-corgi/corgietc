import os
import re
import warnings
from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
import astropy.constants as const
import numpy as np
import scipy.interpolate
from scipy.optimize import minimize, root_scalar
from tqdm import tqdm

from EXOSIMS.OpticalSystem.Nemati import Nemati
from EXOSIMS.util._numpy_compat import copy_if_needed
from cgi_noise import cginoiselib as fl
from cgi_noise.tsnr_core import corePhotonRates

import yaml
from eetc.cgi_eetc import CGIEETC
import eetc
from astropy.table import Table
from eetc.load import load_excam_config, load_flux_grid, _unpack_excam_config
from eetc.excam_tools import _ENF

eetc_root = eetc.lib_dir
pointer_path_sci = os.path.join(eetc_root,'pointer_cont_frames.yaml')
pointer_path_cal = os.path.join(eetc_root,'pointer.yaml')
if os.path.exists(pointer_path_sci):
    print('specified pointer file exists!\n')

with open(pointer_path_sci, 'r') as file:
    pointer_data = yaml.safe_load(file)

for fi in pointer_data:
    config_fi = os.path.join(eetc_root, pointer_data[fi])

excam_config_file = os.path.join(eetc_root, pointer_data['excam_config'])

with open(excam_config_file, 'r') as file:
    excam_data = yaml.safe_load(file)

sequence_file = os.path.join(eetc_root,pointer_data['sequences'])

with open(sequence_file, 'r') as file:
    data = yaml.safe_load(file)

rows = []

for config_name, params in data.items():
    if isinstance(params, dict):
        row = {'config_name': config_name}
        row.update(params)
        rows.append(row)
config_table = Table(rows=rows)
config_table.add_index('config_name')

flux_grid_file = os.path.join(eetc_root, pointer_data['flux_grid'])
eetc_flux_grid = load_flux_grid(flux_grid_file)
eetc_valid_spts_v = list(eetc_flux_grid['v'][1])

_SPT_CLASS_ORDER = "OBAFGKM"

def rough_wall_clock(n_frames,
                     exp_t,
                     mode='continuous',
                     readout_t=0.264,
                     transfer_time=3.):
    """
    Compute the rough wall clock time for a given observation sequence.
    
    Computes the total exposure time including readout/transfer overheads
    for a given observation sequence.

    Parameters
    ----------
    n_frames : int or float
        Number of exposures.
    exp_t : float
        Individual frame exposure time in seconds.
    mode : str, optional
        Observation mode - if set to 'continuous' (default), the
        total wall clock time is computed assuming continuous reads 
        with frame transfer at the end. This is the standard for 'science'
        observations. If not set to continuous, it is assumed that each
        frame takes a frame transfer time hit. This is generally the 
        assumed mode for unocculted calibration observations.
    readout_t : float, optional
        Single EXCAM exposure readout time in seconds. Default is 0.264.
    transfer_time : float, optional
        EXCAM Science frame transfer time in seconds. Default is 3.0.

    Returns
    -------
    float
        Total wall clock time in seconds.

    Notes
    -----
    EXCAM has a 0.264 second read time for SCI-sized frames.
    In continuous exposure mode, EXCAM frames can only be transferred once
    every 3 seconds due to data rate limitations. Very short exposures
    will have large overheads!
    
    EETC has a more accurate version built in for fixed-integration-time
    calculation for using EXCAM in "burst" mode.
    
    For SCI observations, transfer occurs while exposing (no n_frames * transfer_time penalty). 
    
    For CAL observations, frames are read out after each exposure and transferred in bulk at the end
    (includes n_frames * readout_time and n_frames * transfer_time penalty).
    """
    # SCI observations are reading out while exposing
    # (no n_frames * readout_time penalty)
    if mode == 'continuous':
        if exp_t < transfer_time:
            frame_t = transfer_time + readout_t 
        else:
            frame_t = exp_t + readout_t
        return frame_t * n_frames + transfer_time
    # CALIBRATION observations take the n_frames * (readout_time + transfer_time) penalty
    else:
        frame_t = exp_t + readout_t
        return n_frames * (frame_t + transfer_time)
    
def _parse_spt(spt):
    """Roughly parse a spectral type string into (class_index, subclass,
    luminosity_class) for nearest-match purposes.
 
    This is a best-effort parser, not a rigorous spectral classification
    tool, it exists only to find the closest available entry in EETC's
    discrete flux grid for a given TargetList spectral type string.
    """
    m = re.match(r"\s*([OBAFGKM])\s*(\d+(?:\.\d+)?)?\s*([IV]+.*)?", spt.upper())
    if not m:
        return (len(_SPT_CLASS_ORDER), 5.0, "V")
    cls_letter, subclass, lumclass = m.groups()
    cls_idx = (
        _SPT_CLASS_ORDER.index(cls_letter)
        if cls_letter in _SPT_CLASS_ORDER
        else len(_SPT_CLASS_ORDER)
    )
    sub = float(subclass) if subclass else 5.0
    return (cls_idx, sub, (lumclass or "V").strip())

def match_spt_to_eetc_grid(spt, valid_spts=eetc_valid_spts_v):
    """Find the closest available EETC flux-grid spectral type to `spt`.
 
    EETC's CGIEETC class requires an exact (case-insensitive) match against
    its internal flux grid's spectral-type list (`valid_spts`). Most
    catalog/TargetList spectral types will not match exactly.  This performs
    a simple nearest-neighbor match in (spectral-class, subclass) space,
    preferring an exact luminosity-class match where available and falling
    back to any luminosity class otherwise.
 
    Args:
        spt (str):
            Spectral type string to match, e.g. "G2V", "K3.5 III".
        valid_spts (list(str)):
            Candidate spectral types to match against (default:
            `eetc_valid_spts_v`, the types available for phot="v").
 
    Returns:
        str: The closest entry in `valid_spts` (verbatim, as found in the
        grid -- pass this directly to `CGIEETC(..., spt=...)`).
    """
    target = _parse_spt(spt)
    parsed = [(_parse_spt(s), s) for s in valid_spts]
 
    # prefer matches with the same luminosity class, if any exist
    same_lum = [p for p in parsed if p[0][2] == target[2]]
    candidates = same_lum if same_lum else parsed
 
    def _dist(p):
        (cls_idx, sub, _), _ = p
        return abs((cls_idx * 10 + sub) - (target[0] * 10 + target[1]))
 
    return min(candidates, key=_dist)[1]

class corgietc(Nemati):
    r"""corgietc Optical System class

    Optical System Module based on cgi_noise model.

    Args:
        CritLam (float)
            Default critical wavelength (Nyquist sampling) in nm. Only used if not set
            in scienceInstrument input specification definition. Defaults to 500
        compbeamD (float)
            Default compressed beam diameter in m. Only used if not set
            in scienceInstrument input specification definition. Defaults to 0.005
        fnlFocLen (float)
            Default final focal length in m. Only used if not set in scienceInstrument
            input specification definition. Defaults to 0.26
        PSF_x_lamD (float)
            Default PSF core x extent in lam/D. Only used if not set
            in scienceInstrument input specification definition. Defaults to 0.942
        PSF_y_lamD (float)
            Default PSF core y extent in lam/D. Only used if not set
            in scienceInstrument input specification definition.  Defaults to 0.45
        Rlamsq (float)
            Quadratic term of resolving power at PSF model. Only used if not set
            in scienceInstrument input specification definition. Defaults to 0.000854964
        Rlam (float)
            Linear term of resolving power at PSF model. Only used if not set in
            scienceInstrument input specification definition. Defaults to -1.513136232
        Rconst (float)
            Constant term of resolving power at PSF model. Only used if not set
            in scienceInstrument input specification definition. Defaults to 707.8977209
        pp_Factor_CBE (float)
            Post-processing factor (e.g., 30 for 30x speckle suppression). Only used if
            not set in scienceInstrument input specification definition. Defaults to 2.0
        RefStar_SpectralType (str)
            Spectral type of the reference star (eg a0v, b3v, a5v, f5v, g0v, g5v, k0v,
            k5v, m0v, m5v)
        RefStar_V_mag (float)
            Visual Magnitude of the reference star
        TimeonRefStar_tRef_per_tTar (float)
            Time on a reference star per target
        contrast_degradation (float)
            Multiplier for rawcontrast (e.g. 0.5 represents 50% rawcontrast). Defaults
            to 1.0
        desiredRate (float)
            Target value for e-/pix/frame. Defaults to 0.1
        tfmin (float)
            Minimum frame time in seconds. Defaults to 3
        tfmax (float)
            Maximum frame time in seconds. Defaults to 100
        frameThresh (float)
            Threshold value at which to switch from photon counting to analog mode in
            e-/pix/frame.  If the approximated value is above the threshold, analog mode
            is used in calculating frame time and effective QE.  Ignored if
            forcePhotonCounting is set to True. Defaults to 0.5
        forcePhotonCounting (float)
            If True, always use photon counting mode regardless of frame counts.
            Defaults to False
        **specs:
            EXOSIMS input specification dictionary

    Attributes:
        default_vals_extra2 (dict):
            Dictionary of local default values.
        desiredRate (float)
            Target value for e-/pix/frame.
        frameThresh (float)
            Threshold value at which to switch from photon counting to analog mode in
            e-/pix/frame.  If the approximated value is above the threshold, analog mode
            is used in calculating frame time and effective QE.  Ignored if
            forcePhotonCounting is set to True.
        forcePhotonCounting (float)
            If True, always use photon counting mode regardless of frame counts.
        hc (float):
            h * c in m^3 kg s^-2
        radas (float):
            Conversion factor from arcsec to radians
        SPECTRA_Data (cgi_noise.loadCSVrow.loadCSVrow):
            Spectral data for reference stars
        SPECTRA_deltaLambda (float):
            Wavelength step (in m) of SPECTRA_Data
        tfmin (float)
            Minimum frame time in seconds.
        tfmax (float)
            Maximum frame time in seconds.

    """

    def __init__(
        self,
        CritLam=500,
        compbeamD=0.005,
        fnlFocLen=0.26,
        PSF_x_lamD=0.942,
        PSF_y_lamD=0.45,
        Rlamsq=0.000854964,
        Rlam=-1.513136232,
        Rconst=707.8977209,
        pp_Factor_CBE=2.0,
        desiredRate=0.1,
        tfmin=3,
        tfmax=100,
        frameThresh=0.5,
        RefStar_SpectralType="a0v",
        RefStar_V_mag=2.26,
        TimeonRefStar_tRef_per_tTar=0.25,
        contrast_degradation=1.0,
        forcePhotonCounting=False,
        **specs,
    ):

        # useful conversion factors
        self.radas = ((1 * u.arcsec).to(u.rad)).value
        self.hc = (const.h * const.c).to_value(u.m**3 * u.kg / u.s**2)

        # load cgi_noise spec data
        spectra_path = (
            Path(os.environ["CGI_NOISE_DATA_DIR"]) / "Spectra" / "SPECTRA_ALL_BPGS.csv"
        )
        self.SPECTRA_Data = fl.loadCSVrow(spectra_path)
        self.SPECTRA_deltaLambda = (
            self.SPECTRA_Data.df.at[2, "Wavelength_m"]
            - self.SPECTRA_Data.df.at[1, "Wavelength_m"]
        )

        # set frame threshold values
        self.tfmin = tfmin
        self.tfmax = tfmax
        self.desiredRate = desiredRate
        self.frameThresh = frameThresh
        self.forcePhotonCounting = forcePhotonCounting

        self._eetc_star_cache = {}
        self.eetc_excam_config = load_excam_config(excam_config_file)
        # package inputs for use in popoulate*_extra
        self.default_vals_extra2 = {
            "CritLam": CritLam,
            "compbeamD": compbeamD,
            "fnlFocLen": fnlFocLen,
            "PSF_x_lamD": PSF_x_lamD,
            "PSF_y_lamD": PSF_y_lamD,
            "Rlamsq": Rlamsq,
            "Rlam": Rlam,
            "Rconst": Rconst,
            "pp_Factor_CBE": pp_Factor_CBE,
            "RefStar_SpectralType": RefStar_SpectralType,
            "RefStar_V_mag": RefStar_V_mag,
            "TimeonRefStar_tRef_per_tTar": TimeonRefStar_tRef_per_tTar,
            "contrast_degradation": contrast_degradation,
        }

        Nemati.__init__(self, **specs)

        # add local defaults to outspec
        for k in self.default_vals_extra2:
            self._outspec[k] = self.default_vals_extra2[k]

        for k in [
            "desiredRate",
            "tfmin",
            "tfmax",
            "frameThresh",
            "forcePhotonCounting",
        ]:
            self._outspec[k] = getattr(self, k)

    def populate_starlightSuppressionSystems_extra(self):

        # add PSFPeak and contrast stability values
        if "PSFpeak" not in self.allowed_starlightSuppressionSystem_kws:
            self.allowed_starlightSuppressionSystem_kws.append("PSFpeak")

        cstability_params = [
            "AvgRawContrast",
            "ExtContStab",
            "IntContStab",
            "SystematicC",
            "InitStatContrast",
        ]
        for param_name in cstability_params:
            self.allowed_starlightSuppressionSystem_kws.append(param_name)

        for nsyst, syst in enumerate(self.starlightSuppressionSystems):
            # process PSFPeak
            syst = self.get_coro_param(
                syst,
                "PSFpeak",
                expected_ndim=2,
                expected_first_dim=2,
                min_val=0.0,
            )

            for param_name in cstability_params:
                # SystematicC is optional so check for parameter presence in input
                if param_name in syst:
                    # load the data
                    dat, hdr = self.get_param_data(
                        syst[param_name],
                        left_col_name=syst["csv_angsep_colname"],
                        param_name=syst["csv_names"][param_name],
                        expected_ndim=2,
                        expected_first_dim=2,
                    )
                    WA, D = dat[0].astype(float), dat[1].astype(float)

                    # if the first entry is larger than the IWA, update it to the IWA
                    if WA[0] * syst["input_angle_unit_value"] > syst["IWA"]:
                        WA[0] = (
                            (syst["IWA"] / syst["input_angle_unit_value"])
                            .decompose()
                            .value
                        )
                    WA = WA * syst["input_angle_unit_value"]

                    # generate previous entry lookup
                    Dinterp = scipy.interpolate.interp1d(
                        WA,
                        D,
                        kind="previous",
                        fill_value="extrapolate",
                        bounds_error=False,
                    )

                    # create a callable lambda function. for coronagraphs, we need to
                    # scale the angular separation by wavelength
                    syst[param_name] = lambda lam, s, Dinterp=Dinterp, lam0=syst[
                        "lam"
                    ]: np.array(Dinterp((s * lam0 / lam).to_value("arcsec")), ndmin=1)

            # ensure that CGintSamp is in the system
            syst["CGintSamp"] = syst.get("CGintSamp", 0.1)

            # load Throughput Data
            syst["Throughput_Data"] = fl.loadCSVrow(
                os.path.normpath(os.path.expandvars(syst["Throughput_Data"]))
            )

    def populate_scienceInstruments_extra(self):
        """Additional setup for science instruments."""

        # specify dictionary of keys and units
        kws = {
            "CritLam": u.nm,  # critical wavelength
            "compbeamD": u.m,  # compressed beam diameter
            "fnlFocLen": u.m,  # final focal length
            "PSF_x_lamD": None,  # PSF x extent in lambda/D
            "PSF_y_lamD": None,  # PSF y extent in lambda/D
            "Rlamsq": None,
            "Rlam": None,
            "Rconst": None,
        }
        self.allowed_scienceInstrument_kws += list(kws.keys())

        for ninst, inst in enumerate(self.scienceInstruments):

            # load all additional detector specifications
            for kw in kws:
                inst[kw] = float(inst.get(kw, self.default_vals_extra2[kw]))
                self._outspec["scienceInstruments"][ninst][kw] = inst[kw]
                if kws[kw] is not None:
                    inst[kw] *= kws[kw]

            # compute fnumber
            if "imager" in inst["name"].lower():
                pass
            elif "spec" in inst["name"].lower():
                inst["fnumber"] = (
                    (inst["fnlFocLen"] / inst["compbeamD"]).decompose().value
                )
            else:
                raise Exception("Instrument name must contain IMAGER or SPEC")

            # compute pixel scale
            inst["pixelScale"] = (inst["CritLam"] / self.pupilDiam / 2).to(
                u.arcsec, equivalencies=u.dimensionless_angles()
            )

            # load detector and QE data
            inst["DET_CBE_Data"] = fl.loadCSVrow(
                os.path.normpath(os.path.expandvars(inst["DET_CBE_Data"]))
            )
            inst["DET_QE_Data"] = fl.loadCSVrow(
                os.path.normpath(os.path.expandvars(inst["DET_QE_Data"]))
            )
            inst["matrix"] = np.genfromtxt(
                os.path.normpath(os.path.expandvars(inst["matrix"])), delimiter=","
            )

    def populate_observingModes_extra(self):
        """Add specific observing mode keywords"""

        self.allowed_observingMode_kws.append("Scenario")
        self.allowed_observingMode_kws.append("StrayLight_Data")
        self.allowed_observingMode_kws.append("pp_Factor_CBE")
        self.allowed_observingMode_kws.append("RefStar_SpectralType")
        self.allowed_observingMode_kws.append("RefStar_V_mag")
        self.allowed_observingMode_kws.append("TimeonRefStar_tRef_per_tTar")
        self.allowed_observingMode_kws.append("contrast_degradation")

        for nmode, mode in enumerate(self.observingModes):
            assert "Scenario" in mode and isinstance(
                mode["Scenario"], str
            ), "All observing modes must have key 'Scenario'."

            if "imager" in mode["instName"].lower():
                mode["f_SR"] = 1
            elif "spec" in mode["instName"].lower():
                mode["f_SR"] = 1 / (mode["inst"]["Rs"] * mode["BW"])

                # compute mpix:
                pixPerlamD = (
                    (mode["lam"] * mode["inst"]["fnumber"] / mode["inst"]["pixelSize"])
                    .decompose()
                    .value
                )
                xpixPerCor = mode["inst"]["PSF_x_lamD"] * 2 * pixPerlamD
                ypixPerCor = mode["inst"]["PSF_y_lamD"] * 2 * pixPerlamD
                lamnm = mode["lam"].to_value(u.nm)
                ResPowatPSF = (
                    mode["inst"]["Rconst"]
                    + mode["inst"]["Rlam"] * lamnm
                    + mode["inst"]["Rlamsq"] * lamnm**2
                )
                dpix_dlam = ResPowatPSF * xpixPerCor / mode["lam"]
                xpixPerSpec = dpix_dlam * mode["lam"] / mode["inst"]["Rs"]
                mode["mpix"] = xpixPerSpec * ypixPerCor
            else:
                raise Exception("Instrument name must contain IMAGER or SPEC")

            # stray light values
            assert "StrayLight_Data" in mode and isinstance(
                mode["StrayLight_Data"], str
            ), "All observing modes must have key 'StrayLight_Data'."

            mode["stray_ph_s_mm2"] = fl.getStrayLightfromfile(
                mode["Scenario"],
                "CBE",
                fl.loadCSVrow(
                    os.path.normpath(os.path.expandvars(mode["StrayLight_Data"]))
                ),
            )

            mode["stray_ph_s_pix"] = (
                mode["stray_ph_s_mm2"] * (mode["inst"]["pixelSize"].to_value(u.mm)) ** 2
            )

            # generate inBandFlux0_sum object
            lam_m = mode["lam"].to_value(u.m)
            bandRange = self.SPECTRA_Data.df[
                abs(self.SPECTRA_Data.df["Wavelength_m"] - lam_m)
                <= (0.5 * mode["BW"] * lam_m)
            ]
            onlySpec = bandRange.drop(["Wavelength_m", "E_ph_J"], axis=1)

            Ephot = self.hc / lam_m
            onlySpecEphot = onlySpec.apply(
                lambda x: x / Ephot, axis=1, result_type="broadcast"
            )
            mode["inBandFlux0_sum"] = (
                onlySpecEphot.sum(axis=0) * self.SPECTRA_deltaLambda
            )

            # ensure pp_Factor_CBE is in the mode
            mode["pp_Factor_CBE"] = mode.get(
                "pp_Factor_CBE", self.default_vals_extra2["pp_Factor_CBE"]
            )

            # ensure RefStar_SpectralType is in the mode
            mode["RefStar_SpectralType"] = mode.get(
                "RefStar_SpectralType", self.default_vals_extra2["RefStar_SpectralType"]
            )

            # ensure RefStar_V_mag is in the mode
            mode["RefStar_V_mag"] = mode.get(
                "RefStar_V_mag", self.default_vals_extra2["RefStar_V_mag"]
            )

            # ensure TimeonRefStar_tRef_per_tTar is in the mode
            mode["TimeonRefStar_tRef_per_tTar"] = mode.get(
                "TimeonRefStar_tRef_per_tTar",
                self.default_vals_extra2["TimeonRefStar_tRef_per_tTar"],
            )

            # ensure contrast_degradation is in the mode
            mode["contrast_degradation"] = mode.get(
                "contrast_degradation", self.default_vals_extra2["contrast_degradation"]
            )

    def construct_cg(self, mode, WA):
        "Repackage values at a single WA into CGParameters object"

        syst = mode["syst"]
        lam = mode["lam"]

        # NB: syst["core_area"] returns a value in (l/D)^2 units, which needs to be
        # converted to square arcseconds by dividing through by the input angle unit
        # value squared
        omegaPSF = (
            (syst["core_area"](lam, WA) / syst["input_angle_unit_value"] ** 2)
            .decompose()
            .value
        )
        CGintmpix = (
            (
                omegaPSF
                * self.radas**2
                / ((syst["CGintSamp"] * syst["lam"] / self.pupilDiam) ** 2)
            )
            .decompose()
            .value
        )

        WAl = np.repeat(WA, 1)
        PSFpeakI = syst["PSFpeak"](lam, WAl)[0]
        if "IMG_NFB1_HLC" in mode["Scenario"]:
            PSFpeakI /= CGintmpix

        CG_PSFarea_sqlamD = omegaPSF / (syst["lam"].to_value(u.m) / self.radas) ** 2

        out = fl.CGParameters(
            CGcoreThruput=syst["core_thruput"](lam, WAl)[0],
            PSFpeakI=PSFpeakI,
            omegaPSF=omegaPSF,
            CGintSamp=syst["CGintSamp"],
            CGradius_arcsec=None,
            CGdesignWL=lam.to_value(u.m),
            CGintmpix=CGintmpix,
            CG_PSFarea_sqlamD=CG_PSFarea_sqlamD,
            CGintensity=syst["core_mean_intensity"](lam, WAl)[0],
            CG_occulter_transmission=syst["occ_trans"](lam, WAl)[0],
            CGcontrast=syst["core_contrast"](lam, WAl)[0],
        )

        return out
    
    def _get_eetc_sequence_name(self, mode):
        """Resolve the EETC `sequence_name` to use for star-flux calibration
        for a given corgietc observing `mode`.
 
        corgietc's `mode` dict (unlike EETC's own examples) doesn't carry
        the FPAM/FSAM/CFAM/DPAM/SPAM_LSAM optical-configuration keys that
        EETC's `sequences.yaml`/`config_table` are indexed by -- only a
        human-readable `mode["Scenario"]` string. There is no general,
        reliable way to derive one from the other without a lookup table
        that maps corgietc Scenario names to EETC optical configurations
        (I don't have visibility into corgietc's Scenario definitions, so I
        can't build that table for you here). Instead, this looks for one
        of the following, in order:
 
          1. `mode["eetc_sequence_name"]` -- an explicit EETC sequence
             config_name string, if you already know it.
          2. `mode["eetc_search"]` -- a dict of config_table column/value
             pairs to search for (e.g.
             `{"mode": "excam_imaging", "cfam": "1F", "fpam": "open_12", ...}`),
             mirroring the `search_dict` pattern in the EETC tutorial
             notebook.
 
        Use `search_eetc_sequences()` (below) to explore `config_table` and
        find the right values for your mode, matching the EETC tutorial's
        guidance: for coronagraphic observations, select the "open" FPAM
        setting (e.g. "open_12" for bands 1/2), *not* the coronagraph mask
        setting -- the flux computed here is the *unocculted* reference
        flux, matching the role of `TL.starFlux()` in the rest of
        `Cp_Cb_Csp`. Coronagraph-specific attenuation (contrast, PSF core
        throughput, etc.) is still applied afterwards using cgi_noise's own
        `throughput_rates`/`rawContrast`/`cg` values, same as before.
 
        Args:
            mode (dict): Selected observing mode.
 
        Returns:
            str: EETC sequence `config_name`.
        """
        if "eetc_sequence_name" in mode:
            return mode["eetc_sequence_name"]
 
        if "eetc_search" in mode:
            mask = np.ones(len(config_table["mode"]), dtype=bool)
            for key, val in mode["eetc_search"].items():
                mask &= config_table[key] == val
            matches = config_table[mask]["config_name"].value
            if len(matches) == 0:
                raise ValueError(
                    "No EETC sequence matches mode['eetc_search'] = "
                    f"{mode['eetc_search']!r}. Use search_eetc_sequences() "
                    "to explore available config_table entries."
                )
            return matches[0]
 
        raise ValueError(
            "Cp_Cb_Csp(..., fluxModel='eetc') requires either "
            "mode['eetc_sequence_name'] (an explicit EETC sequence name) "
            "or mode['eetc_search'] (a dict of config_table column/value "
            "pairs to search for, e.g. {'mode': 'excam_imaging', "
            "'cfam': '1F', 'fpam': 'open_12', 'dpam': 'imaging_lens', "
            "'spam_lsam': 'open_spam_nfov_lsam', 'fsam': 'open'}) to be "
            "set. See search_eetc_sequences() to find the right values "
            "for your mode."
        )
 
    def search_eetc_sequences(self, **kwargs):
        """Search EETC's `config_table` (all available sequence
        configurations) for entries matching the given column/value pairs.
 
        Convenience wrapper mirroring the search pattern used in the EETC
        tutorial notebook (`eetc_tutorial_CPP_obs.ipynb`, cell 12). Useful
        for finding the right value to put in
        `mode["eetc_search"]`/`mode["eetc_sequence_name"]`.
 
        Args:
            **kwargs: config_table column/value pairs to filter on, e.g.
                `search_eetc_sequences(mode="excam_imaging", cfam="1F")`.
                With no arguments, returns the full table.
 
        Returns:
            astropy.table.Table: Matching rows of `config_table`.
        """
        mask = np.ones(len(config_table["mode"]), dtype=bool)
        for key, val in kwargs.items():
            mask &= config_table[key] == val
        return config_table[mask]
 
    def _get_eetc_star_flux(self, TL, sInds, mode):
        """Compute per-star unocculted flux rates using EETC's
        `CGIEETC.calc_flux_rate`, as a replacement for `TL.starFlux()`.
 
        For each star, builds (or reuses a cached) `CGIEETC` object from
        the star's V magnitude (`TL.Vmag`) and a best-match spectral type
        (`TL.Spec`, matched against EETC's discrete flux grid via
        `match_spt_to_eetc_grid`), then evaluates `calc_flux_rate()` for the
        sequence resolved by `_get_eetc_sequence_name(mode)`.
 
        Important: the returned flux is in **photoelectrons/second**,
        integrated over the full CFAM bandpass at the focal plane, with the
        full open/unocculted OTA+CGI throughput *and* detector QE already
        applied (see `CGIEETC.calc_flux_rate` docstring) -- this is a very
        different quantity from `TL.starFlux()`'s ph/m^2/s flux density at
        the telescope aperture. `Cp_Cb_Csp` accounts for this by not re-applying the
        collecting area (`Acol`) to EETC-sourced fluxes, and by dividing
        out an estimated QE (see `qe_backout` below) so the rest of the
        Nemati noise pipeline -- which expects pre-QE photon rates -- stays
        self-consistent regardless of which `detectorModel` is chosen.
 
        Args:
            TL (:ref:`TargetList`): TargetList class object.
            sInds (~numpy.ndarray(int)): Integer indices of the stars of
                interest.
            mode (dict): Selected observing mode.
 
        Returns:
            tuple:
                flux_star (~numpy.ndarray(float)):
                    Total unocculted flux rate per star [e-/s].
                flux_star_peak (~numpy.ndarray(float)):
                    Peak-pixel unocculted flux rate per star [e-/s].
                matched_spts (list(str)):
                    The EETC flux-grid spectral type actually used for each
                    star (for diagnostics -- may differ from the catalog's
                    raw spectral type, since EETC only supports a discrete
                    grid; see `match_spt_to_eetc_grid`).
        """
        sequence_name = self._get_eetc_sequence_name(mode)
 
        flux_star = np.zeros(len(sInds))
        flux_star_peak = np.zeros(len(sInds))
        matched_spts = []
 
        for jj, ss in enumerate(sInds):
            mag = float(TL.Vmag[ss])
            raw_spt = str(TL.Spec[ss]).strip()
            spt = match_spt_to_eetc_grid(raw_spt)
            matched_spts.append(spt)
 
            cache_key = (spt, round(mag, 3))
            cgieetc = self._eetc_star_cache.get(cache_key)
            if cgieetc is None:
                cgieetc = CGIEETC(
                    mag=mag, phot="v", spt=spt, pointer_path=pointer_path_sci
                )
                self._eetc_star_cache[cache_key] = cgieetc
 
            total_flux_rate, peak_flux_rate = cgieetc.calc_flux_rate(
                sequence_name=sequence_name
            )
            flux_star[jj] = total_flux_rate
            flux_star_peak[jj] = peak_flux_rate
 
        return flux_star, flux_star_peak, matched_spts

    def Cp_Cb_Csp(self, TL, sInds, fZ, JEZ, dMag, WA, mode, returnExtra=False, TK=None):
        """Calculates electron count rates for planet signal, background noise,
        and speckle residuals.

        Args:
            TL (:ref:`TargetList`):
                TargetList class object
            sInds (~numpy.ndarray(int)):
                Integer indices of the stars of interest
            fZ (~astropy.units.Quantity(~numpy.ndarray(float))):
                Surface brightness of local zodiacal light in units of 1/arcsec2
            JEZ (~astropy.units.Quantity(~numpy.ndarray(float))):
                Intensity of exo-zodiacal light in units of ph/s/m2/arcsec2
            dMag (~numpy.ndarray(float)):
                Differences in magnitude between planets and their host star
            WA (~astropy.units.Quantity(~numpy.ndarray(float))):
                Working angles of the planets of interest in units of arcsec
            mode (dict):
                Selected observing mode
            returnExtra (bool):
                Optional flag, default False, set True to return additional rates for
                validation
            TK (:ref:`TimeKeeping`, optional):
                Optional TimeKeeping object (default None), used to model detector
                degradation effects where applicable.


        Returns:
            tuple:
                C_p (~astropy.units.Quantity(~numpy.ndarray(float))):
                    Planet signal electron count rate in units of 1/s
                C_b (~astropy.units.Quantity(~numpy.ndarray(float))):
                    Background noise electron count rate in units of 1/s
                C_sp (~astropy.units.Quantity(~numpy.ndarray(float))):
                    Residual speckle spatial structure (systematic error)
                    in units of 1/s

        """

        # cast sInds to array
        sInds = np.array(sInds, ndmin=1, copy=copy_if_needed)

        # Star fluxes (ph/m^2/s)
        flux_star = TL.starFlux(sInds, mode).flatten()

        # check if stars identified have vmag 9 or greater, must be before the loop
        vmag = TL.Vmag  # create array of VMag
        vmag_greater_than_9 = vmag > 9
        names_greater_than_9 = TL.Name[vmag_greater_than_9]

        if np.any(vmag_greater_than_9):
            warnings.warn(
                "Integration times for these targets may not be accurate: "
                f"{names_greater_than_9}"
            )

        # get mode elements
        syst = mode["syst"]
        inst = mode["inst"]
        lam_m = mode["lam"].to_value(u.m)
        QE_img = (
            inst["DET_QE_Data"]
            .df.loc[
                inst["DET_QE_Data"].df["lambda_nm"] <= mode["lam"].to_value(u.nm),
                "QE_at_neg100degC",
            ]
            .iloc[-1]
        )

        # set default degredation time if TimeKeeping object not provided
        if TK is None:
            monthsAtL2 = 21
        else:
            monthsAtL2 = TK.currentTimeNorm.to_value(u.d) / 30.4375  # convert to months

        if monthsAtL2 > 63:
            warnings.warn(
                f"You have specified a time at L2 of {monthsAtL2} months.  "
                "The detector degradation model is not valid beyond 63 "
                "months, and may produce anomolous values."
            )

        # allocate outputs
        C_p = np.zeros(len(sInds))
        C_b = np.zeros(len(sInds))
        C_sp = np.zeros(len(sInds))
        if returnExtra:
            extra = {
                "dQE": np.zeros(len(sInds)),
                "mpix": np.zeros(len(sInds)),
                "throughput_rates": np.zeros(len(sInds), dtype=object),
                "cphrate": np.zeros(len(sInds), dtype=object),
                "ENF": np.zeros(len(sInds)),
                "effReadnoise": np.zeros(len(sInds)),
                "frameTime": np.zeros(len(sInds)),
                "QE_img": np.zeros(len(sInds)),
                "nvRatesCore": np.zeros(len(sInds), dtype=object),
                "detNoiseRate": np.zeros(len(sInds), dtype=object),
                "photonCounting": np.zeros(len(sInds), dtype=bool),
            }

        # loop through all values
        for jj, ss in enumerate(sInds):
            if WA.size == 1:
                planetWA = WA[0]
            else:
                planetWA = WA[jj]

            # check for out of bounds WA
            if (planetWA < mode["IWA"]) or (planetWA > mode["OWA"]):
                C_p[jj] = 0
                C_b[jj] = 0
                C_sp[jj] = 0
                continue

            if isinstance(dMag, (int, float)):
                dMagi = dMag
            elif len(dMag) == 1:
                dMagi = dMag[0]
            else:
                dMagi = dMag[jj]

            if len(fZ) == 1:
                fZi = fZ[0]
            else:
                fZi = fZ[jj]

            if len(JEZ) == 1:
                JEZi = JEZ[0]
            else:
                JEZi = JEZ[jj]

            # package up coronagraph values
            cg = self.construct_cg(mode, planetWA)

            # grab relevant detector values
            if "mpix" in mode:
                mpix = mode["mpix"]
            else:
                mpix = (
                    cg.omegaPSF
                    * self.radas**2
                    * (lam_m / cg.CGdesignWL) ** 2
                    * (2 * self.pupilDiam / inst["CritLam"]).decompose().value ** 2
                )

            # get throughput values
            _, throughput_rates = fl.compute_throughputs(
                syst["Throughput_Data"], cg, "uniform"
            )

            # get contrast stability values (all are ppb in the interpolants)
            rawContrast = (
                syst["AvgRawContrast"](mode["lam"], planetWA)[0]
                * 1e-9
                * mode["contrast_degradation"]
            )
            if "SystematicC" in syst:
                SystematicCont = syst["SystematicC"](mode["lam"], planetWA)[0] * 1e-9
            else:
                SystematicCont = 0
            ExtContStab = syst["ExtContStab"](mode["lam"], planetWA)[0] * 1e-9
            IntContStab = syst["IntContStab"](mode["lam"], planetWA)[0] * 1e-9
            selDeltaC = np.sqrt(
                (ExtContStab**2) + (IntContStab**2) + (SystematicCont**2)
            )

            # get count rates for star, planet, speckle
            starFlux = flux_star[jj].value
            planetFlux = starFlux * 10.0 ** (-0.4 * dMagi)
            Acol = self.pupilArea.to_value(u.m**2)
            planet_rate = planetFlux * throughput_rates["planet"] * Acol
            speckle_rate = (
                starFlux
                * rawContrast
                * cg.PSFpeakI
                * cg.CGintmpix
                * throughput_rates["speckle"]
                * Acol
            )

            # get zodi rates
            F0_ph_s_m2 = mode["F0"].to_value(u.ph / u.s / u.m**2)
            locZodiAngFlux = F0_ph_s_m2 * fZi.to_value(1 / u.arcsec**2)
            exoZodiAngFlux = JEZi.to_value(u.ph / u.s / u.m**2 / u.arcsec**2)
            locZodi = (
                locZodiAngFlux * cg.omegaPSF * throughput_rates["local_zodi"] * Acol
            )
            exoZodi = exoZodiAngFlux * cg.omegaPSF * throughput_rates["exo_zodi"] * Acol

            cphrate = corePhotonRates(
                planet=planet_rate,
                speckle=speckle_rate,
                locZodi=locZodi,
                exoZodi=exoZodi,
                straylt=mode["stray_ph_s_pix"] * mpix,
            )
            cphrate.total = sum(
                [
                    cphrate.planet,
                    cphrate.speckle,
                    cphrate.locZodi,
                    cphrate.exoZodi,
                    cphrate.straylt,
                ]
            )

            # pre-compute frame time
            if self.forcePhotonCounting:
                photonCounting = True
            else:
                frameTime = round(
                    min(
                        self.tfmax,
                        max(
                            self.tfmin,
                            self.desiredRate / (cphrate.total * QE_img / mpix),
                        ),
                    ),
                    1,
                )
                approxPerPixelPerFrame = frameTime * cphrate.total * QE_img / mpix
                if approxPerPixelPerFrame <= self.frameThresh:
                    photonCounting = True
                else:
                    photonCounting = False

            ENF, effReadnoise, frameTime, dQE, QE_img = fl.compute_frame_time_and_dqe(
                self.desiredRate,
                self.tfmin,
                self.tfmax,
                photonCounting,
                inst["DET_QE_Data"],
                inst["DET_CBE_Data"],
                lam_m,
                mpix,
                cphrate.total,
            )

            detNoiseRate = fl.detector_noise_rates(
                inst["DET_CBE_Data"], monthsAtL2, frameTime, mpix, True
            )

            rdi_penalty = fl.rdi_noise_penalty(
                mode["inBandFlux0_sum"],
                starFlux,
                mode["TimeonRefStar_tRef_per_tTar"],
                mode["RefStar_SpectralType"],
                mode["RefStar_V_mag"],
            )
            k_sp = rdi_penalty["k_sp"]
            k_det = rdi_penalty["k_det"]
            k_lzo = rdi_penalty["k_lzo"]
            k_ezo = rdi_penalty["k_ezo"]

            nvRatesCore, residSpecSdevRate = fl.noiseRates(
                cphrate,
                QE_img,
                dQE,
                ENF,
                detNoiseRate,
                k_sp,
                k_det,
                k_lzo,
                k_ezo,
                mode["f_SR"],
                starFlux,
                selDeltaC,
                mode["pp_Factor_CBE"],
                cg,
                throughput_rates["speckle"],
                Acol,
            )

            # check for pol mode
            if ("polfraction" in mode) and not (np.isnan(mode["polfraction"])):
                assert (
                    0 <= mode["polfraction"] <= 1
                ), "Polarization fraction must be in [0,1]"

                # if we're doing a pol calculation, need to double detector noise rates
                nvRatesCore.detDark *= 2
                nvRatesCore.detCIC *= 2
                nvRatesCore.detRead *= 2

            # populate outputs
            C_p[jj] = mode["f_SR"] * cphrate.planet * dQE
            C_b[jj] = nvRatesCore.total
            C_sp[jj] = residSpecSdevRate

            if returnExtra:
                extra["dQE"][jj] = dQE
                extra["frameTime"][jj] = frameTime
                extra["mpix"][jj] = mpix
                extra["throughput_rates"][jj] = throughput_rates
                extra["cphrate"][jj] = cphrate
                extra["ENF"][jj] = ENF
                extra["effReadnoise"][jj] = effReadnoise
                extra["QE_img"][jj] = QE_img
                extra["nvRatesCore"][jj] = nvRatesCore
                extra["detNoiseRate"][jj] = detNoiseRate
                extra["photonCounting"][jj] = photonCounting

            # end loop through values
        if returnExtra:
            return C_p << self.inv_s, C_b << self.inv_s, C_sp << self.inv_s, extra

        return C_p << self.inv_s, C_b << self.inv_s, C_sp << self.inv_s

    def Cp_Cb_Csp_eetc(
        self,
        TL,
        sInds,
        fZ,
        JEZ,
        dMag,
        WA,
        mode,
        returnExtra=False,
        TK=None,
        detectorModel="eetc",
        fluxModel="eetc",
    ):
        """Calculates electron count rates for planet signal, background noise,
        and speckle residuals, sourcing star flux and/or the detector model
        from EETC instead of cgi_noise.
 
        This is a standalone counterpart to `Cp_Cb_Csp` -- that function is
        unmodified and always uses cgi_noise for everything. This function
        defaults to EETC for both the flux and detector model (see
        `detectorModel`/`fluxModel` below), but either can be set back to
        `"cgi_noise"` for side-by-side comparisons.
 
        Args:
            TL (:ref:`TargetList`):
                TargetList class object
            sInds (~numpy.ndarray(int)):
                Integer indices of the stars of interest
            fZ (~astropy.units.Quantity(~numpy.ndarray(float))):
                Surface brightness of local zodiacal light in units of 1/arcsec2
            JEZ (~astropy.units.Quantity(~numpy.ndarray(float))):
                Intensity of exo-zodiacal light in units of ph/s/m2/arcsec2
            dMag (~numpy.ndarray(float)):
                Differences in magnitude between planets and their host star
            WA (~astropy.units.Quantity(~numpy.ndarray(float))):
                Working angles of the planets of interest in units of arcsec
            mode (dict):
                Selected observing mode
            returnExtra (bool):
                Optional flag, default False, set True to return additional rates for
                validation
            TK (:ref:`TimeKeeping`, optional):
                Optional TimeKeeping object (default None), used to model detector
                degradation effects where applicable.
            detectorModel (str):
                Either "eetc" (default here) or "cgi_noise". When "eetc",
                the frame-time/ENF/QE and detector-noise-rate calculations
                are sourced from EETC's excam_config.yaml instead of
                cgi_noise's DET_QE_Data/DET_CBE_Data curves.
            fluxModel (str):
                Either "eetc" (default here, uses
                `CGIEETC.calc_flux_rate()`) or "cgi_noise" (uses
                `TL.starFlux()`, matching `Cp_Cb_Csp`). When "eetc", `mode`
                must contain `"eetc_sequence_name"` or `"eetc_search"` (see
                `_get_eetc_sequence_name`).
 
 
        Returns:
            tuple:
                C_p (~astropy.units.Quantity(~numpy.ndarray(float))):
                    Planet signal electron count rate in units of 1/s
                C_b (~astropy.units.Quantity(~numpy.ndarray(float))):
                    Background noise electron count rate in units of 1/s
                C_sp (~astropy.units.Quantity(~numpy.ndarray(float))):
                    Residual speckle spatial structure (systematic error)
                    in units of 1/s
 
        """
        if detectorModel not in ("cgi_noise", "eetc"):
            raise ValueError(
                f"detectorModel must be 'cgi_noise' or 'eetc', got {detectorModel!r}"
            )
        if fluxModel not in ("cgi_noise", "eetc"):
            raise ValueError(
                f"fluxModel must be 'cgi_noise' or 'eetc', got {fluxModel!r}"
            )
 
        # cast sInds to array
        sInds = np.array(sInds, ndmin=1, copy=copy_if_needed)
 
        # check if stars identified have vmag 9 or greater, must be before the loop
        vmag = TL.Vmag  # create array of VMag
        vmag_greater_than_9 = vmag > 9
        names_greater_than_9 = TL.Name[vmag_greater_than_9]
 
        if np.any(vmag_greater_than_9):
            warnings.warn(
                "Integration times for these targets may not be accurate: "
                f"{names_greater_than_9}"
            )
 
        # get mode elements
        syst = mode["syst"]
        inst = mode["inst"]
        lam_m = mode["lam"].to_value(u.m)
 
        detQE_estimate = (
            inst["DET_QE_Data"]
            .df.loc[
                inst["DET_QE_Data"].df["lambda_nm"] <= mode["lam"].to_value(u.nm),
                "QE_at_neg100degC",
            ]
            .iloc[-1]
        )

        if detectorModel == "eetc":
            QE_img = 1.0
        else:
            QE_img = detQE_estimate
 
        if fluxModel == "eetc":
            (
                flux_star_eetc,
                flux_star_peak_eetc,
                eetc_matched_spts,
            ) = self._get_eetc_star_flux(TL, sInds, mode)
        else:
            flux_star = TL.starFlux(sInds, mode).flatten()
 
        # set default degredation time if TimeKeeping object not provided
        if TK is None:
            monthsAtL2 = 21
        else:
            monthsAtL2 = TK.currentTimeNorm.to_value(u.d) / 30.4375  # convert to months
 
        if monthsAtL2 > 63:
            warnings.warn(
                f"You have specified a time at L2 of {monthsAtL2} months.  "
                "The detector degradation model is not valid beyond 63 "
                "months, and may produce anomolous values."
            )
 
        # allocate outputs
        C_p = np.zeros(len(sInds))
        C_b = np.zeros(len(sInds))
        C_sp = np.zeros(len(sInds))
        if returnExtra:
            extra = {
                "dQE": np.zeros(len(sInds)),
                "mpix": np.zeros(len(sInds)),
                "throughput_rates": np.zeros(len(sInds), dtype=object),
                "cphrate": np.zeros(len(sInds), dtype=object),
                "ENF": np.zeros(len(sInds)),
                "effReadnoise": np.zeros(len(sInds)),
                "frameTime": np.zeros(len(sInds)),
                "QE_img": np.zeros(len(sInds)),
                "nvRatesCore": np.zeros(len(sInds), dtype=object),
                "detNoiseRate": np.zeros(len(sInds), dtype=object),
                "photonCounting": np.zeros(len(sInds), dtype=bool),
            }
            if fluxModel == "eetc":
                extra["eetc_matched_spt"] = eetc_matched_spts
                extra["eetc_flux_star"] = flux_star_eetc
                extra["eetc_flux_star_peak"] = flux_star_peak_eetc
 
        # loop through all values
        for jj, ss in enumerate(sInds):
            if WA.size == 1:
                planetWA = WA[0]
            else:
                planetWA = WA[jj]
 
            # check for out of bounds WA
            if (planetWA < mode["IWA"]) or (planetWA > mode["OWA"]):
                C_p[jj] = 0
                C_b[jj] = 0
                C_sp[jj] = 0
                continue
 
            if isinstance(dMag, (int, float)):
                dMagi = dMag
            elif len(dMag) == 1:
                dMagi = dMag[0]
            else:
                dMagi = dMag[jj]
 
            if len(fZ) == 1:
                fZi = fZ[0]
            else:
                fZi = fZ[jj]
 
            if len(JEZ) == 1:
                JEZi = JEZ[0]
            else:
                JEZi = JEZ[jj]
 
            # package up coronagraph values
            cg = self.construct_cg(mode, planetWA)
 
            # grab relevant detector values
            if "mpix" in mode:
                mpix = mode["mpix"]
            else:
                mpix = (
                    cg.omegaPSF
                    * self.radas**2
                    * (lam_m / cg.CGdesignWL) ** 2
                    * (2 * self.pupilDiam / inst["CritLam"]).decompose().value ** 2
                )
 
            # get throughput values
            _, throughput_rates = fl.compute_throughputs(
                syst["Throughput_Data"], cg, "uniform"
            )
 
            # get contrast stability values (all are ppb in the interpolants)
            rawContrast = (
                syst["AvgRawContrast"](mode["lam"], planetWA)[0]
                * 1e-9
                * mode["contrast_degradation"]
            )
            if "SystematicC" in syst:
                SystematicCont = syst["SystematicC"](mode["lam"], planetWA)[0] * 1e-9
            else:
                SystematicCont = 0
            ExtContStab = syst["ExtContStab"](mode["lam"], planetWA)[0] * 1e-9
            IntContStab = syst["IntContStab"](mode["lam"], planetWA)[0] * 1e-9
            selDeltaC = np.sqrt(
                (ExtContStab**2) + (IntContStab**2) + (SystematicCont**2)
            )
 
            # get count rates for star, planet, speckle
            Acol = self.pupilArea.to_value(u.m**2)
            if fluxModel == "eetc":
                starFlux = flux_star_eetc[jj] 
                planetFlux = starFlux * 10.0 ** (-0.4 * dMagi)
                planet_rate = planetFlux * throughput_rates["planet"] * Acol
                speckle_rate = (
                    starFlux
                    * rawContrast
                    * cg.PSFpeakI
                    * cg.CGintmpix
                    * throughput_rates["speckle"]
                    * Acol
                )
            else:
                starFlux = flux_star[jj].value
                planetFlux = starFlux * 10.0 ** (-0.4 * dMagi)
                planet_rate = planetFlux * throughput_rates["planet"] * Acol
                speckle_rate = (
                    starFlux
                    * rawContrast
                    * cg.PSFpeakI
                    * cg.CGintmpix
                    * throughput_rates["speckle"]
                    * Acol
                )
 
            # get zodi rates
            F0_ph_s_m2 = mode["F0"].to_value(u.ph / u.s / u.m**2)
            locZodiAngFlux = F0_ph_s_m2 * fZi.to_value(1 / u.arcsec**2)
            exoZodiAngFlux = JEZi.to_value(u.ph / u.s / u.m**2 / u.arcsec**2)
            locZodi = (
                locZodiAngFlux * cg.omegaPSF * throughput_rates["local_zodi"] * Acol
            )
            exoZodi = exoZodiAngFlux * cg.omegaPSF * throughput_rates["exo_zodi"] * Acol
 
            cphrate = corePhotonRates(
                planet=planet_rate,
                speckle=speckle_rate,
                locZodi=locZodi,
                exoZodi=exoZodi,
                straylt=mode["stray_ph_s_pix"] * mpix,
            )
            cphrate.total = sum(
                [
                    cphrate.planet,
                    cphrate.speckle,
                    cphrate.locZodi,
                    cphrate.exoZodi,
                    cphrate.straylt,
                ]
            )
 
            # pre-compute frame time
            if self.forcePhotonCounting:
                photonCounting = True
            else:
                frameTime = round(
                    min(
                        self.tfmax,
                        max(
                            self.tfmin,
                            self.desiredRate / (cphrate.total * QE_img / mpix),
                        ),
                    ),
                    1,
                )
                approxPerPixelPerFrame = frameTime * cphrate.total * QE_img / mpix
                if approxPerPixelPerFrame <= self.frameThresh:
                    photonCounting = True
                else:
                    photonCounting = False
 
            if detectorModel == "eetc":
                (
                    ENF,
                    effReadnoise,
                    frameTime,
                    dQE,
                    QE_img,
                    g_used,
                ) = self.compute_frame_time_and_dqe_eetc(
                    self.desiredRate,
                    self.tfmin,
                    self.tfmax,
                    photonCounting,
                    cphrate.total,
                )
 
                detNoiseRate = self.detector_noise_rates_eetc(
                    frameTime, mpix, g_used
                )
            else:
                (
                    ENF,
                    effReadnoise,
                    frameTime,
                    dQE,
                    QE_img,
                ) = fl.compute_frame_time_and_dqe(
                    self.desiredRate,
                    self.tfmin,
                    self.tfmax,
                    photonCounting,
                    inst["DET_QE_Data"],
                    inst["DET_CBE_Data"],
                    lam_m,
                    mpix,
                    cphrate.total,
                )
 
                detNoiseRate = fl.detector_noise_rates(
                    inst["DET_CBE_Data"], monthsAtL2, frameTime, mpix, True
                )
 
            rdi_penalty = fl.rdi_noise_penalty(
                mode["inBandFlux0_sum"],
                starFlux,
                mode["TimeonRefStar_tRef_per_tTar"],
                mode["RefStar_SpectralType"],
                mode["RefStar_V_mag"],
            )
            k_sp = rdi_penalty["k_sp"]
            k_det = rdi_penalty["k_det"]
            k_lzo = rdi_penalty["k_lzo"]
            k_ezo = rdi_penalty["k_ezo"]
 
            nvRatesCore, residSpecSdevRate = fl.noiseRates(
                cphrate,
                QE_img,
                dQE,
                ENF,
                detNoiseRate,
                k_sp,
                k_det,
                k_lzo,
                k_ezo,
                mode["f_SR"],
                starFlux,
                selDeltaC,
                mode["pp_Factor_CBE"],
                cg,
                throughput_rates["speckle"],
                Acol,
            )
 
            # check for pol mode
            if ("polfraction" in mode) and not (np.isnan(mode["polfraction"])):
                assert (
                    0 <= mode["polfraction"] <= 1
                ), "Polarization fraction must be in [0,1]"
 
                # if we're doing a pol calculation, need to double detector noise rates
                nvRatesCore.detDark *= 2
                nvRatesCore.detCIC *= 2
                nvRatesCore.detRead *= 2
 
            # populate outputs
            C_p[jj] = mode["f_SR"] * cphrate.planet * dQE
            C_b[jj] = nvRatesCore.total
            C_sp[jj] = residSpecSdevRate
 
            if returnExtra:
                extra["dQE"][jj] = dQE
                extra["frameTime"][jj] = frameTime
                extra["mpix"][jj] = mpix
                extra["throughput_rates"][jj] = throughput_rates
                extra["cphrate"][jj] = cphrate
                extra["ENF"][jj] = ENF
                extra["effReadnoise"][jj] = effReadnoise
                extra["QE_img"][jj] = QE_img
                extra["nvRatesCore"][jj] = nvRatesCore
                extra["detNoiseRate"][jj] = detNoiseRate
                extra["photonCounting"][jj] = photonCounting
 
            # end loop through values
        if returnExtra:
            return C_p << self.inv_s, C_b << self.inv_s, C_sp << self.inv_s, extra
 
        return C_p << self.inv_s, C_b << self.inv_s, C_sp << self.inv_s


    def compute_frame_time_and_dqe_eetc(
        self, desiredRate, tfmin, tfmax, photonCounting, countRateTotal
    ):
        """EETC-sourced replacement for cgi_noise's
        cginoiselib.compute_frame_time_and_dqe().
 
        Computes the per-frame integration time and the EM-gain excess
        noise factor (ENF) using the EXCAM detector parameters loaded from
        EETC's excam_config.yaml (`self.eetc_excam_config`), instead of
        cgi_noise's DET_QE_Data/DET_CBE_Data curves.
 
        Args:
            desiredRate (float):
                Target total electron count per pixel per frame.
            tfmin (float):
                Caller-specified minimum frame time [s] (e.g. self.tfmin).
            tfmax (float):
                Caller-specified maximum frame time [s] (e.g. self.tfmax).
            photonCounting (bool):
                Whether this frame is being read out in photon-counting mode.
            countRateTotal (float):
                Total incident count rate on the pixel of interest
                [e-/pix/s] (cphrate.total).
 
        Returns:
            tuple:
                ENF (float):
                    EM-gain excess noise factor, from eetc.excam_tools._ENF.
                effReadnoise (float):
                    Read noise [e-/pix], referred to the EM gain used here.
                frameTime (float):
                    Selected frame time [s].
                dQE (float):
                    Set to 1.0 -- EETC's flux calibration already includes
                    QE and excam_config.yaml has no separate QE curve.
                QE_img (float):
                    Set to 1.0, for the same reason as dQE.
                g (float):
                    EM gain used to compute ENF/effReadnoise (1.0 for analog
                    frames, gconst or gmax for photon-counting frames).
        """
        (
            darke,
            cic,
            rn,
            X,
            a,
            Lij,
            alpha0,
            fwc,
            alpha1,
            fwc_em,
            Nmin,
            Nmax,
            tmin,
            tmax,
            gmax,
            gconst,
            n,
            Nem,
            tol,
            delta_constr,
            overhead,
            pc_ecount_max,
            T_factor,
        ) = _unpack_excam_config(self.eetc_excam_config)
 
        # combine EETC-config bounds with the caller-specified bounds
        t_lo = max(tfmin, tmin)
        t_hi = min(tfmax, tmax)
 
        frameTime = float(np.clip(desiredRate / countRateTotal, t_lo, t_hi))
 
        if photonCounting:
            g = gconst if gconst is not None else gmax
        else:
            g = 1.0
 
        ENF = _ENF(g, Nem)
        effReadnoise = rn / g
 
        dQE = 1.0
        QE_img = 1.0
 
        return ENF, effReadnoise, frameTime, dQE, QE_img, g
 
    def detector_noise_rates_eetc(self, frameTime, mpix, g):
        """EETC-sourced replacement for cgi_noise's
        cginoiselib.detector_noise_rates().
 
        Uses EXCAM dark current and CIC values from EETC's
        excam_config.yaml (`self.eetc_excam_config`). Unlike cgi_noise's
        detector_noise_rates(), no mission-elapsed-time ("monthsAtL2")
        radiation-damage degradation is applied here, because
        excam_config.yaml exposes single darke/cic values with no
        analogous time dependence at this interface
 
        Args:
            frameTime (float):
                Frame time [s].
            mpix (float):
                Number of pixels contributing to the photometric aperture.
            g (float):
                EM gain (used to refer read noise, matching
                compute_frame_time_and_dqe_eetc's convention).
 
        Returns:
            SimpleNamespace:
                Object exposing .detDark, .detCIC, .detRead, and .total
                attributes (electron rate contributions, [e-/pix/s]), to
                match the attribute interface cgi_noise.cginoiselib's
                noiseRates() expects from a detNoiseRate object.
        """
        (
            darke,
            cic,
            rn,
            X,
            a,
            Lij,
            alpha0,
            fwc,
            alpha1,
            fwc_em,
            Nmin,
            Nmax,
            tmin,
            tmax,
            gmax,
            gconst,
            n,
            Nem,
            tol,
            delta_constr,
            overhead,
            pc_ecount_max,
            T_factor,
        ) = _unpack_excam_config(self.eetc_excam_config)
 
        dark = darke * mpix
        CIC = cic * mpix / frameTime
        read = (rn / g) ** 2 * mpix / frameTime
        return SimpleNamespace(dark=dark, CIC=CIC, read=read)
        


    def calc_polfrac(self, p_in, _C_p, mode):
        """
        Compute measured polarization fraction p_f for a given intrinsic
        polarization fraction p_in.
        """
        theta = mode["theta"]

        Pol0 = _C_p * 96.2
        Pol45 = _C_p * 96.5

        I0 = Pol0 / 2 * (1 + p_in * np.cos(2 * theta))
        I90 = Pol0 / 2 * (1 - p_in * np.cos(2 * theta))
        I45 = Pol45 / 2 * (1 + p_in * np.sin(2 * theta))
        I135 = Pol45 / 2 * (1 - p_in * np.sin(2 * theta))

        I_in = (I0 + I90 + I45 + I135) / 2
        Q_in = I0 - I90
        U_in = I45 - I135

        mat = np.array([I_in, Q_in, U_in, [0.0]])

        # Instrument Mueller matrix
        in_mat = mode["inst"]["matrix"] @ mat

        I_m = in_mat[0]
        Q_m = -in_mat[1]
        U_m = in_mat[2]

        # Measured polarization fraction
        with np.errstate(divide="ignore", invalid="ignore"):
            p_f = np.sqrt(Q_m**2 + U_m**2) / I_m

        return p_f

    def calc_intTime(self, TL, sInds, fZ, JEZ, dMag, WA, mode, TK=None):
        """Finds integration times of target systems for a specific observing
        mode (imaging or characterization), based on Nemati 2014 (SPIE).

        Args:
            TL (TargetList module):
                TargetList class object
            sInds (integer ndarray):
                Integer indices of the stars of interest
            fZ (astropy Quantity array):
                Surface brightness of local zodiacal light in units of 1/arcsec2
            JEZ (astropy Quantity array):
                Intensity of exo-zodiacal light in units of ph/s/m2/arcsec2
            dMag (float ndarray):
                Differences in magnitude between planets and their host star
            WA (astropy Quantity array):
                Working angles of the planets of interest in units of arcsec
            mode (dict):
                Selected observing mode
            TK (TimeKeeping object):
                Optional TimeKeeping object (default None), used to model detector
                degradation effects where applicable.

        Returns:
            intTime (astropy Quantity array):
                Integration times in units of day

        """

        # electron counts
        C_p, C_b, C_sp = self.Cp_Cb_Csp(TL, sInds, fZ, JEZ, dMag, WA, mode, TK=TK)
        _C_p = C_p.to_value(self.inv_s)
        _C_b = C_b.to_value(self.inv_s)
        _C_sp = C_sp.to_value(self.inv_s)

        # get SNR threshold
        SNR = mode["SNR"]
        # calculate integration time based on Nemati 2014
        # if doing a pol calculation, include polarization fraction
        with np.errstate(divide="ignore", invalid="ignore"):
            if ("polfraction" in mode) and not np.isnan(mode["polfraction"]):
                if ("theta" not in mode) or np.isnan(mode["theta"]):
                    mode["theta"] = 0

                p_in = mode["polfraction"]
                p_f = self.calc_polfrac(p_in, _C_p, mode)
                # theta_f = 0.5*(np.arctan(U_m/Q_m)*180/np.pi)

                if ("Cp_ab" not in mode) or np.isnan(mode["Cp_ab"]):
                    mode["Cp_ab"] = 0

                intTime = (
                    np.true_divide(
                        SNR**2.0 * _C_b,
                        (
                            (_C_p * p_f) ** 2.0
                            - (SNR**2 * (_C_sp**2 + mode["Cp_ab"] ** 2))
                        ),
                    )
                    * self.s2d
                )
            else:
                intTime = (
                    np.true_divide(SNR**2.0 * _C_b, (_C_p**2.0 - (SNR * _C_sp) ** 2.0))
                    * self.s2d
                )
        # infinite and NAN are set to zero
        intTime[np.isinf(intTime) | np.isnan(intTime)] = np.nan
        # negative values are set to zero
        intTime[intTime < 0.0] = np.nan

        return intTime << u.d
    
    def calc_intTime_eetc(
        self,
        TL,
        sInds,
        fZ,
        JEZ,
        dMag,
        WA,
        mode,
        TK=None,
        detectorModel="eetc",
        fluxModel="eetc",
    ):
        """Finds integration times of target systems for a specific observing
        mode (imaging or characterization), based on Nemati 2014 (SPIE),
        sourcing star flux and/or the detector model from EETC instead of
        cgi_noise.
 
        This is a standalone counterpart to `calc_intTime` -- that function
        is unmodified and always uses cgi_noise (via `Cp_Cb_Csp`) for
        everything. This function calls `Cp_Cb_Csp_eetc` instead, which
        defaults to EETC for both the flux and detector model (see
        `detectorModel`/`fluxModel` below), but either can be set back to
        `"cgi_noise"` for side-by-side comparisons.
 
        Args:
            TL (TargetList module):
                TargetList class object
            sInds (integer ndarray):
                Integer indices of the stars of interest
            fZ (astropy Quantity array):
                Surface brightness of local zodiacal light in units of 1/arcsec2
            JEZ (astropy Quantity array):
                Intensity of exo-zodiacal light in units of ph/s/m2/arcsec2
            dMag (float ndarray):
                Differences in magnitude between planets and their host star
            WA (astropy Quantity array):
                Working angles of the planets of interest in units of arcsec
            mode (dict):
                Selected observing mode
            TK (TimeKeeping object):
                Optional TimeKeeping object (default None), used to model detector
                degradation effects where applicable.
            detectorModel (str):
                Either "eetc" (default here) or "cgi_noise". Forwarded to
                Cp_Cb_Csp_eetc()
            fluxModel (str):
                Either "eetc" (default here) or "cgi_noise". Forwarded to
                Cp_Cb_Csp_eetc() 
 
        Returns:
            intTime (astropy Quantity array):
                Integration times in units of day
 
        """
 
        # electron counts
        C_p, C_b, C_sp = self.Cp_Cb_Csp_eetc(
            TL,
            sInds,
            fZ,
            JEZ,
            dMag,
            WA,
            mode,
            TK=TK,
            detectorModel=detectorModel,
            fluxModel=fluxModel,
        )
        _C_p = C_p.to_value(self.inv_s)
        _C_b = C_b.to_value(self.inv_s)
        _C_sp = C_sp.to_value(self.inv_s)
 
        # get SNR threshold
        SNR = mode["SNR"]
        # calculate integration time based on Nemati 2014
        # if doing a pol calculation, include polarization fraction
        with np.errstate(divide="ignore", invalid="ignore"):
            if ("polfraction" in mode) and not np.isnan(mode["polfraction"]):
                if ("theta" not in mode) or np.isnan(mode["theta"]):
                    mode["theta"] = 0
 
                p_in = mode["polfraction"]
                p_f = self.calc_polfrac(p_in, _C_p, mode)
                # theta_f = 0.5*(np.arctan(U_m/Q_m)*180/np.pi)
 
                if ("Cp_ab" not in mode) or np.isnan(mode["Cp_ab"]):
                    mode["Cp_ab"] = 0
 
                intTime = (
                    np.true_divide(
                        SNR**2.0 * _C_b,
                        (
                            (_C_p * p_f) ** 2.0
                            - (SNR**2 * (_C_sp**2 + mode["Cp_ab"] ** 2))
                        ),
                    )
                    * self.s2d
                )
            else:
                intTime = (
                    np.true_divide(SNR**2.0 * _C_b, (_C_p**2.0 - (SNR * _C_sp) ** 2.0))
                    * self.s2d
                )
        # infinite and NAN are set to zero
        intTime[np.isinf(intTime) | np.isnan(intTime)] = np.nan
        # negative values are set to zero
        intTime[intTime < 0.0] = np.nan
 
        return intTime << u.d

    def calc_critical_polfraction(self, TL, sInds, fZ, JEZ, dMag, WA, mode, TK=None):
        """
        Returns the critical measured polarization fraction p_in,crit such that
        integration time transitions from undefined to finite.
        """

        # electron counts
        C_p, C_b, C_sp = self.Cp_Cb_Csp(TL, sInds, fZ, JEZ, dMag, WA, mode, TK=TK)
        _C_p = C_p.to_value(self.inv_s)
        _C_b = C_b.to_value(self.inv_s)
        _C_sp = C_sp.to_value(self.inv_s)

        # get SNR threshold
        SNR = mode["SNR"]
        theta = mode["theta"]

        if ("Cp_ab" not in mode) or np.isnan(mode["Cp_ab"]):
            mode["Cp_ab"] = 0

        # Critical polarization fraction
        p_f_crit = SNR * np.sqrt(_C_sp**2 + mode["Cp_ab"] ** 2) / _C_p

        # Anything >1 is physically impossible
        if np.any((p_f_crit <= 0) | (p_f_crit >= 1)):
            return np.nan

        def f_root(p_in):
            return float(self.calc_polfrac(p_in, _C_p, mode)[0] - p_f_crit[0])

        if f_root(0) >= 0:
            return 0

        if f_root(1) < 0:
            return np.nan

        try:
            sol = root_scalar(
                f_root,
                bracket=[0.0, 1.0],
                method="brentq",
            )
            return sol.root

        except Exception:
            return np.nan

    def int_time_denom_obj(self, dMag, *args):
        """
        Objective function for calc_dMag_per_intTime's calculation of the root
        of the denominator of calc_inTime to determine the upper bound to use
        for minimizing to find the correct dMag. Only necessary for coronagraphs.

        Args:
            dMag (~numpy.ndarray(float)):
                dMag being tested
            *args:
                all the other arguments that calc_intTime needs

        Returns:
            ~astropy.units.Quantity(~numpy.ndarray(float)):
                Denominator of integration time expression
        """
        TL, sInds, fZ, JEZ, WA, mode, TK = args
        C_p, C_b, C_sp = self.Cp_Cb_Csp(TL, sInds, fZ, JEZ, dMag, WA, mode, TK=TK)
        denom = (
            C_p.to_value(self.inv_s) ** 2
            - (mode["SNR"] * C_sp.to_value(self.inv_s)) ** 2
        )
        return denom[0]

    def calc_dMag_per_intTime(
        self,
        intTimes,
        TL,
        sInds,
        fZ,
        JEZ,
        WA,
        mode,
        C_b=None,
        C_sp=None,
        TK=None,
        analytic_only=False,
        singularity_dMags=None,
    ):
        """Finds achievable planet delta magnitude for one integration
        time per star in the input list at one working angle.

        Args:
            intTimes (~astropy.units.Quantity(~numpy.ndarray(float))):
                Integration times in units of day
            TL (:ref:`TargetList`):
                TargetList class object
            sInds (numpy.ndarray(int)):
                Integer indices of the stars of interest
            fZ (~astropy.units.Quantity(~numpy.ndarray(float))):
                Surface brightness of local zodiacal light in units of 1/arcsec2
            JEZ (~astropy.units.Quantity(~numpy.ndarray(float))):
                Intensity of exo-zodiacal light in units of ph/s/m2/arcsec2
            WA (~astropy.units.Quantity(~numpy.ndarray(float))):
                Working angles of the planets of interest in units of arcsec
            mode (dict):
                Selected observing mode
            C_b (~astropy.units.Quantity(~numpy.ndarray(float))):
                Background noise electron count rate in units of 1/s (optional)
            C_sp (~astropy.units.Quantity(~numpy.ndarray(float))):
                Residual speckle spatial structure (systematic error) in units of 1/s
                (optional)
            TK (:ref:`TimeKeeping`, optional):
                Optional TimeKeeping object (default None), used to model detector
                degradation effects where applicable.
            analytic_only (bool):
                If True, return the analytic solution for dMag. Not used by the
                Prototype OpticalSystem.
            singularity_dMags (~numpy.ndarray(float), optional):
                Largest attainable delta Mag values.  If None (default) these are
                computed at runtime.


        Returns:
            numpy.ndarray(float):
                Achievable dMag for given integration time and working angle

        """

        # cast sInds to array
        sInds = np.array(sInds, ndmin=1, copy=copy_if_needed)

        # Return NaNs if user requests analytic_only (not supported here)
        if analytic_only:
            return np.full(len(sInds), np.nan)

        # Initialize result array
        dMags = np.zeros(len(sInds))

        for i, int_time in enumerate(tqdm(intTimes, desc="Computing dMags", delay=2)):
            if int_time == 0:
                warnings.warn(
                    "calc_dMag_per_intTime received intTime=0, returning nan."
                )
                dMags[i] = np.nan
                continue

            if np.isnan(int_time):
                warnings.warn(
                    "calc_dMag_per_intTime receive intTime = Nan, returning nan."
                )
                dMags[i] = np.nan
                continue

            if (WA[i] > mode["OWA"]) or (WA[i] < mode["IWA"]):
                warnings.warn("WA outside [IWA, OWA], returning nan.")
                dMags[i] = np.nan
                continue

            s = [sInds[i]]
            args_denom = (TL, s, fZ[i].ravel(), JEZ[i].ravel(), WA[i].ravel(), mode, TK)
            args_intTime = (*args_denom, int_time.ravel())

            # Find the singularity dMag (limit as intTime approaches infinity)
            if singularity_dMags is None:
                if mode["syst"]["occulter"]:
                    singularity_dMag = np.inf
                else:
                    try:
                        f_a = self.int_time_denom_obj(10, *args_denom)
                        f_b = self.int_time_denom_obj(30, *args_denom)

                        if f_a * f_b > 0:
                            warnings.warn(
                                "No root found in bracket [10, 30], returning nan."
                            )
                            singularity_dMag = np.inf
                            dMags[i] = np.nan
                            continue

                        singularity_res = root_scalar(
                            self.int_time_denom_obj,
                            args=args_denom,
                            bracket=[10, 30],
                            method="brentq",
                        )
                        singularity_dMag = singularity_res.root

                    except Exception as e:
                        warnings.warn(f"Root finding failed: {e}, returning nan.")
                        singularity_dMag = np.inf
                        dMags[i] = np.nan
                        continue
            else:
                singularity_dMag = singularity_dMags[i]

            # If infinite intTime, return singularity value and move on
            if int_time == np.inf:
                dMags[i] = singularity_dMag
                continue

            # Alternatively, we need to minimize time error between predicted and
            # desired intTime

            # First we need to establish bounds
            # Initial upper bound is 1e-6 under the singularity dMag
            # Initial lower bound is 5 magnitudes below that
            bounds = singularity_dMag - np.array([5, 1e-6])
            # compute integration time deltas for initial bounds and ensure that there's
            # a sign flip
            lbtime, ubtime = (
                self.calc_intTime(
                    TL,
                    s * 2,
                    args_denom[2],
                    args_denom[3],
                    bounds,
                    args_denom[4],
                    mode,
                    TK,
                )
                - int_time
            )

            # if the top bound produces a shorter shorter integration time, then we
            # can safely return the saturation dMag
            if np.sign(ubtime) != 1:
                dMags[i] = singularity_dMag
                continue

            # the lower bound should produce an integration time below the requested
            # time. if not, need to lower the bound until we find it. but once we do,
            # that's our bounding box right there
            if np.sign(lbtime) != -1:
                while (np.sign(lbtime) != -1) and (bounds[0] >= -1):
                    bounds[0] -= 1
                    lbtime = (
                        self.calc_intTime(
                            TL,
                            s,
                            args_denom[2],
                            args_denom[3],
                            bounds[0],
                            args_denom[4],
                            mode,
                            TK,
                        )[0]
                        - int_time
                    )

                # if loop terminates without finding solution, we're probably dealing
                # with an inversion in the curve and will need to do a finer search
                if np.sign(lbtime) != -1:
                    dMags[i] = np.nan
                    continue

                bounds = np.array([bounds[0], bounds[0] + 1])
            else:

                # do coarse line search over region
                tmp = np.linspace(
                    bounds[0],
                    bounds[1],
                    10,
                )
                tmp2 = (
                    self.calc_intTime(
                        TL,
                        s * (len(tmp) - 2),
                        args_denom[2],
                        args_denom[3],
                        tmp[1:-1],
                        args_denom[4],
                        mode,
                        TK,
                    )
                    - int_time
                )
                tmp2 = np.hstack((lbtime, tmp2, ubtime))

                # look for sign flip
                dsigns = np.diff(np.sign(tmp2))
                ind = np.where(dsigns == 2)[0][0]
                bounds = tmp[ind : ind + 2]

            dMag_init = np.mean(bounds)

            # run minimization
            dMag_min_res = minimize(
                self.dMag_per_intTime_obj,
                x0=np.array([dMag_init]),
                args=args_intTime,
                bounds=[bounds],
                method="L-BFGS-B",
                tol=1e-6,
            )

            success = dMag_min_res.get("success", False)
            if success:
                dMags[i] = (
                    dMag_min_res["x"][0]
                    if isinstance(dMag_min_res["x"], np.ndarray)
                    else dMag_min_res["x"]
                )
                continue

            # if we're here, we failed.  let's try tightening the bounds and minimizing
            # again
            counter = 0
            while not (success) and (counter < 3):
                tmp = np.linspace(
                    bounds[0],
                    bounds[1],
                    10,
                )
                tmp2 = (
                    self.calc_intTime(
                        TL,
                        s * len(tmp),
                        args_denom[2],
                        args_denom[3],
                        tmp,
                        args_denom[4],
                        mode,
                        TK,
                    )
                    - int_time
                )

                # look for sign flip
                dsigns = np.diff(np.sign(tmp2))
                ind = np.where(dsigns == 2)[0][0]
                bounds = tmp[ind : ind + 2]

                dMag_init = np.mean(bounds)

                dMag_min_res = minimize(
                    self.dMag_per_intTime_obj,
                    x0=np.array([dMag_init]),
                    args=args_intTime,
                    bounds=[bounds],
                    method="L-BFGS-B",
                    tol=1e-6,
                )

                success = dMag_min_res.get("success", False)
                counter += 1

            # if we're still failing at this point, let's just accept an abnormal
            # termination if the residual is ok (1e-4 of integration time)
            if not (success):
                if dMag_min_res.fun / int_time.to_value(u.day) < 1e-4:
                    success = True

            if success:
                dMags[i] = (
                    dMag_min_res["x"][0]
                    if isinstance(dMag_min_res["x"], np.ndarray)
                    else dMag_min_res["x"]
                )
            else:
                dMags[i] = np.nan
            # end main loop

        return dMags
