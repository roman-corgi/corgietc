import os
import warnings
from pathlib import Path

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
        contrast_degradation (float)
            Multiplier for rawcontrast (e.g. 2 represents 50% rawcontrast). Defaults to 1.0
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
        contrast_degradation = 1.0,
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
            "CritLam": u.nm,  # ciritcal wavelength
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

    def populate_observingModes_extra(self):
        """Add specific observing mode keywords"""

        self.allowed_observingMode_kws.append("Scenario")
        self.allowed_observingMode_kws.append("StrayLight_Data")
        self.allowed_observingMode_kws.append("pp_Factor_CBE")
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
            #ensure contrast_degradation is in the mode
            mode["contrast_degradation"] = mode.get(
                "contrast_degradation", self.default_vals_extra2["contrast_degradation"])

    def construct_cg(self, mode, WA):
        "Repackage values at a single WA into CGParameters object"

        syst = mode["syst"]
        lam = mode["lam"]

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
        flux_star = TL.starFlux(sInds, mode)

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
            rawContrast = syst["AvgRawContrast"](mode["lam"], planetWA)[0] * 1e-9 / mode["contrast_degradation"]
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

            # TODO: change to JSON input or computed value or per-target calculation
            TimeonRefStar_tRef_per_tTar = 0.25
            rdi_penalty = fl.rdi_noise_penalty(
                mode["inBandFlux0_sum"],
                starFlux,
                TimeonRefStar_tRef_per_tTar,
                "a0v",
                2.26,
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
