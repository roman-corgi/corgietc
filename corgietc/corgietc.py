import os
import warnings

import astropy.units as u
import numpy as np
from scipy import interpolate
import pandas as pd
from pathlib import Path

from EXOSIMS.OpticalSystem.Nemati import Nemati
from cgi_noise import cginoiselib as fl


class corgietc(Nemati):
    def __init__(
        self,
        CritLam=500,
        compbeamD=0.005,
        fnlFocLen=0.26,
        PSF_x_lamD=0.942,
        PSF_y_lamD=0.45,
        Rlamsq=0.0008549637206953,
        Rlam=-1.5131362317830300,
        Rconst=707.8977209483250000,
        **specs,
    ):

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
        }

        Nemati.__init__(self, **specs)

        # add local defaults to outspec
        for k in self.default_vals_extra2:
            self._outspec[k] = self.default_vals_extra2[k]

    def populate_starlightSuppressionSystems_extra(self):

        # add PSFPeak and contrast stability values
        if "PSFpeak" not in self.allowed_starlightSuppressionSystem_kws:
            self.allowed_starlightSuppressionSystem_kws.append("PSFpeak")
        for param_name in [
            "AvgRawContrast",
            "ExtContStab",
            "IntContStab",
            "SystematicC",
            "InitStatContrast",
        ]:
            self.allowed_starlightSuppressionSystem_kws.append(param_name)

        for nsyst, syst in enumerate(self.starlightSuppressionSystems):
            syst = self.get_coro_param(
                syst,
                "PSFpeak",
                expected_ndim=2,
                expected_first_dim=2,
                min_val=0.0,
            )

            for param_name in [
                "AvgRawContrast",
                "ExtContStab",
                "IntContStab",
                "SystematicC",
                "InitStatContrast",
            ]:
                syst = self.get_coro_param(
                    syst,
                    param_name,
                    expected_ndim=2,
                    expected_first_dim=2,
                    min_val=0.0,
                    interp_kind="nearest",
                    update_WAs=False,
                    fill="extrapolate",
                )

            # ensure that CGintSamp is in the system
            syst["CGintSamp"] = syst.get("CGIintSamp", 0.1)

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

            if "imager" in inst["name"].lower():
                pass
            elif "spec" in inst["name"].lower():
                inst["fnumber"] = (
                    (inst["fnlFocLen"] / inst["compbeamD"]).decompose().value
                )
            else:
                raise Exception("Instrument name must contain IMAGER or SPEC")

            inst["pixelScale"] = (inst["CritLam"] / self.pupilDiam / 2).to(
                u.arcsec, equivalencies=u.dimensionless_angles()
            )

    def populate_observingModes_extra(self):
        """Add specific observing mode keywords"""

        self.allowed_observingMode_kws.append("Scenario")

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
                ResPowatPSF = (
                    mode["inst"]["Rconst"]
                    + mode["inst"]["Rlam"] * mode["lam"]
                    + mode["inst"]["Rlamsq"] * mode["lam"] ** 2
                )
                dpix_dlam = ResPowatPSF * xpixPerCor / mode["lam"]
                xpixPerSpec = dpix_dlam * mode["lam"] / mode["inst"]["Rs"]
                mode["mpix"] = xpixPerSpec * ypixPerCor

            else:
                raise Exception("Instrument name must contain IMAGER or SPEC")

    def construct_cg(self, mode, WA):
        "Repackage values at a single WA into CGParameters object"

        syst = mode["syst"]
        lam = mode["lam"]

        radas = ((1 * u.arcsec).to(u.rad)).value

        omegaPSF = (
            (syst["core_area"](lam, WA) / syst["input_angle_unit_value"] ** 2)
            .decompose()
            .value
        )
        CGintmpix = (
            omegaPSF * radas**2
            / ((syst["CGintSamp"] * syst["lam"] / self.pupilDiam) ** 2)
        ).decompose().value

        PSFpeakI = syst["PSFpeak"](lam, WA)[0]
        if "IMG_NFB1_HLC" in mode["Scenario"]:
            PSFpeakI /= CGintmpix

        CG_PSFarea_sqlamD = omegaPSF / (syst["lam"].to_value(u.m) / radas)**2

        out = fl.CGParameters(
            CGcoreThruput=syst["core_thruput"](lam, WA)[0],
            PSFpeakI=PSFpeakI,
            omegaPSF=omegaPSF,
            CGintSamp=syst["CGintSamp"],
            CGradius_arcsec=None,
            CGdesignWL=lam.to_value(u.m),
            CGintmpix=CGintmpix,
            CG_PSFarea_sqlamD=CG_PSFarea_sqlamD,
            CGintensity=syst["core_mean_intensity"](lam, WA)[0],
            CG_occulter_transmission=syst["occ_trans"](lam, WA),
            CGcontrast=syst["core_contrast"](lam, WA),
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

