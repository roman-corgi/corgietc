import os
import warnings

import astropy.units as u
import numpy as np
from scipy import interpolate
import pandas as pd
from pathlib import Path

from EXOSIMS.OpticalSystem.Nemati import Nemati
from cgi_noise import cginoiselib as fl
from cgi_noise.loadCSVrow import loadCSVrow


class corgietc(Nemati):
    def __init__(
        self,
        CritLam=500,
        compbeamD=0.005,
        fnlFocLen=0.26,
        PSF_x_lamD=0.942,
        PSF_y_lamD=0.45,
        Rlamsq=0.0008549637206953,
        Rlam=-1.51313623178303,
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
        data_dir = Path(os.environ["CGI_NOISE_DATA_DIR"])
        self.FocalPlaneAtt = loadCSVrow(
            data_dir / "instrument" / "CONST_SNR_FPattributes.csv"
        )
        self.AmiciPar = loadCSVrow(
            data_dir / "instrument" / "CONST_Amici_parameters.csv"
        )

        # add local defaults to outspec
        for k in self.default_vals_extra2:
            self._outspec[k] = self.default_vals_extra2[k]

    def populate_starlightSuppressionSystems_extra(self):

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

            dat = syst["PSFpeak"]
            self._outspec["starlightSuppressionSystems"][nsyst]["PSFpeak"] = (
                dat.value if isinstance(dat, u.Quantity) else dat
            )

            for param_name in [
                "AvgRawContrast",
                "ExtContStab",
                "IntContStab",
                "SystematicC",
                "InitStatContrast",
            ]:
                value = self.get_coro_param(
                    syst,
                    param_name,
                    expected_ndim=2,
                    expected_first_dim=2,
                    min_val=0.0,
                    interp_kind="nearest",
                    update_WAs=False,
                    fill="extrapolate",
                )
                self._outspec["starlightSuppressionSystems"][nsyst][param_name] = (
                    value.value if isinstance(value, u.Quantity) else value
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
                if kws[kw] is not None:
                    inst[kw] *= kws[kw]

                if "imager" in inst["name"].lower():
                    inst["f_SR"] = 1
                    inst["pixelScale"] = (inst["CritLam"] / self.pupilDiam / 2).to(
                        u.arcsec, equivalencies=u.dimensionless_angles()
                    )
                elif "spec" in inst["name"].lower():

                else:
                    raise Exception("Instrument name must contain IMAGER or SPEC")

    # def populate_observingModes_extra(self):
    #     """Add Nemati_2019-specific observing mode keywords"""

    #     super().populate_observingModes_extra()
    #     self.allowed_observingMode_kws.append("Scenario")

    #     for nmode, mode in enumerate(self.observingModes):
    #         mode["Scenario"] = mode.get(
    #             "Scenario", self.default_vals_extra2["Scenario"]
    #         )
    #         self._outspec["observingModes"][nmode]["Scenario"] = mode[
    #             "Scenario"
    #         ]
