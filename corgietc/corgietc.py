import os
import warnings

import astropy.units as u
import numpy as np
from scipy import interpolate
import pandas as pd

from EXOSIMS.OpticalSystem.Nemati import Nemati
from cgi_noise import cginoiselib as fl


class corgietc(Nemati):
    def __init__(
        self,
        **specs,
    ):
        Nemati.__init__(self, **specs)

    def populate_starlightSuppressionSystems_extra(self):
        super().populate_starlightSuppressionSystems_extra()

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
