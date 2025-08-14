import unittest
import corgietc
import os
import json
import copy
import sys
import numpy as np
from EXOSIMS.Prototypes.TargetList import TargetList
import astropy.units as u
from synphot import units


class test_corgietc(unittest.TestCase):
    """

    Global OpticalSystem tests.
    Applied to all implementations, for overloaded methods only.

    Any implementation-specific methods, or to test specific new
    method functionality, separate tests are needed.

    """

    def setUp(self):

        self.dev_null = open(os.devnull, "w")
        self.script = os.path.join(
            os.environ["CORGIETC_DATA_DIR"], "scripts", "CGI_Noise.json"
        )

        with open(self.script) as f:
            self.specs = json.loads(f.read())

        self.nStars = 100
        self.nPoints = 100

        self.specs["VmagFill"] = np.linspace(0, 9, self.nStars)
        self.specs["ntargs"] = self.nStars

        with RedirectStreams(stdout=self.dev_null):
            self.TL = TargetList(**copy.deepcopy(self.specs))

        # units
        self.PHOTLAM_sr_decomposed_val = (1 * units.PHOTLAM / u.sr).to_value(
            u.ph / u.s / u.arcsec**2 / u.cm**2 / u.nm
        )
        self.F0_unit = u.ph / u.s / u.cm**2
        self.deltaLam_unit = u.nm
        self.inv_arcsec2 = 1 / u.arcsec**2

    def tearDown(self):
        self.dev_null.close()

    def test_intTime_dMag_roundtrip(self):
        """
        Check calc_intTime to calc_dMag_per_intTime to calc_intTime to
        calc_dMag_per_intTime give equivalent results
        """
        TL = self.TL
        OS = self.TL.OpticalSystem

        # set up randomized inputs
        sInds = np.random.choice(np.arange(self.TL.nStars), size=self.nPoints)
        lon = np.random.rand(self.nPoints) * 180 << u.deg  # converts to 0-180 deg
        lat = np.random.rand(self.nPoints) * 90 << u.deg
        Izod = TL.ZodiacalLight.zodi_intensity_at_location(lon, lat)

        for jj, mode in enumerate(OS.observingModes):

            # compute local zodi and exozodi
            Izod_color = Izod * TL.ZodiacalLight.zodi_color_correction_factor(
                mode["lam"]
            )

            factor = (
                units.convert_flux(mode["lam"], Izod_color * u.sr, units.PHOTLAM).value
                / Izod_color.value
            )[0]

            Izod_photons = Izod_color.value * factor * self.PHOTLAM_sr_decomposed_val

            fZ = (
                Izod_photons
                / (
                    mode["F0"].to_value(self.F0_unit)
                    / mode["deltaLam"].to_value(self.deltaLam_unit)
                )
            ) << self.inv_arcsec2

            JEZ = TL.JEZ0[mode["hex"]][sInds]

            # define WAs
            WAs = (
                np.linspace(
                    mode["IWA"].value * 1.01, mode["OWA"].value * 0.99, self.nPoints
                )
                * mode["IWA"].unit
            )

            # compute saturation dMags for mode and define random dMags for test
            sat_dMags = OS.calc_saturation_dMag(TL, sInds, fZ, JEZ, WAs, mode)
            dMags1 = sat_dMags - np.random.rand(TL.nStars) * 5

            # integration times from dMags1
            intTimes1 = OS.calc_intTime(TL, sInds, fZ, JEZ, dMags1, WAs, mode)
            self.assertTrue(
                not np.any(np.isnan(intTimes1)),
                msg=f"NaN integration time computed for {mode["Scenario"]}",
            )

            # dMags from intTime1
            dMags2 = OS.calc_dMag_per_intTime(intTimes1, TL, sInds, fZ, JEZ, WAs, mode)
            self.assertTrue(
                not np.any(np.isnan(dMags2)),
                msg=f"NaN Delta mag computed for {mode["Scenario"]}",
            )

            # integration times from dMags2
            intTimes2 = OS.calc_intTime(TL, sInds, fZ, JEZ, dMags2, WAs, mode)
            self.assertTrue(
                not np.any(np.isnan(intTimes2)),
                msg=f"NaN integration time computed for {mode["Scenario"]}",
            )

            # compute errors
            dMag_err = np.abs(dMags1 - dMags2)
            intTime_err = (np.abs(intTimes1 - intTimes2) / intTimes1).value

            # ensure majority of dMags are within range
            self.assertTrue(np.where(dMag_err > 1e-3)[0].size / self.nPoints < 0.25)

            # ensure intTimes for larger dMag errors match to within 10 seconds
            self.assertTrue(intTime_err[dMag_err > 1e-3]) < 0.05


class RedirectStreams(object):
    r"""Set stdout and stderr to redirect to the named streams.

    Used for eliminating chatter to stdout upon module creation."""

    def __init__(self, stdout=None, stderr=None):
        self._stdout = stdout or sys.stdout
        self._stderr = stderr or sys.stderr

    def __enter__(self):
        self.old_stdout, self.old_stderr = sys.stdout, sys.stderr
        self.old_stdout.flush()
        self.old_stderr.flush()
        sys.stdout, sys.stderr = self._stdout, self._stderr

    def __exit__(self, exc_type, exc_value, traceback):
        self._stdout.flush()
        self._stderr.flush()
        sys.stdout = self.old_stdout
        sys.stderr = self.old_stderr
