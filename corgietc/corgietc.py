import os
import warnings

import astropy.units as u
import numpy as np
from pathlib import Path

from EXOSIMS.OpticalSystem.Nemati import Nemati
from cgi_noise import cginoiselib as fl
from cgi_noise.main_core import corePhotonRates
from EXOSIMS.util._numpy_compat import copy_if_needed
import astropy.constants as const


class corgietc(Nemati):
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
        **specs,
    ):

        # useful conversion factorс
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
                if param_name in syst:
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
            "pp_Factor_CBE": None,
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

        # loop through all values
        for jj, ss in enumerate(sInds):
            if WA.size == 1:
                planetWA = WA[0]
            else:
                planetWA = WA[jj]

            if len(dMag) == 1:
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
            rawContrast = syst["AvgRawContrast"](mode["lam"], planetWA)[0] * 1e-9
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

            ENF, effReadnoise, frameTime, dQE, QE_img = fl.compute_frame_time_and_dqe(
                0.1,
                3,
                100,
                True,
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

            nvRatesCore, residSpecRate = fl.noiseVarianceRates(
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
                inst["pp_Factor_CBE"],
                cg,
                throughput_rates["speckle"],
                Acol,
            )

            # populate outputs
            C_p[jj] = mode["f_SR"] * cphrate.planet * dQE
            C_b[jj] = nvRatesCore.total
            C_sp[jj] = residSpecRate

            # end loop through values

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

        Returns:
            numpy.ndarray(float):
                Achievable dMag for given integration time and working angle

        .. warning::

            Temporary. To be Updated.

        """
        # cast sInds to array
        sInds = np.array(sInds, ndmin=1, copy=copy_if_needed)

        if (C_b is None) or (C_sp is None):
            _, C_b, C_sp = self.Cp_Cb_Csp(
                TL, sInds, fZ, JEZ, np.zeros(len(sInds)), WA, mode, TK=TK
            )

        C_p = mode["SNR"] * np.sqrt(C_sp**2 + C_b / intTimes)  # planet count rate
        core_thruput = mode["syst"]["core_thruput"](mode["lam"], WA)
        flux_star = TL.starFlux(sInds, mode)

        dMag = -2.5 * np.log10(C_p / (flux_star * mode["losses"] * core_thruput))

        return dMag.value
