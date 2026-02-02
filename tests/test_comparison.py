import unittest
import os
import json
import copy
import math
from pathlib import Path
from dataclasses import dataclass

import numpy as np
import astropy.units as u
import yaml
import pandas as pd
import corgietc

import cgi_noise.cginoiselib as fl
import cgi_noise.unitsConstants as uc
from cgi_noise.loadCSVrow import loadCSVrow

import EXOSIMS.Prototypes.TargetList
import EXOSIMS.Prototypes.TimeKeeping


class TestCGINoiseVsCORGIETC(unittest.TestCase):

    def setUp(self):
        self.scenarios = [
            "OPT_IMG_NFB1_HLC",
            "CON_IMG_NFB1_HLC",
            "OPT_SPEC_NFB3_SPC",
            "CON_SPEC_NFB3_SPC",
            "OPT_IMG_WFB4_SPC",
            "CON_IMG_WFB4_SPC",
        ]

        scriptfile = os.path.join(
            os.environ["CORGIETC_DATA_DIR"],
            "scripts",
            "CGI_Noise.json",
        )
        with open(scriptfile, "r") as f:
            specs = json.loads(f.read())

        self.TL = EXOSIMS.Prototypes.TargetList.TargetList(
            **copy.deepcopy(specs)
        )
        self.OS = self.TL.OpticalSystem
        self.sInds = 0

        self.TK = EXOSIMS.Prototypes.TimeKeeping.TimeKeeping(missionLife=5.25)
        self.TK.allocate_time(21 * 30.4375 * u.d)

        self.DATA_DIR = Path(os.environ["CGI_NOISE_DATA_DIR"])
        self.SCEN_DIR = self.DATA_DIR / "Scenarios"

        self.results = {}

        for scenario in self.scenarios:

            obs_params = {
                "scenario": f"{scenario}.yml",
                "target_params": {
                    "v_mag": 5.0,
                    "dist_pc": 10.0,
                    "specType": "g0v",
                    "phaseAng_deg": 65,
                    "sma_AU": 4.1536,
                    "radius_Rjup": 5.6211,
                    "geomAlb_ag": 0.44765,
                    "exoZodi": 1,
                },
                "snr": 5.0,
            }

            if obs_params["scenario"] in (
                "OPT_IMG_WFB4_SPC.yml",
                "CON_IMG_WFB4_SPC.yml",
            ):
                obs_params["target_params"]["sma_AU"] = 8

            mode = list(
                filter(
                    lambda m: m["Scenario"] == scenario,
                    self.OS.observingModes,
                )
            )[0]

            syst = mode["syst"]
            inst = mode["inst"]

            fZ = np.repeat(self.TL.ZodiacalLight.fZ0, 1)
            JEZ = self.TL.JEZ0[mode["hex"]] / (
                obs_params["target_params"]["sma_AU"] ** 2
            )

            scenario_path = self.SCEN_DIR / obs_params["scenario"]
            with open(scenario_path, "r") as file:
                config = yaml.safe_load(file)

            target_params = obs_params["target_params"]
            SNRdesired = obs_params["snr"]

            # === corePhotonRates ===
            @dataclass
            class corePhotonRates:
                planet: float
                speckle: float
                locZodi: float
                exoZodi: float
                straylt: float
                total: float = 0.0

            ObservationCase = config["DataSpecification"]["ObservationCase"]

            # === instrument params ===
            DPM = config["instrument"]["Diam"]
            lam = config["instrument"]["wavelength"]
            lamD = lam / DPM
            intTimeDutyFactor = config["instrument"]["dutyFactor"]
            opMode = config["instrument"]["OpMode"]
            bandWidth = config["instrument"]["bandwidth"]

            target = fl.Target(**target_params)
            sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
            target.albedo = target.albedo_from_geomAlbedo(
                target.phaseAng_deg, target.geomAlb_ag
            )

            filenameList = fl.getScenFileNames(config, self.DATA_DIR)
            CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = (
                fl.loadCSVs(filenameList)
                )

            CS_Type = config["DataSpecification"]["CS_Type"]

            IWA, OWA = fl.workingAnglePars(CG_Data, CS_Data)

            planetWA = 6.1
            if planetWA < IWA:
                planetWA = IWA

            WA = np.array([planetWA]) * (
                mode["lam"] / self.OS.pupilDiam
            ).to(u.arcsec, equivalencies=u.dimensionless_angles())

            tol = 0.05
            if (IWA - tol) <= planetWA <= IWA:
                planetWA = IWA
            elif OWA <= planetWA <= (OWA + tol):
                planetWA = OWA
            elif planetWA < (IWA - tol) or planetWA > (OWA + tol):
                raise ValueError(
                    f"Planet WA={planetWA:.1f} outside of IWA={IWA:.1f}, OWA={OWA:.1f}."
                )

            selDeltaC, AvgRawC, SystematicC, initStatRaw, IntContStab, ExtContStab = (
                fl.contrastStabilityPars(CS_Type, planetWA, CS_Data)
            )

            cg = fl.coronagraphParameters(
                CG_Data.df, config, planetWA, DPM
            )

            f_SR, CritLam, detPixSize_m, mpix = fl.getFocalPlaneAttributes(
                    opMode,
                    config,
                    DET_CBE_Data,
                    lam,
                    bandWidth,
                    DPM,
                    cg.CGdesignWL,
                    cg.omegaPSF,
                    self.DATA_DIR,
                )

            inBandFlux0_sum, inBandZeroMagFlux, starFlux = fl.getSpectra(
                target, lam, bandWidth, self.DATA_DIR
            )
            TimeonRefStar_tRef_per_tTar = 0.25

            fluxRatio = (
                target.albedo
                * (target.radius_Rjup * uc.jupiterRadius / (target.sma_AU * uc.AU)) ** 2
            )
            planetFlux = fluxRatio * starFlux

            magLocalZodi = config["instrument"]["LocZodi_magas2"]
            magExoZodi_1AU = config["instrument"]["ExoZodi_magas2"]
            absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10)
            locZodiAngFlux = inBandZeroMagFlux * 10 ** (-0.4 * magLocalZodi)
            exoZodiAngFlux = (
                inBandZeroMagFlux
                * 10 ** (-0.4 * (absMag - uc.sunAbsMag + magExoZodi_1AU))
                / target.sma_AU**2
                * target.exoZodi
            )

            thput, throughput_rates = fl.compute_throughputs(THPT_Data, cg, "uniform")
            Acol = (np.pi / 4) * DPM**2
            stray_ph_s_mm2 = fl.getStrayLightfromfile(ObservationCase, "CBE", STRAY_FRN_Data)
            stray_ph_s_pix = stray_ph_s_mm2 * (1 / uc.mm**2) * detPixSize_m**2

            cphrate = corePhotonRates(
                planet=planetFlux * throughput_rates["planet"] * Acol,
                speckle=starFlux
                * AvgRawC
                * cg.PSFpeakI
                * cg.CGintmpix
                * throughput_rates["speckle"]
                * Acol,
                locZodi=locZodiAngFlux * cg.omegaPSF * throughput_rates["local_zodi"] * Acol,
                exoZodi=exoZodiAngFlux * cg.omegaPSF * throughput_rates["exo_zodi"] * Acol,
                straylt=stray_ph_s_pix * mpix,
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
                0.1, 3, 100, True, QE_Data, DET_CBE_Data, lam, mpix, cphrate.total
            )

            detNoiseRate = fl.detector_noise_rates(DET_CBE_Data, 21, frameTime, mpix, True)

            rdi_penalty = fl.rdi_noise_penalty(
                inBandFlux0_sum, starFlux, TimeonRefStar_tRef_per_tTar, "a0v", 2.26
            )
            k_sp = rdi_penalty["k_sp"]
            k_det = rdi_penalty["k_det"]
            k_lzo = rdi_penalty["k_lzo"]
            k_ezo = rdi_penalty["k_ezo"]

            nvRatesCore, residSpecRate = fl.noiseRates(
                cphrate,
                QE_img,
                dQE,
                ENF,
                detNoiseRate,
                k_sp,
                k_det,
                k_lzo,
                k_ezo,
                f_SR,
                starFlux,
                selDeltaC,
                config["instrument"]["pp_Factor_CBE"],
                cg,
                throughput_rates["speckle"],
                Acol,
            )

            planetSignalRate = f_SR * cphrate.planet * dQE
            timeToSNR, criticalSNR = fl.compute_tsnr(
                SNRdesired, planetSignalRate, nvRatesCore, residSpecRate
            )

            csfilename = None
            for filepath in filenameList:
                base = os.path.basename(filepath)
                if base.startswith("CS_"):
                    csfilename = os.path.splitext(base)[0]
                    break


            omegaPSF_corgi = (
                syst["core_area"](mode["lam"], WA)
                / syst["input_angle_unit_value"] ** 2
            ).decompose().value

            CGintmpix = (
                omegaPSF_corgi
                * self.OS.radas**2
                / (
                    (syst["CGintSamp"] * syst["lam"] / self.OS.pupilDiam) ** 2
                )
            ).decompose().value

            PSFpeakI = syst["PSFpeak"](mode["lam"], WA)
            if "IMG_NFB1_HLC" in mode["Scenario"]:
                PSFpeakI /= CGintmpix

            per_occtrans = (
                abs(
                    syst["occ_trans"](mode["lam"], WA)
                    - cg.CG_occulter_transmission
                )
                / cg.CG_occulter_transmission
                * 100
            )

            per_thruput = (
                abs(
                    syst["core_thruput"](mode["lam"], WA)
                    - cg.CGcoreThruput
                )
                / cg.CGcoreThruput
                * 100
            )
            per_intensity = (
                np.abs(syst["core_mean_intensity"](mode["lam"], WA)[0] - cg.CGintensity)
                / (cg.CGintensity)
                * 100
            )
            per_area = (
                np.abs(
                    syst["core_area"](mode["lam"], WA).value
                    - cg.omegaPSF * lamD**2 * 4.255 * 1e10
                )
                / (cg.omegaPSF * lamD**2 * 4.255 * 1e10)
                * 100
            )
            per_contrast = (
                abs(
                    syst["core_contrast"](mode["lam"], WA)
                    - cg.CGcontrast
                )
                / cg.CGcontrast
                * 100
            )
            per_lam = (
                np.abs(syst["lam"].value - cg.CGdesignWL * 1e9) / (cg.CGdesignWL * 1e9) * 100
            )
            per_rawcontrast = (
                np.abs(syst["AvgRawContrast"](mode["lam"], WA)[0] - AvgRawC * 1e9)
                / (AvgRawC * 1e9)
                * 100
            )
            per_initStatRawContrast = (
                np.abs(syst["InitStatContrast"](mode["lam"], WA)[0] - initStatRaw * 1e9)
                / (initStatRaw * 1e9)
                * 100
            )
            per_IntContStab = (
                np.abs(syst["IntContStab"](mode["lam"], WA)[0] - IntContStab * 1e9)
                / (IntContStab * 1e9)
                * 100
            )
            per_ExtContStab = (
                np.abs(syst["ExtContStab"](mode["lam"], WA)[0] - ExtContStab * 1e9)
                / (ExtContStab * 1e9)
                * 100
            )


            self.results[scenario] = {
                "percent_diff": {
                    "occulter_transmission": per_occtrans,
                    "core_thruput": per_thruput,
                    "core_contrast": per_contrast,
                    "core_mean_intensity": per_intensity,
                    "core_area": per_area,
                    "Lam": per_lam,
                    "avg_raw_contrast": per_rawcontrast,
                }
            }

    def tests(self):
        """
        Compare between cgi_noise and corgietc.
        """
        for scenario, res in self.results.items():
            with self.subTest(scenario=scenario):

                self.assertLess(
                    res["percent_diff"]["occulter_transmission"],
                    5.0,
                    msg=f"{scenario}: occulter transmission mismatch",
                )

                self.assertLess(
                    res["percent_diff"]["core_thruput"],
                    5.0,
                    msg=f"{scenario}: core throughput mismatch",
                )

                self.assertLess(
                    res["percent_diff"]["core_contrast"],
                    15.0,
                    msg=f"{scenario}: core contrast mismatch",
                )


if __name__ == "__main__":
    unittest.main()
