import json
import yaml
from pathlib import Path
import os
import cgi_noise
import corgietc
import numpy as np

DATA_DIR = Path(os.environ["CGI_NOISE_DATA_DIR"])
SCEN_DIR = DATA_DIR / "Scenarios"

# establishing the github links of the .yml files
opt_b1 = SCEN_DIR / "OPT_IMG_NFB1_HLC.yml"
opt_b3 = SCEN_DIR / "OPT_SPEC_NFB3_SPC.yml"
opt_b4 = SCEN_DIR / "OPT_IMG_WFB4_SPC.yml"
con_b1 = SCEN_DIR / "CON_IMG_NFB1_HLC.yml"
con_b3 = SCEN_DIR / "CON_SPEC_NFB3_SPC.yml"
con_b4 = SCEN_DIR / "CON_IMG_WFB4_SPC.yml"

# storing the data in the .yml files as dictionaries
with open(opt_b1, "r") as f:
    data_opt_b1 = yaml.safe_load(f)
with open(opt_b3, "r") as f:
    data_opt_b3 = yaml.safe_load(f)
with open(opt_b4, "r") as f:
    data_opt_b4 = yaml.safe_load(f)
with open(con_b1, "r") as f:
    data_con_b1 = yaml.safe_load(f)
with open(con_b3, "r") as f:
    data_con_b3 = yaml.safe_load(f)
with open(con_b4, "r") as f:
    data_con_b4 = yaml.safe_load(f)

yml_scenarios = [
    data_opt_b1,
    data_opt_b4,
    data_opt_b3,
    data_con_b1,
    data_con_b4,
    data_con_b3,
]

##########################################################################################
scriptfile = os.path.join(os.environ["CORGIETC_DATA_DIR"], "scripts", "CGI_Noise.json")
with open(scriptfile, "r") as f:
    json_data = json.loads(f.read())

errors = []
scenario_count = 0

# Relationship between scenario count and the observing case

scenario_conversion = {
    0: "OPT_IMG_NFB1_HLC",
    1: "OPT_IMG_WFB4_SPC",
    2: "OPT_SPEC_NFB3_SPC",
    3: "CON_IMG_NFB1_HLC",
    4: "CON_IMG_WFB4_SPC",
    5: "CON_SPEC_NFB3_SPC",
}


for scenario in yml_scenarios:

    # Coronagraph_Data
    occ_trans_name = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count]["occ_trans"]
        )[-1]
    )[0]
    core_mean_intensity_name = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count][
                "core_mean_intensity"
            ]
        )[-1]
    )[0]
    core_area_name = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count]["core_area"]
        )[-1]
    )[0]
    core_contrast_name = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count]["core_contrast"]
        )[-1]
    )[0]
    if (
        scenario["DataSpecification"]["Coronagraph_Data"] != occ_trans_name
        or scenario["DataSpecification"]["Coronagraph_Data"] != core_mean_intensity_name
        or scenario["DataSpecification"]["Coronagraph_Data"] != core_area_name
        or scenario["DataSpecification"]["Coronagraph_Data"] != core_contrast_name
    ):
        errors.append(
            str(scenario_conversion[scenario_count])
            + " Coronagraph_Data file not matching"
        )

    # QE_Curve_Data
    if "B1" in scenario["DataSpecification"]["ObservationCase"]:
        qe_0 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][0]["QE"])[-1]
        )[0]
        det_qe_0 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][0]["DET_QE_Data"])[-1]
        )[0]
        if (
            qe_0 != scenario["DataSpecification"]["QE_Curve_Data"]
            or det_qe_0 != scenario["DataSpecification"]["QE_Curve_Data"]
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " QE_Curve_Data not matching"
            )
    elif "B4" in scenario["DataSpecification"]["ObservationCase"]:
        qe_1 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][1]["QE"])[-1]
        )[0]
        det_qe_0 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][0]["DET_QE_Data"])[-1]
        )[0]
        if (
            qe_1 != scenario["DataSpecification"]["QE_Curve_Data"]
            or det_qe_0 != scenario["DataSpecification"]["QE_Curve_Data"]
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " QE_Curve_Data not matching"
            )
    else:
        qe_2 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][2]["QE"])[-1]
        )[0]
        det_qe_0 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][0]["DET_QE_Data"])[-1]
        )[0]
        if (
            qe_2 != scenario["DataSpecification"]["QE_Curve_Data"]
            or det_qe_0 != scenario["DataSpecification"]["QE_Curve_Data"]
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " QE_Curve_Data not matching"
            )

    # Detector_Data
    if "B1" in scenario["DataSpecification"]["ObservationCase"]:
        det_data_0 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][0]["DET_CBE_Data"])[-1]
        )[0]
        if det_data_0 != scenario["DataSpecification"]["Detector_Data"]:
            errors.append(
                str(scenario_conversion[scenario_count]) + " Detector_Data not matching"
            )
    elif "B4" in scenario["DataSpecification"]["ObservationCase"]:
        det_data_1 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][1]["DET_CBE_Data"])[-1]
        )[0]
        if det_data_1 != scenario["DataSpecification"]["Detector_Data"]:
            errors.append(
                str(scenario_conversion[scenario_count]) + " Detector_Data not matching"
            )
    else:
        det_data_2 = os.path.splitext(
            os.path.split(json_data["scienceInstruments"][2]["DET_CBE_Data"])[-1]
        )[0]
        if det_data_2 != scenario["DataSpecification"]["Detector_Data"]:
            errors.append(
                str(scenario_conversion[scenario_count]) + " Detector_Data not matching"
            )

    # ContrastStability_Data
    avg_raw = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count]["AvgRawContrast"]
        )[-1]
    )[0]
    ext_stab = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count]["ExtContStab"]
        )[-1]
    )[0]
    int_stab = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count]["IntContStab"]
        )[-1]
    )[0]
    init_stat = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count]["InitStatContrast"]
        )[-1]
    )[0]
    if (
        scenario["DataSpecification"]["ContrastStability_Data"] != avg_raw
        or scenario["DataSpecification"]["ContrastStability_Data"] != ext_stab
        or scenario["DataSpecification"]["ContrastStability_Data"] != int_stab
        or scenario["DataSpecification"]["ContrastStability_Data"] != init_stat
    ):
        errors.append(
            str(scenario_conversion[scenario_count])
            + " ContrastStability_Data file not matching"
        )

    # CS_Type
    if (
        scenario["DataSpecification"]["CS_Type"]
        != json_data["starlightSuppressionSystems"][scenario_count]["csv_names"][
            "AvgRawContrast"
        ][:4]
        or scenario["DataSpecification"]["CS_Type"]
        != json_data["starlightSuppressionSystems"][scenario_count]["csv_names"][
            "ExtContStab"
        ][:4]
        or scenario["DataSpecification"]["CS_Type"]
        != json_data["starlightSuppressionSystems"][scenario_count]["csv_names"][
            "IntContStab"
        ][:4]
        or scenario["DataSpecification"]["CS_Type"]
        != json_data["starlightSuppressionSystems"][scenario_count]["csv_names"][
            "SystematicC"
        ][:4]
        or scenario["DataSpecification"]["CS_Type"]
        != json_data["starlightSuppressionSystems"][scenario_count]["csv_names"][
            "InitStatContrast"
        ][:4]
    ):
        errors.append(
            str(scenario_conversion[scenario_count]) + " CS_Type not matching"
        )

    # Throughput_Data
    throughput_name = os.path.splitext(
        os.path.split(
            json_data["starlightSuppressionSystems"][scenario_count]["Throughput_Data"]
        )[-1]
    )[0]
    if scenario["DataSpecification"]["Throughput_Data"] != throughput_name:
        errors.append(
            str(scenario_conversion[scenario_count])
            + " Throughput Data file not matching"
        )

    # StrayLight_Data
    stray_name = os.path.splitext(
        os.path.split(json_data["observingModes"][scenario_count]["StrayLight_Data"])[
            -1
        ]
    )[0]
    if scenario["DataSpecification"]["StrayLight_Data"] != stray_name:
        errors.append(
            str(scenario_conversion[scenario_count])
            + " StrayLight_Data file not matching"
        )

    # Diameter
    if scenario["instrument"]["Diam"] != json_data["pupilDiam"]:
        errors.append(
            str(scenario_conversion[scenario_count]) + " Diameter not matching."
        )

    # Area
    if scenario["instrument"]["Acol"] != 4.385:
        errors.append(str(scenario_conversion[scenario_count]) + " Area not matching.")

    # Wavelength
    if float(scenario["instrument"]["wavelength"] * 10**9) != float(
        json_data["starlightSuppressionSystems"][scenario_count]["lam"]
    ):
        errors.append(
            str(scenario_conversion[scenario_count]) + " Wavelength not matching."
        )

    # Bandwdith
    if float(scenario["instrument"]["bandwidth"]) != float(
        json_data["starlightSuppressionSystems"][scenario_count]["BW"]
    ):
        errors.append(
            str(scenario_conversion[scenario_count]) + " bandwidth not matching."
        )

    # CGtype
    if scenario_count == 0 or scenario_count == 3:
        if (
            scenario["instrument"]["CGtype"]
            != json_data["starlightSuppressionSystems"][0]["name"][-3:]
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " CGtype not matching."
            )
    elif scenario_count == 1 or scenario_count == 4:
        if (
            scenario["instrument"]["CGtype"]
            != json_data["starlightSuppressionSystems"][1]["name"][-3:]
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " CGtype not matching."
            )
    else:
        if (
            scenario["instrument"]["CGtype"]
            != json_data["starlightSuppressionSystems"][2]["name"][-3:]
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " CGtype not matching."
            )

    # OpMode
    if scenario_count == 0 or scenario_count == 3:
        if (
            scenario["instrument"]["OpMode"] == "IMG"
            and json_data["scienceInstruments"][0]["name"][-6:] != "Imager"
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " OpMode not matching"
            )
    elif scenario_count == 1 or scenario_count == 4:
        if (
            scenario["instrument"]["OpMode"] == "IMG"
            and json_data["scienceInstruments"][1]["name"][-6:] != "Imager"
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " OpMode not matching"
            )
    else:
        if (
            scenario["instrument"]["OpMode"] == "SPEC"
            and json_data["scienceInstruments"][2]["name"][-4:] != "Spec"
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " OpMode not matching"
            )

    # pp_Factor_CBE
    if (
        scenario["instrument"]["pp_Factor_CBE"]
        != json_data["observingModes"][scenario_count]["pp_Factor_CBE"]
    ):
        errors.append(
            str(scenario_conversion[scenario_count])
            + " pp_Factor_CBE not matching."
        )

    # Kappa_c_HLCB1
    if "B1" in scenario["DataSpecification"]["ObservationCase"]:
        if "OPT" in scenario["DataSpecification"]["ObservationCase"]:
            if (
                np.abs(
                    scenario["TVACmeasured"]["Kappa_c_HLCB1"]
                    * scenario["TVACmeasured"]["CoreThput_HLCB1"]
                    - json_data["starlightSuppressionSystems"][0]["PSFpeak"]
                )
                > 1e-15
            ):

                errors.append(
                    str(scenario_conversion[scenario_count]) + " PSF peak not matching"
                )
        if "CON" in scenario["DataSpecification"]["ObservationCase"]:
            if (
                np.abs(
                    scenario["TVACmeasured"]["Kappa_c_HLCB1"]
                    * scenario["TVACmeasured"]["CoreThput_HLCB1"]
                    - json_data["starlightSuppressionSystems"][3]["PSFpeak"]
                )
                > 1e-15
            ):
                errors.append(
                    str(scenario_conversion[scenario_count]) + " PSF peak not matching"
                )

    # CoreThput_HLCB1
    if scenario_count == 0 or scenario_count == 3:
        if (
            scenario["TVACmeasured"]["CoreThput_HLCB1"]
            != json_data["starlightSuppressionSystems"][scenario_count]["core_thruput"]
        ):
            errors.append(
                str(scenario_conversion[scenario_count]) + " core thruput not matching"
            )

    scenario_count += 1

for e in errors:
    print(e)
