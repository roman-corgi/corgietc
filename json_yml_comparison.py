import json
import yaml
from pathlib import Path
import os
import cgi_noise
import corgietc

# Load directories from environment
DATA_DIR = Path(os.environ["CGI_NOISE_DATA_DIR"])
SCEN_DIR = DATA_DIR / "Scenarios"
SCRIPTFILE = Path(os.environ["CORGIETC_DATA_DIR"]) / "scripts" / "CGI_Noise.json"

def basename(path): 
    return Path(path).stem

def compare_values(yaml_val, json_val, scenario_name, field_name, errors):
    """Helper function to compare values and record errors"""
    if yaml_val != json_val:
        errors.append(f"{scenario_name} {field_name} not matching (YAML: {yaml_val}, JSON: {json_val})")

# Load YAML scenarios into a dictionary keyed by scenario name
yml_scenarios = {}
for yml_path in SCEN_DIR.glob("*.yml"):
    if yml_path.stem.startswith("REFERENCE_"):
        continue  # Skip reference scenario
    with open(yml_path, "r") as f:
        scenario_data = yaml.safe_load(f)
        scenario_name = scenario_data["DataSpecification"]["ObservationCase"]
        yml_scenarios[scenario_name] = scenario_data

# Load JSON file
with open(SCRIPTFILE, "r") as f:
    json_data = json.load(f)

# Create mappings for quick lookup
json_scenarios = {mode["Scenario"]: mode for mode in json_data["observingModes"]}
json_instruments = {inst["name"]: inst for inst in json_data["scienceInstruments"]}
json_systems = {sys["name"]: sys for sys in json_data["starlightSuppressionSystems"]}

errors = []

for scenario_name, yml_scenario in yml_scenarios.items():
    if scenario_name not in json_scenarios:
        errors.append(f"Scenario {scenario_name} not found in JSON")
        continue
        
    json_mode = json_scenarios[scenario_name]
    json_inst = json_instruments[json_mode["instName"]]
    json_sys = json_systems[json_mode["systName"]]
    
    # DataSpecification comparisons
    ds = yml_scenario["DataSpecification"]
    
    # Coronagraph_Data
    coro_files = [
        json_sys["occ_trans"],
        json_sys["core_mean_intensity"],
        json_sys["core_area"],
        json_sys["core_contrast"]
    ]
    if not all(ds["Coronagraph_Data"] == basename(f) for f in coro_files if f):
        errors.append(f"{scenario_name} Coronagraph_Data file not matching")
    
    # QE_Curve_Data
    qe_files = [json_inst["QE"], json_inst["DET_QE_Data"]]
    if not all(ds["QE_Curve_Data"] == basename(f) for f in qe_files if f):
        errors.append(f"{scenario_name} QE_Curve_Data not matching")
    
    # Detector_Data
    compare_values(ds["Detector_Data"], basename(json_inst["DET_CBE_Data"]), 
                  scenario_name, "Detector_Data", errors)
    
    # ContrastStability_Data
    cs_files = [
        json_sys["AvgRawContrast"],
        json_sys["ExtContStab"],
        json_sys["IntContStab"],
        json_sys["InitStatContrast"]
    ]
    if not all(ds["ContrastStability_Data"] == basename(f) for f in cs_files if f):
        errors.append(f"{scenario_name} ContrastStability_Data file not matching")
    
    # CS_Type
    csv_names = json_sys["csv_names"]
    if not all(ds["CS_Type"] == csv_names[key][:4] 
               for key in ["AvgRawContrast", "ExtContStab", "IntContStab", "SystematicC", "InitStatContrast"]):
        errors.append(f"{scenario_name} CS_Type not matching")
    
    # Throughput_Data
    compare_values(ds["Throughput_Data"], basename(json_sys["Throughput_Data"]), 
                  scenario_name, "Throughput_Data", errors)
    
    # StrayLight_Data
    compare_values(ds["StrayLight_Data"], basename(json_mode["StrayLight_Data"]), 
                  scenario_name, "StrayLight_Data", errors)
    
    # Instrument comparisons
    inst = yml_scenario["instrument"]
    
    # Diameter
    compare_values(inst["Diam"], json_data["pupilDiam"], scenario_name, "Diameter", errors)
    
    # Wavelength (convert YAML meters to nm for comparison)
    compare_values(inst["wavelength"] * 1e9, json_sys["lam"], scenario_name, "Wavelength", errors)
    
    # Bandwidth
    compare_values(inst["bandwidth"], json_sys["BW"], scenario_name, "Bandwidth", errors)
    
    # CGtype (last 3 chars of system name)
    compare_values(inst["CGtype"], json_sys["name"][-3:], scenario_name, "CGtype", errors)
    
    # OpMode
    op_mode = inst["OpMode"]
    sci_name = json_inst["name"]
    if (op_mode == "IMG" and "Imager" not in sci_name) or (op_mode == "SPEC" and "Spec" not in sci_name):
        errors.append(f"{scenario_name} OpMode not matching (YAML: {op_mode}, Instrument: {sci_name})")
    
    # pp_Factor_CBE
    compare_values(inst["pp_Factor_CBE"], json_mode["pp_Factor_CBE"], 
                  scenario_name, "pp_Factor_CBE", errors)
    
    # TVACmeasured comparisons
    if "TVACmeasured" in yml_scenario:
        tvac = yml_scenario["TVACmeasured"]
        
        # For B1 scenarios, check PSF peak product
        if "B1" in scenario_name:
            product = tvac["Kappa_c_HLCB1"] * tvac["CoreThput_HLCB1"]
            compare_values(product, json_sys["PSFpeak"], scenario_name, "PSF peak", errors)
        
        # Core throughput for B1 HLC
        if scenario_name.endswith("NFB1_HLC"):
            compare_values(tvac["CoreThput_HLCB1"], json_sys["core_thruput"], 
                          scenario_name, "CoreThput_HLCB1", errors)

# Print collected errors
if errors:
    print("Comparison errors found:")
    for error in errors:
        print(f" - {error}")
else:
    print("All scenarios match successfully!")