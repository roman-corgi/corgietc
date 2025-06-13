from corgietc.corgietc import corgietc
import astropy.units as u
import os, json, copy
from pathlib import Path
import yaml
import math
import numpy as np
from dataclasses import dataclass
import cgi_noise.cginoiselib as fl
import cgi_noise.unitsConstants as uc
from cgi_noise.loadCSVrow import loadCSVrow
from prettytable import PrettyTable
import pandas as pd
import EXOSIMS.Prototypes.TargetList
import EXOSIMS.Prototypes.TimeKeeping
import matplotlib.pyplot as plt

#corgietc
scriptfile = os.path.join(os.environ["CORGIETC_DATA_DIR"], "scripts", "CGI_Noise.json")
with open(scriptfile, "r") as f:
    specs = json.loads(f.read())

OS = corgietc(**copy.deepcopy(specs))

# === edit obs_params["scenario"] for each scenario ===

#possible scenarios: OPT_IMG_NFB1_HLC: science instrument [0], starlight suppression[0], observing modes [0]
#OPT_IMG_WFB4_SPC: science instrument [1], starlight suppression[1], observing modes [1]
#OPT_SPEC_NFB3_SPC: science instrument [2], starlight suppression[2], observing modes [2]
#CON_IMG_NFB1_HLC: science instrument [0], starlight suppression[3], observing modes [3]
#CON_IMG_WFB4_SPC: science instrument [1], starlight suppression[4], observing modes [4]
#CON_SPEC_NFB3_SPC: science instrument [2], starlight suppression[5], observing modes [5]

obs_params = {
        "scenario": "OPT_IMG_NFB1_HLC.yml",
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
 
if obs_params["scenario"] == 'OPT_IMG_WFB4_SPC.yml' or obs_params["scenario"] == 'CON_IMG_WFB4_SPC.yml':
    obs_params["target_params"]['sma_AU'] = 8

DATA_DIR = Path(os.environ["CGI_NOISE_DATA_DIR"])
SCEN_DIR = DATA_DIR / "Scenarios"

scenario_filename = obs_params["scenario"]
scenario_path = SCEN_DIR / scenario_filename
with open(SCEN_DIR / scenario_filename, "r") as file:
    config = yaml.safe_load(file)
target_params = obs_params["target_params"]
SNRdesired = obs_params["snr"]

for i in range(6):
    if OS.observingModes[i]["Scenario"]  == os.path.splitext(obs_params["scenario"])[0]:
        scenario_idx = i


#cgi_noise
@dataclass
class corePhotonRates:
    planet: float
    speckle: float
    locZodi: float
    exoZodi: float
    straylt: float
    total: float = 0.0


ObservationCase = config['DataSpecification']['ObservationCase']
    
DPM = config['instrument']['Diam']
lam = config['instrument']['wavelength']
lamD = lam / DPM
intTimeDutyFactor = config['instrument']['dutyFactor']
opMode = config['instrument']['OpMode']
bandWidth = config['instrument']['bandwidth']

print(f"Central wavelength: {lam / uc.nm:.1f} nm, with {bandWidth * 100:.0f}% BW\n")
print(f"Lambda/D: {lamD / uc.mas:.3f} mas")

target = fl.Target(**target_params)
sep_mas = target.phaseAng_to_sep(target.sma_AU, target.dist_pc, target.phaseAng_deg)
target.albedo = target.albedo_from_geomAlbedo(target.phaseAng_deg, target.geomAlb_ag)

print(f"Separation: {sep_mas:.0f} mas")
print(f"Albedo: {target.albedo:.3f}")

filenameList = fl.getScenFileNames(config, DATA_DIR)
CG_Data, QE_Data, DET_CBE_Data, STRAY_FRN_Data, THPT_Data, CAL_Data, CS_Data = fl.loadCSVs(filenameList)
CS_Type = config['DataSpecification']['CS_Type']

IWA, OWA = fl.workingAnglePars(CG_Data, CS_Data)
planetWA = sep_mas * uc.mas / lamD

# from prettytable import PrettyTable
table = PrettyTable()
table.field_names = ['planet WA', 'phase', 'dist', 'sma', 'sep', 'lam/D', "IWA", "OWA"]
table.add_row([f"{planetWA:.2f}", f"{target.phaseAng_deg:.2f}", f"{target.dist_pc:.2f}", f"{target.sma_AU:.2f}", f"{sep_mas:.2f}", f"{lamD/uc.mas:.2f}", f"{IWA:.2f}", f"{OWA:.2f}"])
print(table)

tol = 0.05
if (IWA - tol) <= planetWA <= IWA:
    planetWA = IWA
elif OWA <= planetWA <= (OWA + tol):
    planetWA = OWA
elif planetWA < (IWA - tol) or planetWA > (OWA + tol):
    raise ValueError(f"Planet WA={planetWA:.1f} outside of IWA={IWA:.1f}, OWA={OWA:.1f}.")

selDeltaC, AvgRawC, SystematicC, initStatRaw, IntContStab, ExtContStab = fl.contrastStabilityPars(CS_Type, planetWA, CS_Data)

cg = fl.coronagraphParameters(CG_Data.df, config, planetWA, DPM)
f_SR, CritLam, detPixSize_m, mpix = fl.getFocalPlaneAttributes(opMode, config, DET_CBE_Data, lam, bandWidth, DPM, cg.CGdesignWL, cg.omegaPSF, DATA_DIR)

inBandFlux0_sum, inBandZeroMagFlux, starFlux = fl.getSpectra(target, lam, bandWidth, DATA_DIR)
print(f"Star Flux = {starFlux:.3e} ph/s/m^2")
TimeonRefStar_tRef_per_tTar = 0.25

fluxRatio = target.albedo * (target.radius_Rjup * uc.jupiterRadius / (target.sma_AU * uc.AU))**2
planetFlux = fluxRatio * starFlux
print(f"Planet Flux Ratio = {fluxRatio:.2e}\nPlanet Flux = {planetFlux:.3f} ph/s/m^2")

magLocalZodi = config['instrument']['LocZodi_magas2']
magExoZodi_1AU = config['instrument']['ExoZodi_magas2']
absMag = target.v_mag - 5 * math.log10(target.dist_pc / 10)
locZodiAngFlux = inBandZeroMagFlux * 10**(-0.4 * magLocalZodi)
exoZodiAngFlux = inBandZeroMagFlux * 10**(-0.4 * (absMag - uc.sunAbsMag + magExoZodi_1AU)) / target.sma_AU**2 * target.exoZodi

thput, throughput_rates = fl.compute_throughputs(THPT_Data, cg, "uniform")
Acol = (np.pi / 4) * DPM**2
stray_ph_s_mm2 = fl.getStrayLightfromfile(ObservationCase, 'CBE', STRAY_FRN_Data)
stray_ph_s_pix = stray_ph_s_mm2 * (1 / uc.mm**2) * detPixSize_m**2

cphrate = corePhotonRates(
    planet=planetFlux * throughput_rates["planet"] * Acol,
    speckle=starFlux * AvgRawC * cg.PSFpeakI * cg.CGintmpix * throughput_rates["speckle"] * Acol,
    locZodi=locZodiAngFlux * cg.omegaPSF * throughput_rates["local_zodi"] * Acol,
    exoZodi=exoZodiAngFlux * cg.omegaPSF * throughput_rates["exo_zodi"] * Acol,
    straylt=stray_ph_s_pix * mpix
)
cphrate.total = sum([cphrate.planet, cphrate.speckle, cphrate.locZodi, cphrate.exoZodi, cphrate.straylt])

ENF, effReadnoise, frameTime, dQE, QE_img = fl.compute_frame_time_and_dqe(0.1, 3, 100, True, QE_Data, DET_CBE_Data, lam, mpix, cphrate.total)
print(f"Calculated Frame Time: {frameTime:.2f} s")
print(f"QE in the image area: {QE_img:.3f}")
print(f"Detected Quantum Efficiency (dQE): {dQE:.3f}")
print(f"Excess Noise Factor (ENF): {ENF:.2f}")
print(f"Core fraction used in the SNR region for mode {config['DataSpecification']['ObservationCase']}: f_SR: {f_SR:.3f}")

detNoiseRate = fl.detector_noise_rates(DET_CBE_Data, 21, frameTime, mpix, True)

rdi_penalty = fl.rdi_noise_penalty(inBandFlux0_sum, starFlux, TimeonRefStar_tRef_per_tTar, 'a0v', 2.26)
k_sp = rdi_penalty['k_sp']
k_det = rdi_penalty['k_det']
k_lzo = rdi_penalty['k_lzo']
k_ezo = rdi_penalty['k_ezo']

nvRatesCore, residSpecRate = fl.noiseVarianceRates(
    cphrate, QE_img, dQE, ENF, detNoiseRate,
    k_sp, k_det, k_lzo, k_ezo,
    f_SR, starFlux, selDeltaC,
    config['instrument']['pp_Factor_CBE'], cg,
    throughput_rates['speckle'], Acol
)

# SNR threshold from caller
# SNRdesired is now passed in as an argument
planetSignalRate = f_SR * cphrate.planet * dQE
timeToSNR, criticalSNR = fl.compute_tsnr(SNRdesired, planetSignalRate, nvRatesCore, residSpecRate)

csfilename = None
for filepath in filenameList:
    base = os.path.basename(filepath)
    if base.startswith("CS_"):
        csfilename = os.path.splitext(base)[0]
        break


print("\nTotal noise variance rate beakdown:")
table = PrettyTable()
table.field_names = ['planet', 'speckle', 'local Zodi', 'exo Zodi', 'Stray']
table.add_row([f"{nvRatesCore.planet:.4f}",f"{nvRatesCore.speckle:.4f}",f"{nvRatesCore.locZodi:.4f}",f"{nvRatesCore.exoZodi:.4f}",f"{nvRatesCore.straylt:.3e}",])
print(table)

print("\nContrast Stability Numbers:") 
table = PrettyTable()
table.field_names = ['CS case', 'DeltaC',  'External CS', 'Internal CS', 'Systematic', 'Avg Raw', 'initStatRaw']
table.add_row([f"{csfilename}", f"{selDeltaC:.2e}", f"{ExtContStab:.2e}", f"{IntContStab:.2e}", f"{SystematicC:.2e}", f"{AvgRawC:.2e}", f"{initStatRaw:.2e}"])
print(table)

print(f"\nCalculation ingredients for Time to SNR = {SNRdesired:.1f}:")
table = PrettyTable()
table.field_names = [ 'planet rate', 'total noise rate', 'resid specle rate']
table.add_row([ f"{planetSignalRate:.3f}", f"{nvRatesCore.total:.3f}", f"{residSpecRate:.4f}"])
print(table)

print(f"\nTarget SNR = {SNRdesired:.1f} \nCritical SNR = {criticalSNR:.2f}")
print(f"Time to SNR = {timeToSNR:.1f} seconds or {timeToSNR/uc.hour:.3f} hours\n")

if timeToSNR > intTimeDutyFactor * 100 * uc.hour:
    print(f"Warning: Time to SNR ({timeToSNR/uc.hour:.2f} hrs) exceeds usable integration time ({(intTimeDutyFactor * 100 * uc.hour)/uc.hour:.2f} hrs).")
elif timeToSNR <= 0:
    print(f"Warning: Calculated Time to SNR is not positive ({timeToSNR/uc.hour:.2f} hrs). Check input parameters and intermediate calculations.")
WA = planetWA*lamD*u.rad  #Working angle in LAMBDA/D units


# === core_thruput, PSFpeak, core_area, lam, core_mean_intensity, occ_trans, core_contrast ===
#core area is converted to arcsec^2
table = PrettyTable()
table.field_names = ['System', 'Occulter Transmission', 'Core Thruput', 'Core Mean Intensity', 'Core Area', 'Core Contrast', 'PSF Peak', 'lam']
table.add_row( ['cgi_noise',f"{cg.CG_occulter_transmission}", f"{cg.CGcoreThruput}", f"{cg.CGintensity}", f"{cg.omegaPSF*lamD**2*4.255*1E10}", f"{cg.CGcontrast}", f"{cg.PSFpeakI}", f"{cg.CGdesignWL*1E9}"] )
table.add_row( ['corgietc',f"{OS.starlightSuppressionSystems[scenario_idx]["occ_trans"](OS.observingModes[scenario_idx]["lam"], WA)}",
                f"{OS.starlightSuppressionSystems[scenario_idx]["core_thruput"](OS.observingModes[scenario_idx]["lam"], WA)}", 
                f"{OS.starlightSuppressionSystems[scenario_idx]["core_mean_intensity"](OS.observingModes[scenario_idx]["lam"], WA)[0]}", 
                f"{OS.starlightSuppressionSystems[scenario_idx]["core_area"](OS.observingModes[scenario_idx]["lam"], WA)}", 
                f"{OS.starlightSuppressionSystems[scenario_idx]["core_contrast"](OS.observingModes[scenario_idx]["lam"], WA)}", 
                f"{OS.starlightSuppressionSystems[scenario_idx]["PSFpeak"](OS.observingModes[scenario_idx]["lam"], WA)}",
                f"{OS.starlightSuppressionSystems[scenario_idx]["lam"]}"] )
print(table)

# === AvgRawContrast, SystematicC, InitStatContrast, IntContStab, ExtContStab ===
#cgi_noise is in ppb
try:
    corgi_systematicC = OS.starlightSuppressionSystems[scenario_idx]["SystematicC"](OS.observingModes[scenario_idx]["lam"], WA)
except KeyError:
    corgi_systematicC = 0

table = PrettyTable()
table.field_names = ['System', 'Raw Contrast', 'Systematic Cont', 'initStatRawContrast', 'IntContStab', 'ExtContStab']
table.add_row( ['cgi_noise',f"{AvgRawC*1E9}", f"{SystematicC*1E9}", f"{initStatRaw*1E9}", f"{IntContStab*1E9}", f"{ExtContStab*1E9}"] )
table.add_row( ['corgietc',f"{OS.starlightSuppressionSystems[scenario_idx]["AvgRawContrast"](OS.observingModes[scenario_idx]["lam"], WA)[0]}",
                f"{corgi_systematicC}", 
                f"{OS.starlightSuppressionSystems[scenario_idx]["InitStatContrast"](OS.observingModes[scenario_idx]["lam"], WA)[0]}", 
                f"{OS.starlightSuppressionSystems[scenario_idx]["IntContStab"](OS.observingModes[scenario_idx]["lam"], WA)[0]}", 
                f"{OS.starlightSuppressionSystems[scenario_idx]["ExtContStab"](OS.observingModes[scenario_idx]["lam"], WA)[0]}"] )
print(table)

# === f_SR, CritLam, pixelSize ===
# CritLam in m and pixel scale in mas for cgi_noise
table = PrettyTable()
table.field_names = ['System', 'f_SR', 'CritLam', 'pixelSize']
table.add_row( ['cgi_noise',f"{f_SR}", f"{CritLam*1E9}", f"{detPixSize_m}"] )
table.add_row( ['corgietc', f"{OS.observingModes[scenario_idx%3]['f_SR']}",
                f"{OS.scienceInstruments[scenario_idx%3]['CritLam']}", 
                f"{OS.scienceInstruments[scenario_idx%3]['pixelSize']}"] )
print(table)

# === pupilArea, pupilDiam, compbeamD, fnlFocLen, Rlamsq, Rlam, Rconst, fnumber, Rs, strayLight, stray_ph_s_pix ===
FocalPlaneAtt = loadCSVrow(DATA_DIR / 'instrument' / 'CONST_SNR_FPattributes.csv')
AmiciPar = loadCSVrow(DATA_DIR / 'instrument' / 'CONST_Amici_parameters.csv')
compbeamD_m = AmiciPar.df.loc[0, 'compressd_beam_diamtr_m']
fnlFocLen_m = AmiciPar.df.loc[0, 'final_focal_length_m']
Rlamsq = AmiciPar.df.loc[0, 'lam_squared']
Rlam = AmiciPar.df.loc[0, 'lam']
Rconst = AmiciPar.df.loc[0, 'constant']
Fno = fnlFocLen_m / compbeamD_m
if config['instrument']['R_required'] != 'NaN':
    resolution = config['instrument']['R_required'] 
else:
    resolution = 1

table = PrettyTable()
straylight_path = os.path.expandvars(OS.observingModes[scenario_idx]['StrayLight_Data'])
straylight = pd.read_csv(straylight_path, comment = "#")

table.field_names = ['System', 'pupilArea', 'pupilDiam', 'compbeamD', 'fnlFocLen', 'Rlamsq', 'Rlam', 'Rconst', 'fnumber', 'Rs', 'strayLight', 'stray_ph_s_pix']
table.add_row( ['cgi_noise',f"{Acol}", f"{DPM}", f"{compbeamD_m}", f"{fnlFocLen_m}", f"{Rlamsq}", f"{Rlam}", f"{Rconst}", f"{Fno}", f"{resolution}", f"{stray_ph_s_mm2}", f"{stray_ph_s_pix}"])
table.add_row( ['corgietc', f"{OS.pupilArea}",
                f"{OS.pupilDiam}", 
                f"{OS.observingModes[scenario_idx]['inst']['compbeamD']}", 
                f"{OS.observingModes[scenario_idx]['inst']['fnlFocLen']}",
                f"{OS.observingModes[scenario_idx]['inst']['Rlamsq']}",
                f"{OS.observingModes[scenario_idx]['inst']['Rlam']}",
                f"{OS.observingModes[scenario_idx]['inst']['Rconst']}",
                f"{OS.observingModes[scenario_idx]['inst']['fnumber']}",
                f"{OS.observingModes[scenario_idx]['inst']['Rs']}", 
                f"{straylight[OS.observingModes[scenario_idx]["Scenario"]][0]}", 
                f"{OS.observingModes[scenario_idx]['stray_ph_s_pix']}"])
print(table)

# === Percentage differences ===
per_occtrans = np.abs(OS.starlightSuppressionSystems[scenario_idx]["occ_trans"](OS.observingModes[scenario_idx]["lam"], WA) - cg.CG_occulter_transmission)/(cg.CG_occulter_transmission)*100
per_thruput = np.abs(OS.starlightSuppressionSystems[scenario_idx]["core_thruput"](OS.observingModes[scenario_idx]["lam"], WA) - cg.CGcoreThruput)/(cg.CGcoreThruput)*100
per_intensity = np.abs(OS.starlightSuppressionSystems[scenario_idx]["core_mean_intensity"](OS.observingModes[scenario_idx]["lam"], WA)[0] - cg.CGintensity)/(cg.CGintensity)*100
per_area = np.abs(OS.starlightSuppressionSystems[scenario_idx]["core_area"](OS.observingModes[scenario_idx]["lam"], WA).value - cg.omegaPSF*lamD**2*4.255*1E10)/(cg.omegaPSF*lamD**2*4.255*1E10)*100
per_contrast = np.abs(OS.starlightSuppressionSystems[scenario_idx]["core_contrast"](OS.observingModes[scenario_idx]["lam"], WA) - cg.CGcontrast)/(cg.CGcontrast)*100
per_lam = np.abs(OS.starlightSuppressionSystems[scenario_idx]["lam"].value - cg.CGdesignWL*1E9)/(cg.CGdesignWL*1E9)*100
per_rawcontrast = np.abs(OS.starlightSuppressionSystems[scenario_idx]["AvgRawContrast"](OS.observingModes[scenario_idx]["lam"], WA)[0] - AvgRawC*1E9)/(AvgRawC*1E9)*100
per_initStatRawContrast = np.abs(OS.starlightSuppressionSystems[scenario_idx]["InitStatContrast"](OS.observingModes[scenario_idx]["lam"], WA)[0] - initStatRaw*1E9)/(initStatRaw*1E9)*100
per_IntContStab = np.abs(OS.starlightSuppressionSystems[scenario_idx]["IntContStab"](OS.observingModes[scenario_idx]["lam"], WA)[0] - IntContStab*1E9)/(IntContStab*1E9)*100
per_ExtContStab = np.abs(OS.starlightSuppressionSystems[scenario_idx]["ExtContStab"](OS.observingModes[scenario_idx]["lam"], WA)[0] - ExtContStab*1E9)/(ExtContStab*1E9)*100


print("Occulter Transmission % difference: " + f"{per_occtrans}")
print("Core Thruput % difference: " + f"{per_thruput}")
print("Core Mean Intensity % difference: " + f"{per_intensity}")
print("Core Area % difference: " + f"{per_area}")
print("Core Contrast % difference: " + f"{per_contrast}")
print("Lam % difference: " + f"{per_lam}")
print("Avg Raw Contrast % difference: " + f"{per_rawcontrast}")
print("initStatRawContrast % difference (can ignore): " + f"{per_initStatRawContrast}")
print("IntContStab % difference (can ignore): " + f"{per_IntContStab}")
print("ExtContStab % difference (can ignore): " + f"{per_ExtContStab}")

## === Integration time ===

TL = EXOSIMS.Prototypes.TargetList.TargetList(**copy.deepcopy(specs))
OS = TL.OpticalSystem

# create a Timekeeping object and advance the mission time a bit
TK = EXOSIMS.Prototypes.TimeKeeping.TimeKeeping(missionLife = 5.25)   # 63 months in years is 5.25, 21 months is 1.75
TK.allocate_time(21*30.4375*u.d)
sInds = 0
fZ = np.repeat(TL.ZodiacalLight.fZ0, 1)
mode = OS.observingModes[scenario_idx]
JEZ = TL.JEZ0[mode["hex"]]/(obs_params["target_params"]['sma_AU']**2)
dMag = np.array([17.5])
WA = np.array([7.5]) * (mode["lam"]/OS.pupilDiam).to(u.arcsec, equivalencies=u.dimensionless_angles())
intTimes = OS.calc_intTime(TL, sInds, fZ, JEZ, dMag, WA, mode, TK=TK)
print("Integration times: " + f"{intTimes.to(u.s)}")
dMags = np.linspace(16, 23, 100)
intTimes = OS.calc_intTime(TL, [sInds]*len(dMags), fZ, JEZ, dMags, WA, mode, TK=TK)
plt.figure()
plt.semilogy(dMags, intTimes)
plt.xlabel(r"Target $\Delta$mag")
plt.ylabel(f"Integration Time ({intTimes.unit})")
plt.title(fr"Maximum Achievable $\Delta$mag = {np.min(dMags[np.isnan(intTimes)]) :.2f}")
plt.show()