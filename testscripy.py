import corgietc
import os
import json
import EXOSIMS.Prototypes.TargetList
import EXOSIMS.Prototypes.TimeKeeping
import copy
import astropy.units as u
import numpy as np

scriptfile = os.path.join(os.environ["CORGIETC_DATA_DIR"], "scripts", "CGI_Noise.json")
with open(scriptfile, "r") as f:
    specs = json.loads(f.read())
    
    
TL = EXOSIMS.Prototypes.TargetList.TargetList(**copy.deepcopy(specs))
OS = TL.OpticalSystem

TK = EXOSIMS.Prototypes.TimeKeeping.TimeKeeping(missionLife = 5.25)   # 63 months in years is 5.25, 21 months is 1.75
TK.allocate_time(21*30.4375*u.d);

sInds = 0
fZ = np.repeat(TL.ZodiacalLight.fZ0, 1)

mode = OS.observingModes[0]

JEZ = TL.JEZ0[mode["hex"]]/(4.1536**2)

dMag = np.array([17.5])

WA = np.array([7.5]) * (mode["lam"]/OS.pupilDiam).to(u.arcsec, equivalencies=u.dimensionless_angles())

intTimes = OS.calc_intTime(TL, sInds, fZ, JEZ, dMag, WA, mode, TK=TK)
intTimes.to(u.s)

dMags = np.linspace(16, 23, 100)
intTimes = OS.calc_intTime(TL, [sInds]*len(dMags), fZ, JEZ, dMags, WA, mode, TK=TK)

soln = OS.int_time_denom_obj(dMag, TL, sInds, fZ, JEZ, WA, mode, None)

print(soln)

print(OS.int_time_denom_obj.__qualname__)
print(OS.int_time_denom_obj.__module__)

