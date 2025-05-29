import os
import warnings

import astropy.units as u
import numpy as np
from scipy import interpolate
import pandas as pd

from EXOSIMS.OpticalSystem.Nemati import Nemati


class CGI_Noise(Nemati):

    def __init__(self, **specs):
        # call upstream init
        Nemati.__init__(self, **specs)
