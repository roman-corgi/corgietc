import os
import importlib.resources

name = "corgietc"
__version__ = "1.0.1"

# identify data directory and add to environment variables for this session
datapath = importlib.resources.files("corgietc").joinpath("data")
assert datapath.exists(), (
    "Could not identify corgietc datapath. Check that the "
    "corgietc installation was fully successful."
)
os.environ["CORGIETC_DATA_DIR"] = str(datapath)
