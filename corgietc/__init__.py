import os
import importlib.resources

name = "corgietc"
__version__ = "0.1.0"

# identify data directory and add to environment variables for this session
datapath = importlib.resources.files("corgietc").joinpath("data")
assert datapath.exists(), (
    "Could not identify corgietc datapath. Check that the "
    "corgietc installation was fully successful."
)
os.environ["CGI_PERF_DIR"] = str(datapath)
