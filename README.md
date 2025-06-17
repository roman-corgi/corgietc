# corgietc
Roman Coronagraph Instrument Exposure Time Calculator

## Using corgietc

The ETC can be run entirely via your browser, with no local installation required.  This is recommended for anyone who wants to use the ETC, but not edit it or do other development tasks. Alternatively, you can install everything locally on your own machine, which allows you to run the ETC without an internet connection (although an internet connection is still required for initial runs in order to download additional data required for various calculations).

### Method 1: Google Colab (No Local Installation)

Navigate to:

https://colab.research.google.com/github/roman-corgi/corgietc/blob/main/Notebooks/00_Google_Colab_Setup.ipynb

Ensure that you are logged in with the Google account you wish to use (data will be written to the Google Drive associated with this account).  You can check which account you are logged into by clicking on the user icon in the top right-hand corner of the page. 

Execute all of the cells in the notebook, responding to any pop-up prompts along the way (see the notebook for more detailed instructions). Note that you only need to run this notebook **once** (even if you log out/close the browser instance, the files written to your Google Drive will be persistent).

After successfully executing the setup notebook, navigate to any of the other notebooks.  You may do so directly via Colab by clicking on the `File` menu (at the top left of the page) and selecting `Open Notebook`. Select `GitHub` in the left-hand pane of the dialog that appears, make sure that `roman-corgi` is entered in the top text box and select the `roman-corgi/corgietc` repository from the dropdown menu.  A list of all available notebooks should appear. You may see a `Leave Page` (or equivalent) prompt - it is okay to proceed. 

Alternatively, you may access any notebooks directly by pre-pending `https://colab.research.google.com/github/roman-corgi/corgietc/blob/main/Notebooks/` to the name of any of the notebooks.  So, for example, you can directly access the first notebook by navigating to:

https://colab.research.google.com/github/roman-corgi/corgietc/blob/main/Notebooks/01_Anatomy_of_an_Integration_Time_Calculation.ipynb 02_Getting_to_Know_Your_Observing_Mode.ipynb

Executing the setup notebook will create a directory called `corgietc` in your Google Drive.  At any point, you may remove this directory entirely from your Google drive and re-run the setup notebook to get a fresh installation of `corgietc`. This should not be necessary in ordinary usage, but may be helpful if the contents of this directory become corrupted for any reason. 

### Method 2: Local Installation

If you wish to run the ETC on your own computer, you must first install all required packages and download all relevant data.  

We **strongly** recommend use of a dedicated Python virtual environment.  The instructions below assume that you have Python (version 3.10 or higher) and pip installed and working on your machine. For help with that, start here: https://wiki.python.org/moin/BeginnersGuide/. We'll assume that Python and pip are executable as `python` and `pip`, respectively, but they might be called `python3` and `pip3` (or something else) on your particular system. These instructions are based on working in a terminal (macOS/Linux) or command prompt/PowerShell (Windows). You also need to have git installed (https://github.com/git-guides/install-git) and executable as `git`. 

1. Clone the corgietc repository (https://github.com/roman-corgi/corgietc) and the cgi_noise repository (https://github.com/roman-corgi/cgi_noise) to your computer (see here for help on this: https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository).

2. Create a python virtual environment (we'll call ours `corgietc` but you can replace this with any name you like). In a terminal/command prompt/powershell/etc, navigate to a directory where you want to create the virtual environment and run:
   
   ```python -m venv corgietc```
   
3. Activate the environment. On macOS/Linux:

    ```source corgietc/bin/activate```

   For Windows, see https://docs.python.org/3/library/venv.html

5. In the same terminal with the active virtual environment, navigate to the cloned/downloaded `cgi_noise` repository.  From the top level directory of the repository (the one that contains the file `pyproject.toml`) run:

    ```pip install .```

   Do not omit the period at the end - it is important. 

7. Once `cgi_noise` successfully installs, navigate to the cloned/downlaoded `corgietc` repository.  From the top level directory of the repository (the one that contains the file `pyproject.toml`) run:

    ```
    pip install .
    ```
    
    This will install all of the python packages required by the demo notebooks.
 
8. You will also need a Jypyter environment to execute the notebooks.  We recommend jupyter-lab, which can be installed by running (from any directory):

    ```pip install jupyter-lab```


9. Navigate to the `Notebooks` subdirectory of the repository (this should just be `cd Notebooks` from where you were last) and then start JupyterLab by running `jupyter-lab`

10. Skip the code blocks of any notebook that are marked as for Collab execution only.

11. To stop JupyterLab, type `ctrl+c` in the terminal where it is running and then hit `ctrl+c` again (or type `y` at the prompt). To deactivate the virtual environment just type `deactivate` at the prompt.  Next time you want to run the tutorials again, you just activate the environment again, navigate to the Notebooks directory and run `jupyter-lab`

>**Warning**
>There appears to be an issue (at least on macOS) where if you already have jupyter-lab installed in a system path, it will be executed initially instead of the one you install in your virtual environment.  A simple fix is to deactivate and re-activate the virtual environment after you run the initial pip installation (i.e., between steps 6 and 7).  Environments can be deactivated by running `deactivate` in the terminal where the environment is active, and then re-activated by sourcing the same file as in step 3. 

### Updating corgietc

If using `corgietc` via Google Colab, updates occur automatically - no user action is required. Each time you run one of the notebooks, you will automatically be using the latest code version.  Occassionally, Google Colab may hold on to a stale notebook version.  You can fix this by deleting the current runtime (in the `Runtime` menu) and clearing your browser cache. 

For local installations, run `git pull` in the repository directory, and then upgrade the installed version in your virtual environment by executing 

```pip install --upgrade .```

from the repository top level.


