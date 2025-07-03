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

https://colab.research.google.com/github/roman-corgi/corgietc/blob/main/Notebooks/01_Anatomy_of_an_Integration_Time_Calculation.ipynb

Executing the setup notebook will create a directory called `corgietc` in your Google Drive.  At any point, you may remove this directory entirely from your Google drive and re-run the setup notebook to get a fresh installation of `corgietc`. This should not be necessary in ordinary usage, but may be helpful if the contents of this directory become corrupted for any reason. 

### Method 2: Local Installation

If you wish to run the ETC on your own computer, you must first install all required packages and download all relevant data.  

We **strongly** recommend use of a dedicated Python virtual environment.  The instructions below assume that you have Python (version 3.10 or higher) and pip installed and working on your machine. For help with that, start here: https://wiki.python.org/moin/BeginnersGuide/. We'll assume that Python and pip are executable as `python` and `pip`, respectively, but they might be called `python3` and `pip3` (or something else) on your particular system. You also need to have git installed (https://github.com/git-guides/install-git). These instructions are based on working in a terminal (macOS/Linux) or command prompt/PowerShell (Windows).

The steps below assume that you are on a POSIX system (e.g. macOS or Linux) and running a bash/zsh-style shell.  Where relevant, we have included links to documentation on how certain steps may differ on other systems/shells. 

1. Clone the corgietc repository (https://github.com/roman-corgi/corgietc) and the cgi_noise repository (https://github.com/roman-corgi/cgi_noise) to your computer (see here for help on this: https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository). Place these in any directory on disk - we will refer to the directory where you cloned these two repositories as `$CODE_ROOT`.

2. Create a python virtual environment (we'll call ours `corgietc_env` but you can replace this with any name you like). We will place our virtual environment directory in the same place as our cloned repositories, but you can put the virtual environment anywhere on disk (just keep track of where things are). In a terminal/command prompt/powershell/etc run:
   
   ```
   cd $CODE_ROOT
   python -m venv corgietc_env
   ```
   
   >**Warning** The `cd $CODE_ROOT` command will work as written only on POSIX systems (not Windows, which has a different environment variable syntax) and *only* if you have exported the environment variable `CODE_ROOT` to point at the directory where you wish to store things. You can do so in a bash/zsh shell by first executing:
   
   ```
   export CODE_ROOT="/full/path/to/CODE_ROOT"
   ```

   This step is entirely optional - if you prefer not to set the environment variable, simply replace all instances of `$CODE_ROOT` in these instructions with the full path to your code root directory. If you'd like to use environment variables in a different shell/system, or make them persistent across multiple sessions, see the relevant documentation for your particular system configuration. 
   
3. Activate the environment. On macOS/Linux (if using a bash/zsh shell):

    ```source corgietc_env/bin/activate```

   For Windows, see https://docs.python.org/3/library/venv.html.  If using a different shell on macOS/Linux, see here: https://docs.python.org/3/library/venv.html#how-venvs-work.  To determine which shell you're using on a POSIX system, use the command `echo $SHELL`. For tcsh/csh shells, the relevant command is:

   ```source corgietc_env/bin/activate.csh```

4. In the same terminal with the active virtual environment, navigate to the cloned/downloaded `cgi_noise` repository and run the install command:

    ```
    cd $CODE_ROOT/cgi_noise
    pip install .
    ```

   Do not omit the period at the end - it is important. To check whether you are in a session with the virtual environment active, look at the left-hand side of your command prompt - if the environment is active, it should show the environment name in parentheses (e.g., `(corgietc_env)`). 

5. Once `cgi_noise` successfully installs, navigate to the cloned/downlaoded `corgietc` repository and run the install command:

    ```
    cd $CODE_ROOT/corgietc
    pip install .
    ```
     
6. You will also need a Jupyter environment (we recommend jupyter-lab) to execute the notebooks and matplotlib to plot results. These can be installed by running (still in the same session with the active virtual environment):

    ```pip install jupyterlab matplotlib```

7. Restart your virtual environment.  While possibly not necessary, on many systems, the JupyterLab executable will not be automatically found unless you restart the environment or otherwise rescan your path.  The simplest thing to do is to restart the virtual environment:

    ```
    deactivate
    source $CODE_ROOT/corgietc_env/bin/activate.csh
    ```

8. In the terminal with the active virtual environment session, navigate to the `Notebooks` subdirectory of the repository and start JupyterLab:

    ```
    cd $CODE_ROOT/corgietc/Notebooks
    jupyter-lab
    ```

   >**Warning** Note that the JupyterLab executable has a different from the package.  The JupyterLab package is `jupyterlab` (no hyphen) whereas the executable is `jupyter-lab` (with hyphen). The command `jupyter lab` (space instead of hyphen) is an equivalent variant.
    
9. A JupyterLab instance should start up in your default browser.  In the left-hand pane double-click on any of the Notebooks (`ipynb` files) to open it.  If you do not see any notebooks listed, make sure you have the file viewer up by clicking on the Folder icon on the left-hand side of the screen. Skip the code blocks of any notebook that are marked as for Collab execution only.

10. To stop JupyterLab, type `ctrl+c` in the terminal where it is running and then hit `ctrl+c` again (or type `y` at the prompt). To deactivate the virtual environment type `deactivate` at the prompt.  Next time you want to run the Notebooks again, you activate the environment again, navigate to the Notebooks directory and run `jupyter-lab`:

    ```
    source $CODE_ROOT/corgietc_env/bin/activate.csh
    cd $CODE_ROOT/corgietc/Notebooks
    jupyter-lab
    ```

   >**Warning** If you set the `CODE_ROOT` environment variable in one terminal session, it will not be persist across different terminal sessions.  You must define it each time, or add it to your shell configuration file (see your system's documentation for details) or just use the full path to the code root directory in all commands. 

### Updating corgietc

If using `corgietc` via Google Colab, updates occur automatically - no user action is required. Each time you run one of the notebooks, you will automatically be using the latest code version.  Occassionally, Google Colab may hold on to a stale notebook version.  You can fix this by deleting the current runtime (in the `Runtime` menu) and clearing your browser cache. 

For local installations, if using a virtual environment, make sure that it is active (see installation instruction, above, on how to activate/deactivate environments). Then:

```
cd $CODE_ROOT/cgi_noise
git pull
pip install --upgrade .

cd ../corgietc
git pull
pip install --upgrade .
```

The same caveats apply as in the installation instructions - if you have not defined the `$CODE_ROOT` environment variable, simply replace it with the full path to the folder in which you cloned the github repositories. 


