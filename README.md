# TraceAnalysis

Software for single-molecule trace analysis

## Prerequisites

### Miniconda python 3.7 64-bit

Install miniconda python version 3.7 64-bit from https://docs.conda.io/en/latest/miniconda.html  
    IMPORTANT: Under advanced options tick the box "Add Anaconda to my PATH enivronment variable"

### Git

Install git from https://git-scm.com/download/win  
    Use all standard settings  
    IMPORTANT: On the page "Adjusting your PATH enivronment" make sure that you choose "Git from the command line and also from 3rd party software"

## Installation

*Advanced installation instructions can be found below*

- Get the `Setup.cmd` file for Windows or the `Setup.command` file for Mac (and probably Linux)
- Place this file in a directory where you want the code to be located
- Double click on the setup file
  *Note for Mac: If you download the .command file from the internet you have to give it execute permission using `chmod 755 Setup.command`*  
- Give your netid and password to log into gitlab and download the code
- The code will be in a subfolder called "traceanalysis", you can move this folder anywhere you like
- The original `Setup.cmd` or `Setup.command` file can be deleted.

## Using the GUI
Start the user interface by double clicking the file `Start GUI.cmd` for Windows or `Start GUI.command` for Mac

## Updating
Double click on the file `Update.cmd` for Windows or `Update.command` for Mac

---

## Installation - advanced

Clone the git repository:

```batch
clone https://gitlab.tudelft.nl/ivoseverins/traceanalysis.git
```

Create a conda environment named `trace_analysis` based on the `environment.yml` file. (Make sure you are in the main folder)

```batch
conda env create --file environment.yml
```

## Using the software - advanced

Either start the GUI by (use `pythonw` instead of `python` on Mac)

```batch
python GUI.py
```

Or start python and import the `Experiment` and `InteractivePlot` classes from the `trace_analysis` library located in the `traceanalysis` folder:

```python
from trace_analysis import Experiment
from trace_analysis import InteractivePlot
```

An example script using these classes can be found in the `example.py` file. 

## Updating - advanced

Pull the latest version from the git repository and update the `trace_analysis` environment if necessary:

``` batch
conda env update --file environment.yml --prune
```


