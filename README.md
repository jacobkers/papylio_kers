# Papilio

Software for single-molecule trace analysis

## Prerequisites

### Anaconda

Install Anaconda 64-bit from https://www.anaconda.com/distribution/  
    IMPORTANT: For Windows under advanced options tick the box "Add Anaconda to my PATH enivronment variable"

*Note: Advanced users may prefer the more lightweight Miniconda without any packages preinstalled and without a graphical user interface, i.e. you can only use it through the command line. Miniconda is available from https://docs.conda.io/en/latest/miniconda.html*

### Git

GitHub desktop provides an easy way to use git.
Download and install GitHub Desktop from https://desktop.github.com.

*Note: Advanced users may prefer using git via the command line. Git may already be preinstalled, but otherwise you can get it at https://git-scm.com/downloads.*

## Installation

### Clone the repository to your computer
- Open GitHub Desktop (Skip any login screens)
- Click File > Clone repository
- Go to the URL tab, paste the HTTPS URL of the repository and select the local path where you want the code to be located.
- Click clone

### Add a Conda environment
- Open Anaconda-Navigator
- Go to Environments (located on the left)
- Click Import (located on the bottom) and import the environment.yml file in the code. Please don't change the default environment name.

## Using the GUI
Start the user interface by double clicking the file `Start GUI.cmd` for Windows or `Start GUI.command` for Mac

---

## Installation - advanced

Clone the git repository:

```batch
clone https://gitlab.tudelft.nl/ivoseverins/traceanalysis.git
```

Create a conda environment named `papilio` based on the `environment.yml` file. (Make sure you are in the main folder)

```batch
conda env create --file environment.yml
```



