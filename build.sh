#!/bin/bash

### Build script for TVGuide-Aniso - 
### This build script will download the requisite Python modules in order to appropriately setup TVGuide-Aniso 
### This task is accomplished by first checking to see if there is an appropriate local "tvguide" conda/mamba environment available to use
### If there is not, the script will download Mamba and create a local Mamba "tvguide" environment

### Local Mamba Installation - should be located in root of tvguide-aniso



export MINIFORGE_LOCATION=$(dirname ${BASH_SOURCE[0]})

### Checking for "tvguide" environment - if it exists, skip
VENV_NAME=$(dirname ${BASH_SOURCE[0]})/tvguide
if [ -d $VENV_NAME ]; then
    echo "$VENV_NAME already exists. Skipping environment creation"
    exit 0 
fi 


### Checking for mamba/conda commands - if the command exists, user has own mamba/conda installation, and can use that instead
if [[ -z `command -v conda` || -z `command -v mamba` ]]; then
    echo "Conda/Mamba commands not found - initiating installation"
    wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
    MINIFORGE_PREFIX="$MINIFORGE_LOCATION/.miniforge3"
    bash Miniforge3-$(uname)-$(uname -m).sh -b -p $MINIFORGE_PREFIX
    rm Miniforge3-$(uname)-$(uname -m).sh
    CONDA_PROFILE="$MINIFORGE_PREFIX/etc/profile.d/conda.sh"
else
    CONDA_PROFILE=$CONDA_PREFIX/etc/profile.d/conda.sh
    echo "Conda/Mamba exists already - continuing"
fi



### We source the appropriate conda profile from either the mamba installation that was just setup, or the users own conda profile
### We then begin environment creation, and install the necessary 

source $CONDA_PROFILE
echo "Using $CONDA_EXE to create $VENV_NAME ."
conda create --prefix $VENV_NAME python=3.11 --yes
conda activate $PWD/$VENV_NAME
mamba install pandas=1.5.3 numpy scipy dash=4.4.0 dash-bootstrap-components gunicorn matplotlib iteration_utilities --yes
exit 0 
