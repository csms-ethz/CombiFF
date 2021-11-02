# Installation Guide

## Introduction

To compile the C++ programs of CombiFF, you'll need CMake (https://cmake.org/). To use the python scripts, we recommend using Anaconda (https://www.anaconda.com/)

## Installation

```
# clone the repository
git clone https://github.com/csms-ethz/CombiFF.git
cd CombiFF

# create and activate the conda environment
conda env create -f dev/conda_envs/combiff.yml 
conda activate combiff

# compile
python scr/compile.py

# alternative compilation without python script
mkdir build
cd build
cmake ..
make
```

After installation, the executables for the different programs are located in the `bin` directory
