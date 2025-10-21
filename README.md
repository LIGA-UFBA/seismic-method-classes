# Seismic Method Classes

This repository contains the Devito notebooks used in the "GEO208 - Seismic Methods" course at Federal University of Bahia.

## How to use

To run the notebooks, it is necessary to have the [GeoSymCodes' Devito](https://github.com/GeoSymCodes/devito) and the Jupyter Notebook packages installed. Devito only works on Linux machines, so if you want to install Devito in a Windows machine, you need to use it through WSL or a virtual machine.

### Devito installation guide

>It is recommended to install Devito under a Conda or Virtualenv environment. This guide assumes the user will use a Conda environment.

1. Create a conda environment with python>=3.10 

        > conda create --name devito python=3.10
        > conda activate devito

2. Clone the devito repository to your machine

        > git clone https://github.com/GeoSymCodes/devito.git

3. Enter the newly created folder through your terminal and install the following packages

        > cd devito
        > pip install "codepy<2025.1"
        > pip install -e .[extras,tests]
        > pip install -3 .
        > pip install notebook

With, Devito and Jupyter Notebook should be installed in your computer under the environment "devito". You can test if everything worked by going back to this folder through your terminal, launching Jupyter Notebook by executing

    jupyter notebook

in this folder. And running of of the notebooks.