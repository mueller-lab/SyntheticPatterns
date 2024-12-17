# SyntheticPatterns

## Introduction
This repository provides custom code used in our paper: "Self-organized traveling waves in a synthetic multicellular reaction-diffusion system"

In this study, we developed a mathematical modeling framwork for self-organized pattern formation via reaction-diffusion mechanism. We combined mathematical modeling with synthetic biology approaches to create a synthetic multicellular system capable of generating self-organized travelling waves. 

### Code organization
- <a href="https://github.com/mueller-lab/SyntheticPatterns/tree/main/MathematicalModels"><b>Mathematical models </b></a>: This directory contains custom Python code to simulate the mathematical models used in this study. FiPy installation is needed to run this code. Detailed instructions are given below.
- <a href="https://github.com/mueller-lab/SyntheticPatterns/tree/main/DataAnalysis"><b>Data analysis </b></a>: This directory contains custom Python code used to analyze experimentally generated data. Python libraries such as numpy, scipy, matplotlib, and seaborn are needed to run this code. Detailed instructions are given below.


## Software requirements

### Operating systems
- Windows11 Pro 
- macOS (Sonoma 14.6.1)

### Python packages
```
fipy
numpy
scipy
matplotlib
pandas
seaborn
```

## Installation
### Creating a virtual environment with all required Python packages
Download the source code as ZIP, unzip and open the directory in the terminal. Use the following commands in the terminal - 

On MacOS
`python3 -m venv synPat_env`

`source synPat_env/bin/activate`

`pip install -r requirements_synPat_env.txt`

On windows OS
`python3 -m venv synPat_env`

`. synPat_env/Scripts/activate`

`pip install -r requirements_synPat_env.txt`

## Code running instructions
### Running the mathematical model code
- After successfully installing the virtual environment, change directory to `MathematicalModels` and run the mathematical model (eg., model1.py) using ... 

`cd MathematicalModels`

`python3 model1.py`

- Runtime (using Apple M2, 16 GB, macOS:Sonoma 14.6.1): 16s (model1), 14s (model2), 2min:11s (model3), 1min:33s (model4), 28s (model5),12min:36s (model5_2D).

### Running the data analysis code

- <a href="https://github.com/mueller-lab/SyntheticPatterns/tree/main/DataAnalysis/ImageAnalysisMacros"><b>ImageAnalysisMacros </b></a>: The custom ImageJ macro uses the BaSiC plugin for shading correction. Please install FiJi (ImageJ 2.9.0/1.53t) and the plugin as detailed on - https://github.com/marrlab/BaSiC. Then, the macro can be tested with the provided demo dataset. To do so, unzip the demo data file. Run the macro in ImageJ and select the correct inputs when prompted - select directories - `flatField`, `Input`, and `Results` in that order.

- <a href="https://github.com/mueller-lab/SyntheticPatterns/tree/main/DataAnalysis/PythonScripts"><b>PythonScripts </b></a>: Custom Python scripts were used in this study for analysis of 1. diffusion coefficient estimates, 2. fluorescence plate reader assay data, 3. qPCR data, and 4. traveling wave data. The code as well as example input data are provided. To test the code, simply run the respective Python scripts in the directory - for example, run `python3 kymograph_plotting.py` in terminal (with synPat_env active) for testing `kymograph_plotting.py`.



