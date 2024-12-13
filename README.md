# SyntheticPatterns

## Paper: Self-organized traveling waves in a synthetic multicellular reaction-diffusion system

# Code organization
- Mathematical models: This directory contains custome python code to simulate the mathematical models used in this study.  FiPy installation is needed to run this code. The detailed instructions are given below.
- Data analysis: This directory contains custome python code used to analyze experimentally generated data. Python libraries such as numpy, scipy, matplotlib, and seaborn are needed to run this code. The detailed instructions are given below.


# Software requirements

## Operating systems
- Windows11 
- macOS (Sonoma 14.6.1)

## Python packages
```
fipy
numpy
scipy
matplotlib
pandas
seaborn
```

# Installation
## Creating a virtual environment with all required python packages
Download the source code as ZIP, unzip and open the directory in the terminal. Use the following commands in the terminal - 
'''
python3 -m venv fipy_env

source fipy_env/bin/activate

pip install -r requirements_fipy_env.txt
'''

## Code running instructions
 - Running the mathemtical model code
After successfully installing the virtual environment, change directory to `MathematicalModels` and run the mathematical model (eg., model1.py) using `python3 model1.py`

Runtime (using Apple M2, 16 GB, macOS:Sonoma 14.6.1): 16 s (model1), 14 s (model2), 2min:11s (model3), 1min:33s (model4)

- Runnning the data analysis code:

- Demo datasets (Expected output and runtime):







