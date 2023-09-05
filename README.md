# GANGLIA
Generation of Automatic Neuron Graph-Like Interconnected Arrangements (GANGLIA) for designing customized neuron circuit patterns

[![DOI](https://zenodo.org/badge/622026352.svg)](https://zenodo.org/badge/latestdoi/622026352)

Ashlee S. Liao, Yongjie Jessica Zhang, Victoria A. Webster-Wood
Liao, A.S., Zhang, Y.J., Webster-Wood, V.A. (2023). GANGLIA: A Tool for Designing Customized Neuron Circuit Patterns. In: Meder, F., Hunt, A., Margheri, L., Mura, A., Mazzolai, B. (eds) Biomimetic and Biohybrid Systems. Living Machines 2023. Lecture Notes in Computer Science(), vol 14158. Springer, Cham. https://doi.org/10.1007/978-3-031-39504-8_14

## Table of Contents
* [Programs Used for Developing and Running GANGLIA](#programs-used-for-developing-and-running-ganglia)
  * [Required Packages](#required-packages)
* [Reproducing the Outputs using GANGLIA](#reproducing-the-outputs-using-ganglia)

## Programs Used for Developing and Running GANGLIA
- Python 3.9.7
- Anaconda 3 (https://www.anaconda.com/)
- Spyder 5.4.2 (Anaconda distribution or https://www.spyder-ide.org/)
```
conda update anaconda
conda install spyder=5.4.2
```

### Required Packages 
(Outside of Default Spyder/Anaconda Installation Packages)
- igraph (0.9.11) (https://python.igraph.org/en/stable/) (0.9.11-py39h4a3397e_0)
```
conda install -c conda-forge python-igraph 
```
- ezdxf (1.0.2) (might not be necessary due to cadquery installation) (https://ezdxf.readthedocs.io/en/stable/setup.html)
```
pip3 install ezdxf
```
- cadquery (2.2.0) (currently installed cadquery-master-py3.9) (https://cadquery.readthedocs.io/en/latest/installation.html)
```
 conda install -c cadquery -c conda-forge cadquery=master
```
- jupyter-cadquery (3.5.2) (https://github.com/bernhard-42/jupyter-cadquery)
  - NOTE: Installing jupyter-cadquery without downgrading nbconvert will ultimately throw an error due to nbformat being a wrong version (5.7.0) when it is looking for 5.4.0
```
conda install nbconvert==6.5.4
pip install jupyter-cadquery==3.5.2 cadquery-massembly==1.0.0 matplotlib
(Windows) conda install pywin32
```
- cairosvg (https://cairosvg.org/)
```
conda install -c conda-forge cairosvg
```
- NOTE: When installing on a new Windows 11 machine - it seems that matplotlib threw an error for importing the library - it is related to Pillow. This was resolved by upgrading pillow (https://pillow.readthedocs.io/en/stable/)
```
conda upgrade Pillow
```

## Reproducing the Outputs using GANGLIA
1. Edit line 117 (exportDir) to the appropriate export dirctory
2. Comment in/out the appropriate connectivity network based on desired output. (Lines 151 - 303)
    - For the half-center oscillator, lines 167-170
    - For the user-designed 9-cell network, lines 232-243
    - For the rat single limb joint control network, lines 185-225
    - For the boolean network of _Aplysia_ feeding, lines 257-299
3. Run GANGLIA.py

## (Work in Progress) GANGLIA - Interactive Python Notebook (GANGLIA.ipynb)
- Reproducibility
  - Python 3.11
  - Use a conda virtual environment to set up the same packages
  - miniconda (conda 23.7.3)
  - Windows 11
  ```
  conda env create --file condaVENV_GANGLIA_v01.yaml --name virtual_environment_name
  conda activate virtual_environment_name
  ```
  - manual creation (not recommended)
  ```
  conda create --name condaVENV_GANGLIA_v01
  conda activate condaVENV_GANGLIA_v01
  conda install matplotlib
  conda install -c conda-forge python-igraph
  conda install -c conda-forge ezdxf
  conda install -c conda-forge cadquery
  conda install -c conda-forge cairosvg
  conda install ipykernel
  ```

- Important ongoing changes - GANGLIA.ipynb
  - Removed the need for jupyter-cadquery as it is unnecessary if already in a Python notebook
  - split the code into individually runnable sections

