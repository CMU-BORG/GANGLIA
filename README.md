# GANGLIA
Generation of Automatic Neuron Graph-Like Interconnected Arrangements (GANGLIA) for designing customized neuron circuit patterns

Ashlee S. Liao, Yongjie Jessica Zhang, Victoria A. Webster-Wood
GANGLIA: A tool for designing customized neuron circuit patterns
Submitted to: Living Machines 2023

## Programs Used for Developing the Tool
- Python 3.9.7
- Anaconda 3 (https://www.anaconda.com/)
- Spyder 5.4.2 (Anaconda distribution or https://www.spyder-ide.org/)
  - conda update anaconda
  - conda install spyder=5.4.2

### Required Packages (Outside of Default Spyder/Anaconda Installation Packages)
- igraph (0.9.11) (https://python.igraph.org/en/stable/)
  - conda install -c conda-forge python-igraph (0.9.11-py39h4a3397e_0)
- ezdxf (1.0.2) (might not be necessary due to cadquery installation)
  - pip3 install ezdxf (https://ezdxf.readthedocs.io/en/stable/setup.html)
- cadquery (2.2.0) (https://cadquery.readthedocs.io/en/latest/installation.html)
  - conda install -c cadquery -c conda-forge cadquery=master (currently installed cadquery-master-py3.9)
- jupyter-cadquery (3.5.2) (https://github.com/bernhard-42/jupyter-cadquery)
  - conda install nbconvert==6.5.4
    - NOTE: Installing jupyter-cadquery without downgrading nbconvert will ultimately throw an error due to nbformat being a wrong version (5.7.0) when it is looking for 5.4.0
  - pip install jupyter-cadquery==3.5.2 cadquery-massembly==1.0.0 matplotlib
  - (Windows) conda install pywin32 (recommended by the jupyter-cadquery github install           instructions to do after installing jupyter-cadquery)
- cairosvg (https://cairosvg.org/)
  - conda install -c conda-forge cairosvg
- NOTE: When installing on a new Windows 11 machine - it seems that matplotlib threw an error for importing the library - it is related to Pillow. This was resolved by upgrading pillow (https://pillow.readthedocs.io/en/stable/)
  - conda upgrade Pillow

## Reproducing the Outputs using GANGLIA.py
1. Edit line 117 (exportDir) to the appropriate export dirctory
2. Comment in/out the appropriate connectivity network based on desired output. (Lines 151 - 303)
    - For the half-center oscillator, lines 167-170
    - For the user-designed 9-cell network, lines 232-243
    - For the rat single limb joint control network, lines 185-225
    - For the boolean network of _Aplysia_ feeding, lines 257-299
3. Run GANGLIA.py

