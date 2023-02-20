# Multi-Aperture-Quantum-Imaging
Quantum super-resolution imaging using optical systems with multiple apertures.

The manuscript can be found [here](https://www.overleaf.com/read/wxxqxsrhhqwz).

# Software Versions and Hardware
- Matlab 2021b 
- GPU optimized (if available)
- backwards compatibility to earlier editions of Matlab is possible with slight changes

## Setup
You can clone this repository with,

`git clone https://github.com/I2SL/Multi-Aperture-Quantum-Imaging.git`

Otherwise, click on the green "CODE" button on the top right and download the zip folder containing the project files. Then extract the zipped folder in your desired directory

## File Structure
The parent folder `Multi-Aperture-Quantum-Imaging` contains scripts for running a constellation estimation task (complimentary scripts exist for running batch arrays on a computing cluster).

The `utils/` folder has a collection of functions mainly for computing/visualizing the different modal bases and running the expectation maximization algorithm.

To run a single estimation task, open the `run.m` file and run the script. You can adjust several optional variables including:

```

```


### Future Code Updates:
- Make X,Y coordinate grid types just vectors (main change will be in EM algorithm)
- Remove unused functions from utils folder
- Create a structure for Data collection and merging accross multiple datasets
- Make aperture and image plane coordinate definitions independent (no need for reference unit) 
