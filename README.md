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

## Running An Estimation Task
To run a single estimation task, open the `run.m` file. An example is already prepared for you to run as is. Otherwise the parameters of the imaging system, the target constellation, and the estimation algorithm are summarized below.

```
% parameters
scene;              % input scene.                         [Nx3 matrix] Columns 1 and 2 are the x and y coordinates of the N sources [units : fractional rl]. Column 3 contains the relative brightness of each source. The sum of values in column 3 should equal 1.
n_pho;              % mean photon number                   [integer]
max_order;          % max modal order                      [integer]
basis;              % basis                                [string] ['Gram-Schmidt','Zernike', 'Direct-Detection']
subap_radius;       % sub-aperture radius                  [double] [units : length]
aper_coords;        % sub-aperture centroid positions      [Mx2 matrix] [units : length] columns 1 and 2 are the kx and ky coordinates respectively of the centroids for M circular apertures
mom_samp;           % sub-aperture sampling density along each axis in k-space [double]
img_samp;           % image-plane samples along each axis  [odd integer] Image plane for source position estimates has dimensions [img_samp,img_samp]
EM_max;             % max number of EM iterations          [integer]
brite_flag;         % brightness estimation flag           [boolean]
visualize;          % visualization flag for seeing output [boolean]
```


### Future Code Updates:
- Make X,Y coordinate grid types just vectors (main change will be in EM algorithm)
- Remove unused functions from utils folder
- Create a structure for Data collection and merging accross multiple datasets
