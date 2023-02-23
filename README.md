# SCIDATA

## Introduction

This my version of a data analysis tool for the code AENUS-ALCARR. It is a _easy-to-use_ package built entirely in `Python`, and it is meant to work alongside with the _soon to be release_ __scivis__ (a tool for visualizing supernovae data from simulations made withthe same code), but works also standalone. All features are implemented, tested and work fairly well for 1D and 2D simulations, while I did not test them in 3D. However, _most_ of them should work also in 3D.<br/>
Please bear in mind that this package will be updated fairly often, so refer always to the latest version.<br/>
Any tips or suggestion are more than welcome! 

## Requirements

This being a `Python` package, it requires a version of [Python3](https://www.python.org/) installed on your computer. Besides, it makes heavy use of three `Python` libraries
 - [NumPy](https://numpy.org/)
 - [SciPy](https://scipy.org/)
 - [h5Py](https://www.h5py.org/)
 - [Pandas](https://pandas.pydata.org/docs/index.html)
 - [f90nml](https://pypi.org/project/f90nml/)

To install their latest version just run
```
pip install -r requirements.txt
```
inside the `scidata` folder.

## Installation
To install this package just clone 
```
git clone https://github.com/MarcoCusinato/scidata.git
```
or download and unzip it into your Python3 `site-packages` folder, that is usually located in
- Windows:
```
C:/Users/<USER>/AppData/Roaming/Python/Python3X/site-packages
```
- Linux server:
```
/home/<USER>/.local/lib/python3.X/site-packages
```
- Anaconda environment:
```
/home/<USER>/anaconda3/envs/<env_name>/lib/python3.X/site-packages/
or
/home/<USER>/.conda/envs/<env_name>/lib/python3.X/site-packages/
```
where X is the `Python` version for which you would like to install `scidata`.<br/>
Finally, navigate to the `scidata` folder and 
```
pip install requirements.txt
```
Congrats you have successfully installed `scidata`!!!

## Getting started

### Add search paths
First of all let's include some paths to search for simulations. Navigate to 
```
scidata/paths
```
and run `add_paths.py` with one or both of the following options:
- `--add-paths-file`, followed by one or more paths to files containing search paths
- `--add-paths`, followed by one or more search paths.
A path file can be every text file or a `.bz2` or `.gz` archive containing a text file. Moreover, it should look like a list of paths. <br/>
Let's name `foo.txt` a file with a list of path that looks like:
```
/first/path/to/search
/second/path/to/search
...
```
Then call 
```
python3 add_paths.py --add-paths-file /path/to/foo.txt
>Paths to search: [/first/path/to/search, /second/path/to/search]
>Path: /first/path/to/search
>Checking...
>   /subpath/one
>   /subpath/two
>   ...
>Added n simulations
```
This process can take from few seconds to several minutes depending on the number of simulations present in the chosen paths. <br/>
In case you are working on windows, the first time it will ask you where your Linu server is mounted. If you are working with simulations downloaded in your current machine just leave it blank. <br/>
PLEASE REMEMBER: every time you add a simulation to a folder you have to re-run the script. 

### Data Analysis 

Please before starting the data analysis bear in mind that this script will create a folder named `Simulations_analysis_quantities` in your home directory if you are working on Linux or Desktop on Windows to store some quantities that are numerically expensive to calculate.<br/>
To start the analysis, first of all, load the _main_ script of the module
```
from scidata.quantities.quantities import SimulationAnalysis
```
then load a simulation (i. e. `s20-064`)
```
sim = SimulationAnalysis('s20-064')
```
If there are two simulations with different dimension (i. e. one is 1D and the other 3D) it will ask you to specify the requested one. <br/> 
However, if you have not added any simulation path, you can supply it directly as
``` 
sim = SimulationAnalysis('sim_name', simulation_folder_path = '/path/to/sim/')
```
Now you have access to every `SimulationAnalys` method to get the most common quantities present in a simulation.
Let's try, for example, to calculate the PNS radius for all simulation timestep. To do this just run:
```
sim.get_PNS_radius()
> PNS radius file not found. Creating one with default settings.
> If different settings are needed please refer to the "PNS_radius(...)"
```
With default settings this method would not return nothing but will create a file named `PNS_radius.h5` inside the local storage folder in your home or desktop directory. <br/>
Now re-run the same method with some options enabled to get the averaged PNS radius,
```
averaged_pns_radius = sim.get_PNS_radius(PNS_radius = False, indices = False, min_max_average = True,
                                        ret_time = False, ghost_cells = False, tob_corrected = True)
```
By plotting it you will get something like [this plot]()

Similarly you can get other quantities and derive many useful stuff by using other `SimulationAnalys` methods. Please refer to its methods description (when avaiable) to find the ones you need.

## Future developments
Probably what I will need to use at the moment :)
