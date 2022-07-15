# Description

This 'readme' explains how to use the spatial version of the Escalator Boxcar Train algorithm, hereafter x-EBT. This spatial algorithm has been developed for sessile populations only. The theoretical development can be found at ???, with an application on trees.

This document is organised as follow:

1. Description of subfolders, listed by alphabetical order
2. Input for x-EBT
3. Output from x-EBT

# List of directories

- **alglib**: contains the C++ library alglib
- **analyses**: contains R scripts to analyse the output of simulations ran by C++
- **cpp**: contains the Escalator Boxcar Train algorithm (x-EBT) in C++
- **createIC**: to create initial condition for x-EBT. This directory is made by the script `helper.sh`
- **createLandscape**: to create the landscape (climate conditions, patch size). This directory is made by the script `helper.sh`
- **createParams**: to create the species used by x-EBT. This folder contains R scripts and uses the parameters estimated in @LeSquin2020. This directory is made by the script `helper.sh`
- **helpers**: contains scripts to help building the initial condition, the environmental conditions, and the directories to save results
- **obj_alglib**: contains the *.o files from the compilation of alglib. Those files are created with the `Makefile`
- **obj_ebt**: contains the *.o files from the compilation of x-EBT (from cpp/ folder). Those files are created with the `Makefile`

# Execution of x-EBT

## Command lines

### Compilation

I compiled the program with **gcc-10.3** with the `Makefile` from the root directory, using the command line:

```sh
make ebt
```

Note that you must have the TBB library. You might have to update the makefile to your own computer, specifically the directory of the TBB library.

### Execution

I use the following command line to run the program:

```sh
./ebt run/simulationParameters.txt
```

where `run/simulationParameters.txt` is a file stored in the folder `run`, and written as follow (you must respect the space before and after equal sign!):

```bash
#### Files
climate_file = landscapeData.txt
species_filenames = Acer_saccharum.txt

#### Paths
species_path = ../createParams/
initPath = ../createIC/ic_1/
initFilenamePattern = ic_
summaryFilePath = ./summary/
popDynFilePath = ./popDyn/

#### Patterns filenames
summaryFilePattern = su_
popDynFilePattern = pd_

#### Time simulation
t0 = 0
tmax = 100
nIter = 500

#### Saving options
saveOnlyLast = false
freqSave = 1

#### Size population
maxCohorts = 500

#### Order landscape
rasterOrder_Rlang = true
```

The name of this file does not matter, as long as it matches your command line. I now describe the use of each entry of this file. To run, the program requires many input. Some of them can be created automatically, see the section [Creating the initial condition](#creating-the-initial-condition,-the-landscape,-and-the-species-files).

**Before executing the program x-EBT, the following things must be ready:**

- Folder to store the results
- Climate conditions
- Initial condition
- Species parameters

There are helping scripts to do these tasks.

**For the folder**, I recommend using the target `mkdir_sim` from the makefile:

```sh
make mkdir_sim
```

This will run the script `code_helper/mkdir_sim.sh` with the default argument `./run/` for the folder, and `simulationParameters.txt` for the file storing the simulation informations. The target `mkdir_sim` also accepts two optional arguments to change the folder and filename:

```shell
make run_dir=path/to/dir sim_file=whatever.txt mkdir_sim
```

**For the climate conditions**, it is possible to use:

1. **downloadClimateData.R**: Download world climate data from Chelsa
2. **prepareRawClimateData.R**: Prepare the raw climate data (subset to USA and Canada, compute bioclimatic variable 17---mean precipitation driest quarter)
3. **landscapeForCpp.R**: Write the files needed for x-EBT for the landscape (one file per patch). Note that the folder `run/data/` (or whatever was chosen) should be created before that (it should be the case if the target `mkdir_sim` was executed properly).

#### climate_file

The structure of the climate file (named `landscapeData.txt` in the above example) contains the following:

```bash
nRow=5
nCol=5
path=../createLandscape/clim_5x5/
delimiter= = 
filenamePattern=climate_
deltaLon=20
deltaLat=20.10989
distance=euclidean
```

You must respect the spacing here too before and after equal signs!

- `nRow` and `nCol` gives the size of the landscape (i.e., number of cells in the lattice)
- `path` gives where to look for the climate data
- `delimiter` is the delimiter used in the climate files (here it is an equal sign with a space before and after)
- `filenamePattern` specifies how the climate files are named. Here it is `climate_0.txt, ..., climate_24.txt`
- `deltaLon` and `deltaLat` give the spatial discretisation steps
- `distance` gives which distance is used (so far, there are only **euclidean** and **orthodromic** distances)

To be honest, `climate_file` is not the best name now, but it used to be dedicated to climate only!

#### species_filenames

This argument can contain several names, separated by **a comma and a space**. Example:

```bash
species_filenames = Acer_saccharum.txt, Abies_balsamea.txt
```

Then, each file will be read in the folder `species_path`. There are seven files per species. Here is an example for Abies balsamea, for any other species, replace `Abies_balsamea` by any species name from `species_filenames` (e.g., `Acer_saccharum`):

1. **Abies_balsamea.txt**: contains the name of the species, its fecundity parameter, and its maximum diameter
2. **Abies_balsamea_M.txt**: contains the mortality parameters
3. **Abies_balsamea_dispersal.txt**: contains the dispersal parameters. Only the parameters specified in `keysToRead` from this file are read. There are different options for `refKernel_doi` (which is incidentally not necessariliy a DOI):
   - **10.1016/j.jtbi.2005.12.019**: Moorcroft2006
   - **laplace**: Laplace kernel, Cousens2008, p. 82
   - **10.2307/176541**: 2Dt Clark1999, Boisvert-Marsh2022
   - **gaussian**: Gaussian kernel, Cousens2008, p. 82
   - **dirac**: Dirac kernel (i.e., no dispersal outside of a patch)
4. **Abies_balsamea_scaling_M.txt**: contains the mortality scaling
5. **Abies_balsamea_G.txt**: contains the growth parameters
6. **Abies_balsamea_allometries.txt**: contains the allometries for tree height and crown diameter (computed from dbh)
7. **Abies_balsamea_scaling_G.txt**: contains the growth scaling

#### species_path

See section above

#### initPath and initFilenamePattern

Specifies where to load the initial densities for trees (i.e., initial condition of the forest structure for each plot initially colonised). If a plot is originally empty, they there should be no file for that plot. Only filenames starting by the pattern specified in `filenamePattern` will be read. In this example, the filenames start by `ic_` (ic stands for initial condition).

#### summaryFilePath, popDynFilePath, summaryFilePattern, popDynFilePattern

The first two are the paths to save the summary results and population dynamics results, respectively. The program will automatically create subfolders for each species that has been initialised. The last two specifies the beginning of each filename (one file per patch).

#### t0, tmax, nIter

These set the starting time, the ending time, and the number of iterations, respectively. The program computes the time step from them.

#### saveOnlyLast, freqSave

Tells if only the last iteration should be save for the population dynamics, or if it should be saved with a given frequency `freqSave`. Note that if `saveOnlyLast` is true, then `freqSave` is disregarded.

#### maxCohorts

The maximal number of cohorts. If during the simulation the number of cohorts overcome `maxCohorts`, then the program tries to merge similar cohorts. If it is impossible to merge them, then the program throw a segmentation fault.

#### rasterOrder_Rlang

It specifies how the cells of the landscape are organised. If `rasterOrder_Rlang` is true, then the default raster order from R language is used. For instance a 3 row x 5 columns landscape would be

```R
0  1  2  3  4
5  6  7  8  9
10 11 12 13 14
```

### Creating the initial condition, the landscape, and the species files

The first thing is to create the file `simulationParameters.txt` (or whatever name you choose, see the content of this file in the section [Execution](#execution)) and then run the script `helpers/helpher.sh` that will create automatically the necessary directories. Then, you can run the following R scripts to help you with the climate data and the initial condition:

1. **landscapeForCpp.R**: lala
2. 
