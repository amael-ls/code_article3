# Description

This 'readme' explains how to reproduce the analytical result in the article *Development of a spatial Escalator Boxcar Train algorithm for sessile species: application to the boreal–temperate forest ecotone* (Le Squin et. al., 2022, LINK). To get a general description of the program, switch to the branch *master* on Github.

This document is organised as follow:

1. Description of subfolders, listed by alphabetical order
2. Input for x-EBT
3. Output from x-EBT

# List of directories

- **analyses**: contains R scripts to analyse the output of simulations ran by C++
- **code_ebt:** contains the C++ scripts for the x-EBT program in different folders:
	- **alglib**: contains the C++ library alglib
	- **cpp**: contains the Escalator Boxcar Train algorithm (x-EBT) in C++
	- **obj_alglib**: contains the *.o files from the compilation of alglib. Those files are created with the `Makefile`
	- **obj_ebt**: contains the *.o files from the compilation of x-EBT (from cpp/ folder). Those files are created with the `Makefile`

- **code_helpers**: contains scripts to help building the initial condition, the environmental conditions, and the directories to save results
- **run:** contains the results and data associated to this branch (*analytical_case*)

# Execution of x-EBT

## Command lines

### Compilation

I compiled the program with **gcc-10.3** with the `Makefile` from the root directory, using the command line:

```sh
make ebt
```

Note that you must have the TBB library. You might have to update the makefile to your own computer, specifically the directory of the TBB library (variable *libs_dir* in the Makefile. I set-up default values for Linux and Mac which worked for me).

### Execution

I use the following command line to run the program:

```sh
./ebt run/simulationParameters.txt
```

where `run/simulationParameters.txt` is a file stored in the folder `run`, and written as follow (you must respect the space before and after equal sign!):

```bash
#### Files
climate_file = ./run/landscape.txt
species_filenames = Acer_saccharum.txt

#### Paths
species_path = ./run/data/speciesParameters/
initPath = ./run/data/initialCondition/
initFilenamePattern = ic_
summaryFilePath = ./run/results/summary/
popDynFilePath = ./run/results/popDyn/

#### Patterns filenames
summaryFilePattern = su_
popDynFilePattern = pd_

#### Time simulation
t0 = 0
tmax = 100
nIter = 751

#### Saving options
saveOnlyLast = false
freqSave = 1

#### Size population
maxCohorts = 1000

#### Order landscape
rasterOrder_Rlang = true
```

To run, the program requires many input and directories that are already prepared for this analytical example (see the branch *master* for a description on how to prepare the inputs):

- Directories to store the results. For the branch *analytical_case*, they are **run/results/popDyn/Acer_saccharum/** and **run/results/summary/Acer_saccharum/**
- Climate conditions. For the branch *analytical_case*, they are in **run/data/landscape_5x5_acsa/**
- Initial condition. For the branch *analytical_case*, it is in **run/data/initialCondition/Acer_saccharum/**
- Species parameters. For the branch *analytical_case*, they are in **run/data/speciesParameters/**

#### climate_file

The structure of the climate file (named `landscapeData.txt` in the above example) contains the following:

```bash
nRow=5
nCol=5
path=./run/data/landscape_5x5_acsa/
delimiter= = 
filenamePattern=climate_
deltaLon=20.43124
deltaLat=20.43124
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

This argument can contain several names, separated by a comma and a space, but for the analytical example, I set it to:

```bash
species_filenames = Acer_saccharum.txt
```

Each file concerning *Acer saccharum* will be read in the folder `species_path`. There are seven files per species, here for *Acer saccharum*:

1. **Acer_saccharum.txt**: contains the name of the species, its fecundity parameter, and its maximum diameter
2. **Acer_saccharum_M.txt**: contains the mortality parameters
3. **Acer_saccharum_dispersal.txt**: contains the dispersal parameters. Only the parameters specified in `keysToRead` from this file are read. For the analytical example, I used a Dirac kernel (i.e., no dispersal outside of a patch).
4. **Acer_saccharum_scaling_M.txt**: contains the mortality scaling
5. **Acer_saccharum_G.txt**: contains the growth parameters
6. **Acer_saccharum_allometries.txt**: contains the allometries for tree height and crown diameter (computed from dbh)
7. **Acer_saccharum_scaling_G.txt**: contains the growth scaling

#### species_path

See section above.

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



# Analysing the results

After each run, you need to process the results with the script **analyses/analytical_results.R**. Once you have done enough run with different $\Delta t$ (i.e., with different values of $n_{\text{iter}}$, while keeping all the other parameters unchanged) then you can use the script **analyses/plot_diff_vs_deltaT.R** to reproduce the figures 5.a and 5.b in the article. The results of this analysis are already stored in **analyses/diffResults/*.rds** 

