# Project Name

Does the Evolution of Ontogenetic Niche Shifts Favor Species Coexistence? An empirical test in Trinidadian streams. 
By Jaime M. Anaya-Rojas, Ronald D. Bassar, Blake Matthews, Joshua F. Goldberg, Leighton King, David Reznick, and Joseph Travis.
Acepted in the Journal Of Animal Ecology


## Overview

This repository contains the R code for analyzing data on the ontogenetic niche shifts of two coexisting fish species, Trinidadian guppies (Poecilia reticulata) and killifish (Rivulus hartii), in three different natural community types and experimental communities. The study aims to test the prediction that competing species can reduce competitive interactions and facilitate coexistence through ontogenetic changes in the trophic niche of one or more of the interacting species. 

The project utilized stable isotopes (δ15N and δ13C) to examine fish from different locations representing the three natural community types, and experimental killifish and guppy communities. The results provide comparative evidence for ontogenetic niche shifts in contributing to species coexistence and comparative and experimental evidence for evolutionary or plastic changes in ontogenetic niche shifts following the formation of new communities.

The R code in `Main_script.R` provides a reproducible and documented analysis of the data and can be used as a template for analyzing similar data sets.

## Data

The `data/data.csv` file contains information on the ontogenetic niche shifts of the guppies and killifish in different natural and experimental communities. The data includes the following variables:

- `Location`: the location where the fish sample was collected (a factor variable with three levels: "predator", "no predator", and "killifish only").
- `Species`: the species of fish (a factor variable with two levels: "guppy" and "killifish").
- `Size`: the size of the fish in millimeters.
- `δ15N`: the stable isotope metric for trophic position.
- `δ13C`: the stable isotope metric for carbon source.

The data set has 1,052 observations and 5 variables. The data set includes missing values that need to be removed or imputed before the analysis.

The data and code for this project can also be found in the Dryad repository [here](https://datadryad.org/stash/dataset/XXXXX).

This project was used to conduct a scientific study titled "Does the Evolution of Ontogenetic Niche Shifts Favor Species Coexistence? An empirical test in Trinidadian streams", which was published in [Insert journal name here].

## Requirements

The R packages required for running the code in `R/Main_script.R` are listed in the `requirements.R` file. You can install the packages by running the following command in your R console:

'''R
install.packages(c("MuMIn","partR2","stringr","plyr","lme4","lmerTest","Hmisc"))
'''

## Usage

To use the code in `R/Main_script.R`, follow these steps:

1. Clone this repository to your local machine.
2. Install the required R packages by running the command in the `requirements.R` file.
3. Open the `R/Main_script.R` file in RStudio or your preferred R editor.
4. Set the working directory to the root directory of the cloned repository.
5. Run the script to perform the analysis on the data in `data/data.csv`.

## Results

After running the script, the results of the analysis will be printed to the R console and saved to a file in the `results/` directory. The R code provides a detailed analysis of the data and includes statistical tests to determine whether ontogenetic niche shifts contribute to species coexistence.

## Contributing

If you would like to contribute to this project, please fork the repository and create a pull request.

