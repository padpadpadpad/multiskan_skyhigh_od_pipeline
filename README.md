
# MultiSkan SkyHigh OD pipeline

This is a pipeline that makes it easy to process data from the ThermoFisher Multiskan Skyhigh plate readers. It relies heavily on the **tidyverse** and **gcplyr**.

## Help improve this repository

If you have suggestions for how the code (or pipeline more generally) can be improved please open an [Issue](https://github.com/padpadpadpad/multiskan_skyhigh_od_pipeline/issues).

## How to use this repository

- Download this repository and place it in the folder where you are collecting data on a new project.
- Rename `multiskan_skyhigh_od_pipeline.Rproj` to better reflect the active project.
- Put raw data files into `data/raw`.
- Go through the scripts to process your data.
- Make edits to the scripts as necessary to suit your needs.
- At the end of **03_analyse.R** you should have a processed dataset of growth metrics to use in downstream analyses.

#### Naming of raw files

How you name raw files is key for analysing data spread across multiple files. Naming in a consistent way allows you to combine data from multiple files, and allows you to describe metadata.

- Use lowercase letters, with underscores separating treatments.
- Start each name with "run**A**" where **A** is the run number. This allows all the raw and processed data through the scripts to be easily referenced.
- For example, if I was analysing plates that were done at different temperatures, I would name the files "run**A**_temp**B**" where **A** is the run number and **B** is the temperature. For the first run at 40ºC, This would return a file name **run1_temp40.xlsx**.

#### Plate layout

To be able to assign meta-data to each plate, you need to go into each raw file and fill in the design template in the **Plate layout** tab. Once this is done, save the file and close it.

#### Where you have to make changes

Each script has a section called "things to change". Specifically:

- **01_process_od.R**: Need to make sure you are reading in the right files and that you set `output` to be the correct name for the output file. 
- **02_visualise_runs.R**: Need to make sure you are reading in the right files. A good rule of thumb is to set `input` to be the same as `output` from **01_process_od.R**.

## Individual scripts

### The pipeline

- **00_functions.R**: contains helper functions that are used in other scripts.
- **01_process_od.R**: processes raw data files from the MultiSkan SkyHigh plate reader. 
- **02_visualise_runs.R**: Visualises the raw data from a run of a MultiSkan Skyhigh plate readers. Creates a pdf of the raw data in `plots/first_look_plots`. At this stage it might be useful to look at each plot and note down specific well by file combinations that have failed/have poor data to remove from downstream analyses.
- **03_analyse.R**

### Extras

- **04_tpcs.R**: a script that demonstrates how to fit thermal performance curves to microbial growth rate data.
