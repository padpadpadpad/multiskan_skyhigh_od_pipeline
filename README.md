
# MultiSkan SkyHigh OD pipeline

This is a pipeline that makes it easy to process data from the ThermoFisher Multiskan Skyhigh plate readers. It relies heavily on the [**tidyverse**](https://tidyverse.tidyverse.org/) and [**gcplyr**](https://mikeblazanin.github.io/gcplyr/). I recommend reading the walkthroughs within **gcplyr** for an in-depth description of how it works and things you need to think about!

## Help improve this repository

If you have suggestions for how the code (or pipeline more generally) can be improved please open an [Issue](https://github.com/padpadpadpad/multiskan_skyhigh_od_pipeline/issues).

## How to use this repository

- Download this repository and place it in the folder where you are collecting data on a new project.
- Rename `multiskan_skyhigh_od_pipeline.Rproj` to better reflect the active project.
- Put raw data files into `data/raw`.
- Go through the scripts to process your data.
- Make edits to the scripts as necessary to suit your needs.
- At the end of **03_analyse.R** you should have a processed dataset of growth metrics to use in downstream analyses.

### Naming of raw files

How you name raw files is key for analysing data spread across multiple files. Naming in a consistent way allows you to combine data from multiple files, and allows you to describe metadata.

- Use lowercase letters, with underscores separating treatments.
- Start each name with "run**A**" where **A** is the run number. This allows all the raw and processed data through the scripts to be easily referenced.
- For example, if I was analysing plates that were done at different temperatures, I would name the files "run**A**_temp**B**" where **A** is the run number (chronologically from the first you do during an experiment) and **B** is the temperature. For the first run at 40ºC, This would return a file name **run1_temp40.xlsx**.

### Plate layout

To be able to assign meta-data to each plate, you need to go into each raw file and fill in the design template in the **Plate layout** tab. Once this is done, save the file and close it. Currently, we assume blank wells start with an **"X"**, but this can be changed in **03_analyse.R**.

### Where you have to make changes

Each script has a section called "things to change". Specifically:

- **01_process_od.R**: Need to make sure you are reading in the right files and that you set `output` to be the correct name for the output file. 
- **02_visualise_runs.R**: Need to make sure you are reading in the right files. A good rule of thumb is to set `input` to be the same as `output` from **01_process_od.R**.
- **03_analyse.R**: 
    - Need to make sure you are reading in the right files. You may want to split up the column named `file` into the different treatments. This can be done using **tidyr::separate()**. If you do this, you may want to change the variables in `groupings` to reflect this.
    - You can set the method for calculating blanks to either `blank_median` or `well_specific`. It is `blank_median` by default, which takes all the blanks per plate and calculates the median and sd of `raw_absorbance`. If there are any contaminated blanks you should check WHY and remove them. The `well_specific` method calculates the mean `raw_absorbance` from the first three time points of each well, representing a well-specific blank. This is based on guidance from [Atolia (2020) mBio](https://journals.asm.org/doi/10.1128/mbio.01378-20). Using the first few timepoints to subtract was used by [Govers (2024) Cell Systems](https://www.cell.com/cell-systems/fulltext/S2405-4712(23)00331-9). This method relies on cultures not yet reaching their detectable range.
    - Set how often the plate reader was taking measurements (assumed to be ubiquitous across all runs). Sets window width as number of measurements in 1 hour.
    - At ~Line 120 under the comment `# quick plot of data ####`, you may need to change the grouping variables to get a nice overall plot. Specifically might want to change `facet_wrap()` and `interaction(well, run)` to reflect your grouping variables.
    - At ~Line 174 there is a check to see if the number of rows for the output is as expected. If not need to do some troubleshooting!
    - From ~Line 189 to ~Line 290, a for loop creates a pdfs of each plate with the estimated growth metrics placed on top of the raw data. This code currently needs a lot of editing to make sure your grouping variables are correct. It currently assumes you have different runs and temperatures within runs.

## Individual scripts

### The pipeline

- **00_functions.R**: contains helper functions that are used in other scripts.
- **01_process_od.R**: processes raw data files from the MultiSkan SkyHigh plate reader. 
- **02_visualise_runs.R**: Visualises the raw data from a run of a MultiSkan Skyhigh plate readers. Creates a pdf of the raw data in `plots/first_look_plots`. At this stage it might be useful to look at each plot and note down specific well by file combinations that have failed/have poor data to remove from downstream analyses. And also check whether any of your blanks have been contaminated.
- **03_analyse.R**: Calculates growth metrics from the processed datasets and creates visualisations of the processed datasets indicating the position of lag time, maximum growth rate, and maximum density in each well of a plate.
    - Currently used **gcplyr** to calculate max growth rate, minimum density, lag time, time at which max growth rate was estimated, density at which max growth rate was calculated, maximum density, time at which max density was reached, doubling time, and area under the curve.
    - Makes visualisation of the growth metrics on the raw and log scale for corrected OD. Needs a lot of changing based on the treatments and grouping variables you have.
    
### Extras

- **04_tpcs.R**: a script that demonstrates how to fit thermal performance curves to microbial growth rate data.
