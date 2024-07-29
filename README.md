# velo_calibration_scripts

## Description
The code in this repository contains the functions and framework necesary for the analysis of data taken from exposure of the VeLo ASICS to a radiation source. 
## Installation
This project has been tested using Python 3.9.16 with the dependencies below specified. To use this project, follow these steps:

1. Clone the repository to your local machine:

2. Install the required dependencies:
!pip install numpy matplotlib pandas scipy glob os csv inspect

3. Navigate to the project directory:

4. Run the run.py script to start the data analysis:
`python run.py`

5. Explore the different modules and functions available in the codebase to customize your data analysis. Refer to the inline comments in the code for detailed explanations of each part.

## Usage
First the dataset specifics corresponding to the data to be analysed need to be added to the code. The flux equation to be used with the dataset can also be specified if needed. The code can then take this information and use it to carry the analysis and plotting of the various figures.

During execution of the code, lines are printed inside the terminal to provide guidance on what is being analysed, as well as providing the location of the figures and files created. Log files are created for the datasets analysed containing information on some relevant parameters found throughout the analysis, together with a `results_collection.csv` file which collects the results from various datasets to then be processed and plotted all together if desired.


### Setting up Data Information for Analysis

1. Open the `utils/datainfo.py` file in your project directory.

2. Create a `type_data` if condition to have your information saved on. I will use `type_data = calibration` as an example. 

3. Update the following variables in the `datainfo.py` file with the appropriate data paths and parameters for your analysis:
   - `datafolder_l`: List of data folders containing the data to be analyzed.
   - `FPTfolder_l`: List of folders for storing results. (e.g.`2*['/data/bfys/apuicerc/N037_new/resultscalibration/']`)
   - `FPTfilename_l`: List of filenames for the data.
   - `ST_l`: List of time values.
   - `nacq_l`: List of acquisition numbers.
   - `fromThr_l`: List of start threshold values.
   - `toThr_l`: List of end threshold values.
   - `stepThr_l`: List of step threshold values.
   - `label_l`: List of labels for the data.
   - `bad_cut_l`: List of bad cut values.
   - `condition_l`: List of conditions for data filtering (e.g. [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']).
   - `p0`: List of initial parameter values.
   - `p0_baseline`: List of initial baseline parameter values.
   - `baseline_opt`: Boolean value indicating whether baseline is included before or after fitting of data. *
   - `separation_opt`: Boolean value indicating whether separation in groups of pixels is enabled. *
   - `predict_opt_l`: List of values for prediction options.
   - `output_opt_l`: List of values for output options.

    *(Not sure if it is used by the code or if it takes this information from input in terminal when task is specified or when prompted by code)

4. Save the `datainfo.py` file after making the necessary changes.

5. Run the main script to start the data analysis. After this, you can access the data thoughout the analysis by using `-dt calibration` (e.g. `type_data == calibration`) when running the script. The script will create the `yourpath/resultscalibration/` and `/resultscalibration/plots/` folder structures for all the analysis files and plots.

See `type_data = calibration` in `utils/datainfo.py` for reference on default inputs.

### Data Processing and Analysis
The data processing and analysis can now be performed with the data information saved (e.g. `calibration`). Now `python run.py analysis` can be used. Using `-h` gives you the various options but the following usage should be followed to generate all the necessary files (You will be prompted if you want to run the process again if the files are found already within the dataset directory):

1. First the Flux files need to be obtained. This can be done by running:
`python run.py analysis flux -dt calibration`
If `-t (e.g. warm, cold)` is not provided, it goes through all temperatures provided in `datainfo.py`

2. Second, the data can be fitted to the flux equation specified in `analyse.py`, thus obtaining a categorisation of pixels depending on the quality of fit, as well as matrices containing the individual pixel fitted parameters involved in the flux equation. This can be done using:
`python run.py analysis fitter -dt calibration` 

You will be prompted if you want to consider baseline shift before fitting or afterwards. (Both are equivalent mathematically)

3. Third and last, the target and gain K matrices containing as well as the ASIC and ASICgood values can be found running:
`python run.py analysis matrices -dt calibration`

You will be prompted if you want to consider baseline shift while fitting ASIC. This also saves all the information regarding mean, ASIC, ASICgood and its uncertainties for target and K, on `results_collection.csv` file which saves all these values for later comparison with other datasets (such as e.g. `calibration`).

4. After running the three options of the analysis section of the code, the analysis of data is done completely.

If a separation in groups of pixels is wanted for the analysis, specify so when writting the task in terminal:

`python run.py analysis flux -dt calibration -sep even`

### Plotting
Now that all files required for the analysis have been found, the `plot` section of this code can be used for plotting of data to provide insights on the data. The code can create:

`-plot_heatmap, -hm `   make a 2D heatmap of parameter for pixels
`-plot_results_collection, -rc `
                        make a plot showing target and K[e/DAC] histograms for warm and cold scans.
`-distribution, -d  `   make a histogram distribution of parameter. Can be specified with --param option.
`-flux_over_thresholdscan, -flux_vs_thr `
                        make a flux [Hits/s] vs threshold scan [DAC] plot (for individual pixels and ASIC). Can calculate residuals, later requested while code is executed.
`-opt [{None,discrepancy,QQ,summary}], -o [{None,discrepancy,QQ,summary}]`
                        addition option to plot, default - None. Discrepancy compares in number of std to a specific value. Summary creates plot showing results for K[e-/DAC].
`--param [{E0,f,s,target,Fittype,totalflux,K}], -p [{E0,f,s,target,Fittype,totalflux,K}] `
                        select parameter to plot, default - None. Analyses all parameters within the dataset. Can be used with: -plot_heatmap, -distribution
`--all, -a `            make all available plots for certain type of scan. (this doesnt work taken from other code)
`--separation [{even,odd,rows}], -sep [{even,odd,rows}] `
                        option for separation in even, odd columns and 16th rows, default - None.

Example of usage:

`python run.py plot -distribution -dt calibration -p E0 -opt discrepancy`

`python run.py plot -d -dt calibration`

`python run.py plot -plot_heatmap -dt calibration -p target`

`python run.py plot -flux_vs_thr -dt calibration -t cold`

`python run.py plot -rc -opt summary`


