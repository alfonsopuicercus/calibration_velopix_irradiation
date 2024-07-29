import os
import glob
import numpy as np
import csv

class FileOrganizer():
    def __init__(self, path = '', prefix = ''):
        self.path = path
        self.prefix = prefix
        self.check_and_create_path(self.path)
        self.csv_files = self.get_csv_files(prefix = '', path = self.path)

    def get_csv_files(self, prefix, path):
        """Get a list of CSV files in a given path starting with a prefix."""
        csv_files = []
        for file in os.listdir(path):
            if file.startswith(prefix) and file.endswith('.csv'):
                csv_files.append(os.path.join(path, file))

        return csv_files
    
    def open_csv(self, prefix, label = '', skiprows = 0, ravel = False, path = ''):
        csv_files = self.get_csv_files(prefix = prefix, path = path)
        
        if not csv_files:
            print('No files found with prefix {}'.format(prefix))
            return None

        data_list = []
        for file in csv_files:
            if label and not label in os.path.basename(file):
          
                continue

            print(f'Reading {file}...')
            data = np.loadtxt(file, skiprows = skiprows, delimiter = ',')
            if ravel:
                #print('Raveling data')
                data = data.ravel()

            data_list.append(data)

        return data_list
    

    def is_file_found(self, prefix, label = ''):
        """Checks if a file is in a list of files."""
        csv_files = self.get_csv_files(prefix = '', path = self.path)
        output = False
        
        #print('csv', csv_files)
        for file in csv_files:
            if label in os.path.basename(file) and os.path.basename(file).startswith(prefix):
                output = True

                break

        return output
    
    
    def check_and_create_path(self, path):
        '''
        Checks if a path exists, if not it creates the path.

        Args:
            path (str): Path to be checked and created.
        '''
        if not os.path.exists(path):
            os.makedirs(path)
            print('Path {} created'.format(path))
        
    

    def create_directory_structure(self,path, list_of_params = None):
        '''
        Creates diretory structure required by the other functions to store data and plots.

        Args:
            base_folder (str): Folder where results are going to be saved after data processing and analysis.
        '''

        plots_folder = os.path.join(path, 'plots')
        heat_map_plots_folder = os.path.join(plots_folder, 'heat_map_plots')
        ind_pixel_plots_folder = os.path.join(plots_folder, 'ind_pixel_plots')
        
        folders = ['target', 'counts', 'Fittype', 'totalflux'] #base subfolders to create
        if list_of_params != None: # if the list of parameters involved in the flux equation chose for the dataset is provided then 
            for parameter in list_of_params:
                folders.append(parameter)

        

        dir_exist = True
        # Create 'plots' directory
        if not os.path.exists(plots_folder):
            dir_exist = False
            os.makedirs(plots_folder)

        # Create 'heat_map_plots' directory
        if not os.path.exists(heat_map_plots_folder):
            dir_exist = False
            os.makedirs(heat_map_plots_folder)

        # Create 'ind_pixel_plots' directory
        if not os.path.exists(ind_pixel_plots_folder):
            dir_exist = False
            os.makedirs(ind_pixel_plots_folder)


        # Create subdirectories inside 'heat_map_plots' directory
        for folder in folders:
            folder_path = os.path.join(heat_map_plots_folder, folder)
            if not os.path.exists(folder_path):
                dir_exist = False
                os.makedirs(folder_path)
        
        if dir_exist == False:
            print('Directory and subfolders created...')
        
    
    def results_csv_file(self, savepath, ST,label, baseline, targetmean=0.0, utargetmean=0.0, targetASIC=0.0, utargetASIC=0.0, targetASICgood=0.0, utargetASICgood=0.0, Kmean=0.0, Kstd = 0.0, uKmean=0.0, KASIC=0.0, uKASIC=0.0,KASICgood=0.0, uKASICgood=0.0, path='/data/bfys/apuicerc/N037_new/'):
        '''
        Creates .csv file to store values of data found throughtout the analysis.
    
        Args:
            savepath (str): Name given to the results file. Varies depending on how data is analysed and thus useful to include to differentiate betweem approaches.
            ST (str): Shutter Time.
            label (str): Refers to temperature of ASIC while data taking.
            baseline (bool): True or False depending on inclusion of baseline or not.
            targetmean (float, optional): Mean of target. Defaults to 0.0.
            utargetmean (float, optional): Uncertainty of targetmean. Defaults to 0.0.
            targetASIC (float, optional): Target ASIC. Defaults to 0.0.
            utargetASIC (float, optional): Uncertainty of target ASIC. Defaults to 0.0.
            targetASICgood (float, optional): Target ASIC considering only good pixels. Defaults to 0.0.
            utargetASICgood (float, optional): Uncertainty of target ASIC considering only good pixels. Defaults to 0.0.
            Kmean (float, optional): Mean of K[e-/DAC] factor. Defaults to 0.0.
            uKmean (float, optional): Uncertainty of mean of K[e-/DAC] factor. Defaults to 0.0.
            KASIC (float, optional): K[e-/DAC] ASIC factor. Defaults to 0.0.
            uKASIC (float, optional): Uncertainty of K[e-/DAC] ASIC factor. Defaults to 0.0.
            KASICgood (float, optional): K[e-/DAC] ASIC factor considering only good pixels. Defaults to 0.0.
            uKASICgood (float, optional): Uncertainty of K[e-/DAC] ASIC factor considering only good pixels. Defaults to 0.0.
            path (str, optional): Main path where both data and results is stored. Defaults to '/data/bfys/apuicerc/N037_new/'.
        
        Returns:
            Creates or updates results_collection.csv file with:
            savepath, label, ST, baseline, targetmean, utargetmean, targetASIC, utargetASIC, targetASICgood, utargetASICgood, Kmean, uKmean, KASIC, uKASIC
        
        '''

        # If path is not provided, use current directory
        baseline = str(baseline)
        
        if not path:
            path = "./"
        
        # Ensure path ends with '/'
        if not path.endswith("/"):
            path += "/"
        
        csv_file_path = path + "results_collection.csv"
        savepath_results = savepath.split('/')[5]
        savepath_other = savepath.split('/')[6]
    
        if savepath_other == 'even' or savepath_other == 'odd' or savepath_other == 'rows':
            savepath = savepath_results +'/'+ savepath_other
        else:
            savepath = savepath_results
        
        
    
        # Check if CSV file exists
        rows = []
        existing_row = None
        if os.path.exists(csv_file_path):
            # Read existing data
            with open(csv_file_path, mode='r') as file:
                reader = csv.reader(file)
                rows = list(reader)

            # Find existing row if it exists
            for row in rows:
                if row and row[0] == savepath and row[1] == label and row[2] == ST and row[3] == baseline:
                    row_before_change = row
                    print('Row before change')
                    print(row_before_change)
                    
                    # Only update non-zero values
                    row[1] = label
                    row[2] = ST
                    row[3] = baseline
                    row[4] = targetmean if targetmean != 0 else row[4]
                    row[5] = utargetmean if utargetmean != 0 else row[5]
                    row[6] = targetASIC if targetASIC != 0 else row[6]
                    row[7] = utargetASIC if utargetASIC != 0 else row[7]
                    row[8] = targetASICgood if targetASICgood != 0 else row[8]
                    row[9] = utargetASICgood if utargetASICgood != 0 else row[9]
                    row[10] = Kmean if Kmean != 0 else row[10]
                    row[11] = Kstd if Kstd != 0 else row[11]
                    row[12] = uKmean if uKmean != 0 else row[12]
                    row[13] = KASIC if KASIC != 0 else row[13]
                    row[14] = uKASIC if uKASIC != 0 else row[14]
                    row[15] = KASICgood if KASICgood != 0 else row[15]
                    row[16] = uKASICgood if uKASICgood != 0 else row[16]

                    print('Row after change')
                    print(row)
                    break
                
            else:
            # Append new row
                rows.append([savepath, label, ST, baseline, targetmean, utargetmean, targetASIC, utargetASIC, 
                            targetASICgood, utargetASICgood, Kmean, Kstd, uKmean, KASIC, uKASIC,KASICgood, uKASICgood])
        

        with open(csv_file_path, mode='w') as file:
            writer = csv.writer(file)
            writer.writerows(rows)
        
        print_content = True
        if print_content:
            print('savepath, label, ST, baseline, targetmean, utargetmean, targetASIC, utargetASIC, targetASICgood, utargetASICgood, Kmean, Kstd, uKmean, KASIC, uKASIC')
            for row in rows:
                if row == existing_row:
                    print(row)





#path = '/data/bfys/apuicerc/N037_new/resultscheck/'
#prefix = 'E0'
#fo = FileOrganizer(path)
#label = 'warm'
#csv_files= FileOrganizer.get_csv_files(self = fo ,prefix = prefix)
##print('csv_files', csv_files, len(csv_files))
#prefix_files = FileOrganizer.open_csv(self = fo ,prefix = prefix)
##print('prefix_files', prefix_files, len(prefix_files))
##print('len(prefix_files)', )

#
##print('Files in {} starting with {}:'.format(path, prefix))
##csv_files = fo.get_csv_files()
###print('csv_files', csv_files)
#
##print('Is "Flux.csv" in the list?')
#print(fo.is_file_found(path,'target'))




