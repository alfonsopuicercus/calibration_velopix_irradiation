#import numpy as np
import glob as glob
import os
from typing import Callable
from utils import datainfo, organizer
from analyse import analyse
from plot import plot, plotcaller
import numpy as np
import inspect

class Run(plotcaller.PlotCaller,organizer.FileOrganizer):
    
    def __init__(
        self,
        task: dict = None
    ):
        """
        task list is built from client parser
        """
        self.task = task
        print('task', task)
        
        
        
    #analysis functions
    def analyse_flux(self, data_info):
            print('Analysing data and generating flux files from {}, saved to {} with temp {}'.format(
                data_info.datafolder_l, data_info.FPTfolder_l, data_info.label_l))
            analyse.AnalyseFlux(**data_info.__dict__)

    def analyse_fitter(self,baseline_info, data_info, s2_opt = None ):
             
        print('For fitting, p0 is necessary. For the data being processed:')
        try:
            p0 = data_info.p0
            print('p0 = {}'.format(p0))
        except: print('Error, p0 not set as property of dataset in datainfo.py')
        try:
            p0_baseline = data_info.p0_baseline
            print('p0_baseline = {}'.format(p0_baseline))
        except: print('Error, p0_baseline not set as property of dataset in datainfo.py')
        
        baseline_question = input('Do you want to consider baseline while fitting? (y/n): ')
        if baseline_question.lower() == 'y':
            print('Noise baseline is considered while fitting.')
            print('Fitting for {} {} {}'.format(data_info.datafolder_l,
            data_info.FPTfolder_l, data_info.label_l))
            analyse.Fitter(baseline_info.FPTfolder_l,s2_opt = s2_opt, **data_info.__dict__)
        else:
            print('Noise baseline not considered while fitting.')
            print('Fitting for {} {} {}'.format(data_info.datafolder_l,
            data_info.FPTfolder_l, data_info.label_l))
            analyse.Fitter(FPTfolder_baseline_l = None,s2_opt = s2_opt, **data_info.__dict__)

    def analyse_target(self,baseline_info, data_info ):
        print('Computing target matrices for {} {} {}'.format(data_info.datafolder_l,
            data_info.FPTfolder_l, data_info.label_l))
        matrix_env = analyse.Matrices(FPTfolder_baseline_l = baseline_info.FPTfolder_l, **data_info.__dict__)

        if baseline_info.FPTfolder_l != None:
            matrix_env.target_finder()
        else:
            print('Error FPTfolder_baseline_l not specified in datainfo.py for {}'.format(data_info.datafolder_l))

    
    def analyse_K(self,baseline_info, data_info,s2_opt = None ):
        print('Computing K[e-/DAC] matrices for {} {} {}'.format(data_info.datafolder_l,
            data_info.FPTfolder_l, data_info.label_l))
        matrix_env = analyse.Matrices(FPTfolder_baseline_l = baseline_info.FPTfolder_l,s2_opt = s2_opt, **data_info.__dict__)

        baseline_question = input('Do you want to consider baseline while fitting ASIC? (y/n): ')
        if baseline_question.lower() == 'y':
            print('Noise baseline is considered while fitting ASIC.')
            matrix_env.K_finder(FPTfolder_baseline = baseline_info.FPTfolder_l, baselineASICbeforefit=True, s2_opt = s2_opt)
        else:
            print('Noise baseline not considered while fitting ASIC.')
            matrix_env.K_finder(FPTfolder_baseline = baseline_info.FPTfolder_l, baselineASICbeforefit=None, s2_opt = s2_opt)

    # search functions
    def find(self, key: str, f: Callable, cond: Callable, id = '', option = '') -> None:
        print('Searching for ' + str(cond) + ' in ' + key)
        f(key, id, cond)

    # plotting functions
    def plot_heatmap(self, param,temp ,ST,opt = None, path = '') -> None:
        if param == 'totalflux':
            print('totalflux')
            self.plot_hmap(self.flux_heat_map, 'totalflux',temp ,ST, path = path)

        if param != None and param != 'Fittype' and param != 'totalflux':
            self.plot_hmap(self.heat_map, param,temp ,ST, path = path)
        if param == 'Fittype':
            self.plot_hmap(self.fittype_heat_map, 'Fittype',temp ,ST, path = path)

        


        if param == None:
            print('Plotting heatmaps for all parameters...')
            try: 
                self.plot_hmap(self.heat_map, 'target',temp ,ST, path = path)
            except: print('Error plotting heatmap of target.')

            try: 
                self.plot_hmap(self.heat_map, 'E0',temp ,ST, path = path)
            except: print('Error plotting heatmap of E0.')

            
            self.plot_hmap(self.heat_map, 'f',temp ,ST, path = path)
            self.plot_hmap(self.heat_map, 's',temp ,ST, path = path)
            self.plot_hmap(self.fittype_heat_map, 'Fittype',temp ,ST, path = path)
            #self.plot_hmap(self.heat_map, 'target',temp ,ST, path = path)

    def plot_results_collection(self, param,temp ,ST,opt = None, path = '') -> None:
        specified_result = path.split('/')[-2]

        categories = ['even', 'odd', 'rows']
        if specified_result in categories:
            specified_result = path.split('/')[-3] +'/'+path.split('/')[-2]
        else:
            specified_result = path.split('/')[-2]

        if param == None:
            if opt == 'summary':
                print("Choose an option to plot its summary:")
                print("1. Gain")
                print("2. Good Pixels")
                print("3. Both Gain and Pixels")
                print("4. Energy/electrons vs target")
                choice = int(input("Enter your choice (1, 2, 3 or 4): "))
      
                if choice == 1:
                    self.plot_summary()
                if choice == 2:
                    self.plot_summary(typeplot = 'pixels')
                if choice == 3:
                    self.plot_summary()
                    self.plot_summary(typeplot = 'pixels')
                if choice == 4:
                    self.plot_energy_vs_target()


           
            if opt == 'discrepancy':
                print('Plotting discrepancy histograms from collection for target and K ...')
                print("Which value should we compare to? Enter 'ASIC' or 'ASICgood':")
                which_to_compare = input()
                opt = which_to_compare
                self.plot_collection(specified_result = specified_result, opt = opt)

            if opt == None:
                print('Plotting histograms from collection for target and K ...')

                self.plot_collection(specified_result = specified_result, opt = opt, temp = temp)

            
    
    def distribution(self, param = None, temp = None, ST = None, opt = None, path = '') -> None:
        print('Plotting distributions...')
        print('temp', temp,'param', param)
  
        if opt == 'discrepancy':
            
            print("Which value should we compare to? Enter 'ASIC' or 'ASICgood':")
            which_to_compare = input()
            opt = which_to_compare
        self.plot_histograms(self.param_hist,param=param,temp = temp, ST = ST, opt = opt, path = path)

    def flux_over_thresholdscan(self,param = None, temp = None, ST = None, opt = None, path = '') -> None:
        print('Plotting distributions...')
        print('path', path)
        print('temp', temp,'param', param)
        print('opt ', opt)
        print('opt.label', opt.FPTfolder_l)
        temp_baseline_info = datainfo.DataInfo('baseline_information', temp)
        matrix_env = analyse.Matrices(FPTfolder_baseline_l = temp_baseline_info.FPTfolder_l, **opt.__dict__)
        # Prompt user for plot choice
        plot_question = input('Do you want to plot ASIC, ASICgood, both, or neither? (asic/asicgood/y/n): ').lower()

        # Initialize variables
        ASIC_information = None
        ASICgood_information = None

        # Determine whether to plot and set up baseline question if needed
        if plot_question in {'asic', 'asicgood', 'y'}:
            # Ask about baseline only if plotting is chosen
            baseline_question = input('Do you want to consider baseline while fitting ASIC? (y/n): ').lower()
            baselineASICbeforefit = True if baseline_question == 'y' else None
            print('baselineASICbeforefit plot', baselineASICbeforefit)
            print(f'Noise baseline {"is" if baselineASICbeforefit else "is not"} considered while fitting ASIC.')

            # Perform the appropriate plotting
            if plot_question == 'asic':
                print('You have chosen to plot ASIC.')
                ASIC_information = matrix_env.find_ASIC(FPTfolder_baseline=temp_baseline_info.FPTfolder_l, baselineASICbeforefit=baselineASICbeforefit, plot_output=True)
            elif plot_question == 'asicgood':
                print('You have chosen to plot ASICgood.')
                ASICgood_information = matrix_env.find_ASICgood(FPTfolder_baseline=temp_baseline_info.FPTfolder_l, baselineASICbeforefit=baselineASICbeforefit, plot_output=True)
            elif plot_question == 'y':
                print('You have chosen to plot both ASIC and ASICgood.')
                ASIC_information = matrix_env.find_ASIC(FPTfolder_baseline=temp_baseline_info.FPTfolder_l, baselineASICbeforefit=baselineASICbeforefit, plot_output=True)
                ASICgood_information = matrix_env.find_ASICgood(FPTfolder_baseline=temp_baseline_info.FPTfolder_l, baselineASICbeforefit=baselineASICbeforefit, plot_output=True)

        elif plot_question == 'n':
            print('You have chosen not to plot either ASIC or ASICgood.')
            

        choice_question = input('Do you want to plot pixel, pixels in area, default pixel or no pixels? (pixel/area/y/n): ').lower()

        if choice_question == 'pixel':
            pixel_x = int(input('Pixel x coordinate? (0-255)'))
            pixel_y = int(input('Pixel y coordinate? (0-255)'))
            radius = 1
        elif choice_question == 'area':
            pixel_x = int(input('Pixel x coordinate? (0-255)'))
            pixel_y = int(input('Pixel y coordinate? (0-255)'))
            radius = int(input('Radius of area to plot (in number pixels) (it plots 5 random pixels in this area)'))
        elif choice_question == 'y':
            pixel_x = 127
            pixel_y = 128
            radius = 0
            print('Plotting pixel {}x{} in center of exposure...'.format(pixel_x,pixel_y))
        else:
            pixel_x = None
            pixel_y = None
            radius = 0

        
        self.ind_pixel_plotter(pixel_row = pixel_x, pixel_column = pixel_y, radii = radius, baseline_opt = False,FPTfolder_l = path, ST_l = ST,label_l = temp,p0 = opt.p0, ASIC_information = ASIC_information, ASICgood_information = ASICgood_information, together ='False')
        
## MD thanks
    # monitoring functions
    def plot_mask(self, __, path = '') -> None:
        self.plot_cv('m_map', self.mask_map, path = path)

    def plot_trim(self, __, path = '') -> None: 
        self.plot_cv('t_hist', self.trim_hist, path = path)
        self.plot_cv('t_map', self.trim_map, path = path)

    def plot_normal(self, __, path = '') -> None:
        self.plot_cv('n_summary', self.n_summary, path = path)

    def plot_control(self, __, path = '') -> None:
        self.plot_cv('eq_summary', self.eq_summary, path = path)

    def plot_nc(self, __, path = '') -> None:
        self.plot_cv('nm_corr', self.n_corr, opt = 'mean', path = path)
        self.plot_cv('nr_corr', self.n_corr, opt = 'rate', path = path)

    def plot_ns(self, scan, path = '') -> None:
        if scan == 'normal' or scan == 'quick':
            self.plot_cv('ns_t0_map', self.nw_map, opt = '0', path = path)
            self.plot_cv('ns_t15_map', self.nw_map, opt = '15', path = path)
            self.plot_cv('ns_t0_hist', self.nw_hist, opt = '0', path = path)
            self.plot_cv('ns_t15_hist', self.nw_hist, opt = '15', path = path)
        elif scan == 'control':
            self.plot_cv('ns_t16_map', self.nw_map, opt = '16', path = path)
            self.plot_cv('ns_t16_hist', self.nw_hist, opt = '16', path = path)

    def plot_nr(self, scan, path = '') -> None:
        if scan == 'normal' or scan == 'quick':
            self.plot_cv('nr_t0_map', self.n_map, opt = 'n_rate_t0', path = path)
            self.plot_cv('nr_t15_map', self.n_map, opt = 'n_rate_t15', path = path)
        elif scan == 'control':
            self.plot_cv('nr_t16_map', self.n_map, opt = 'n_rate_t16', path = path)

    def plot_nm(self, scan, path = '') -> None:
        if scan == 'normal' or scan == 'quick':
            self.plot_cv('nm_t0_map', self.n_map, opt = 'n_mean_t0', path = path)
            self.plot_cv('nm_t15_map', self.n_map, opt = 'n_mean_t15', path = path)
        elif scan == 'control':
            self.plot_cv('nm_t16_map', self.n_map, opt = 'n_mean_t16', path = path)

 
    def plot_pix(self) -> None:
        pix = self.task['pix'][0] + '_' + self.task['pix'][1] + '_' + self.task['pix'][2]
        self.plot('pixel_' + pix, self.pix_scan, id = pix, path = path)
    
    def plot_all(self, scan, path = '') -> None:
        if scan == 'normal' or scan == 'quick':
            self.plot_trim(scan, path = path)
            self.plot_normal(scan, path = path)
            self.plot_nc(scan, path = path)
            self.plot_ns(scan, path = path)
            self.plot_nr(scan, path = path)
            self.plot_nm(scan, path = path)
            self.plot_dist(scan, path = path)
        if scan == 'control':
            self.plot_mask(scan, path = path)
            self.plot_control(scan, path = path)
            self.plot_ns(scan, path = path)
            self.plot_nr(scan, path = path)
    

    # main manager
    def parse(self) -> None:
        # Parse the command
        if not any(self.task.values()) is True: #if  no task is specified, ask for input
            print("Empty task.")
            return
            
        if self.task["which"] == "analysis":
            # Extract information about the dataset and temperature
            data_info = datainfo.DataInfo(self.task['dataset'], self.task['temp'])


            baseline_info = datainfo.DataInfo('baseline_information', self.task['temp'])
            
            # Create a FileOrganizer for the dataset
            fo = organizer.FileOrganizer(data_info.FPTfolder_l[0])
            
            # If no temperature is specified, iterate over all temperature labels
            if self.task['temp'] is None:
                # Extract all temperature labels in case more than one temperature to be analysed
                temp_labels = data_info.label_l
                
                # Loop over each temperature label (e.g. warm and cold in this case)
                for temp_label in temp_labels:
                    # Extract information about the dataset for the current temperature
                    temp_info = datainfo.DataInfo(self.task['dataset'], temp_label, separation =self.task['separation'])
                    temp_baseline_info = datainfo.DataInfo('baseline_information', temp_label)
                    # Create a FileOrganizer for the current dataset

                    temp_fo = organizer.FileOrganizer(path = temp_info.FPTfolder_l)
                    
                    # Depending on the generator, perform analyses on the current dataset
                    if self.task['generator'] == 'flux': 
                        # If the fluxperASIC files do not exist, analyse the dataset
                        if not temp_fo.is_file_found(prefix='FluxperASIC'):
                            self.analyse_flux(temp_info)
                        else: 
                            print('Flux files already exist for {}.'.format(temp_label))
                            again = input('Do you still want to find them? (y/n): ')
                            if again.lower() == 'y':
                                self.analyse_flux(temp_info)

                    elif self.task['generator'] == 'fitter':
                        # Create the directory structure for the current dataset

                        temp_fitfunction = temp_info.fit_function
                        temp_signature = inspect.signature(temp_fitfunction)
                        temp_params_l = list(temp_signature.parameters)[1:]
   
                        temp_fo.create_directory_structure(temp_info.FPTfolder_l, list_of_params = temp_params_l)
                        
                        
                        # If the E0 files do not exist, analyse the dataset
                        if not temp_fo.is_file_found(prefix='E0'):
                            self.analyse_fitter(temp_baseline_info, temp_info, s2_opt = self.task['s2_opt'])
                        
                        else: 
                            print('Fitter matrices already exist for {}.'.format(temp_label))
                            again = input('Do you still want to find them? (y/n): ')
                            if again.lower() == 'y':
                                self.analyse_fitter(temp_baseline_info, temp_info, s2_opt = self.task['s2_opt'])

                    elif self.task['generator'] == 'matrices':
                        # If the target and K matrices do not exist, analyse the dataset
                    
                        if not temp_fo.is_file_found(prefix='target', label = temp_label):
                            self.analyse_target(temp_baseline_info, temp_info)
                        else: 
                            print('Target matrices already exist for {}.'.format(temp_label))
                            again = input('Do you still want to find them? (y/n): ')
                            if again.lower() == 'y':
                                self.analyse_target(temp_baseline_info, temp_info)

                        if not temp_fo.is_file_found(prefix='K', label = temp_label):
                            self.analyse_K(temp_baseline_info, temp_info, s2_opt = self.task['s2_opt'])
                        else: 
                            print('K matrices already exist for {}.'.format(temp_label))
                            again = input('Do you still want to find them? (y/n): ')
                            if again.lower() == 'y':
                                self.analyse_K(temp_baseline_info, temp_info, s2_opt = self.task['s2_opt'])

            else:
                data_info = datainfo.DataInfo(self.task['dataset'], label_to_consider = self.task['temp'], separation =self.task['separation'])
                fo = organizer.FileOrganizer(data_info.FPTfolder_l)

                # Obtaining fit function and its parameters
                fitfunction = data_info.fit_function
                signature = inspect.signature(fitfunction)
                params_l = list(signature.parameters)[1:] # list of parameter names in a list

                if self.task['generator'] == 'flux': 
                    if not fo.is_file_found(prefix='FluxperASIC'):
                        self.analyse_flux(data_info)
                    else: 
                        print('Flux files already exist for {}.'.format(data_info.label_l))
                        again = input('Do you still want to find them? (y/n): ')
                        if again.lower() == 'y':
                            self.analyse_flux(data_info)

                elif self.task['generator'] == 'fitter':
                    
                    #creating directory structure for plots and subfolders for flux scans and heatmap plots
                    fo.create_directory_structure(data_info.FPTfolder_l, list_of_params = params_l)

                    
                    if not fo.is_file_found(prefix='E0'):
                        self.analyse_fitter(baseline_info, data_info, s2_opt = self.task['s2_opt'])
                    else: 
                        print('Fitter matrices already exist')
                        again = input('Do you still want to find them? (y/n): ')
                        if again.lower() == 'y':
                            self.analyse_fitter(baseline_info, data_info, s2_opt = self.task['s2_opt'])

                elif self.task['generator'] == 'matrices':
                    if not all(fo.is_file_found(prefix=pfx, label=data_info.label_l)
                               for pfx in ('target', 'K')):
                        self.analyse_target(baseline_info, data_info)
                        self.analyse_K(baseline_info, data_info, s2_opt = self.task['s2_opt'])
                    else: 
                        print('Target and K matrices already exist')
                        again = input('Do you still want to find them? (y/n): ')
                        if again.lower() == 'y':
                            self.analyse_target(baseline_info, data_info, s2_opt = self.task['s2_opt'])
                            self.analyse_K(baseline_info, data_info, s2_opt = self.task['s2_opt'])

        if self.task["which"] == "plot":
            print('new plot code')
            data_info = datainfo.DataInfo(self.task['dataset'], self.task['temp'], self.task['separation'])
            path = data_info.FPTfolder_l #+ 'plots/'
            ST = data_info.ST_l
            label = None
            if self.task['temp'] is None and self.task['plot_results_collection'] is not None:
                for element in self.task:
                    if self.task[element] != None and element != 'which' and element != 'dataset' and element != 'param' and element != 'temp' and element != 'opt' and element != 'separation':
                        getattr(self, element)(self.task['param'],label, ST, opt = self.task['opt'],path = path)
            else:
                temp_labels = data_info.label_l
                for temp_label in temp_labels:
                    temp_info = datainfo.DataInfo(self.task['dataset'], temp_label, self.task['separation'])
                    temp_fo = organizer.FileOrganizer(path = temp_info.FPTfolder_l)
                    temp_fo.create_directory_structure(temp_info.FPTfolder_l)
                    path = temp_info.FPTfolder_l #+ 'plots/'
                    ST = temp_info.ST_l
                    temp_baseline_info = datainfo.DataInfo('baseline_information', temp_label)
                    for element in self.task:
                        if self.task[element] != None and element != 'which' and element != 'dataset' and element != 'param' and element != 'temp' and element != 'opt' and element != 'separation':
                            if self.task['flux_over_thresholdscan'] is not None: #if plot flux over thrs, provide all info of dataset to the plotting function
                                opt = temp_info # using opt to pass the dataset info at this temperature to the plotting function
                            else:
                                opt = self.task['opt']
                            getattr(self, element)(self.task['param'],temp_label, ST, opt = opt, path = path)

        
        
        
        
        #plot manager
        if self.task["which"] == "plotold":
            data_info = datainfo.DataInfo(self.task['dataset'], self.task['temp'], self.task['separation'])
            fo = organizer.FileOrganizer(data_info.FPTfolder_l)

            # Obtaining fit function and its parameters
            fitfunction = data_info.fit_function
            signature = inspect.signature(fitfunction)
            params_l = list(signature.parameters)[1:] # list of parameter names in a list


            print('path ', data_info.FPTfolder_l )
            if self.task['temp'] is None and self.task['plot_results_collection'] is not None:
                print('-rc is not None path ', data_info.FPTfolder_l )
                path = data_info.FPTfolder_l[0] #+ 'plots/'
                ST = data_info.ST_l[0]
                label = None
                print('path = {}, ST = {}, label = {}'.format(path, ST, label))

                for element in self.task:
                    #print(element)
                    if self.task[element] != None and element != 'which' and element != 'dataset' and element != 'param' and element != 'temp' and element != 'opt' and element != 'separation':
                        #print('true',element, self.task[element])
                        getattr(self, element)(self.task['param'],label, ST, opt = self.task['opt'],path = path)



            elif self.task['temp'] is None and self.task['plot_results_collection'] is None:
                temp_labels = data_info.label_l
                for temp_label in temp_labels:

                    temp_info = datainfo.DataInfo(self.task['dataset'], temp_label, self.task['separation'])
                    temp_fo = organizer.FileOrganizer(path = temp_info.FPTfolder_l)
                    temp_fo.create_directory_structure(temp_info.FPTfolder_l)
                    
                    path = temp_info.FPTfolder_l #+ 'plots/'
                    ST = temp_info.ST_l

                    temp_baseline_info = datainfo.DataInfo('baseline_information', temp_label)

                    #print('its going through here, temp_label = {} {}'.format(temp_label,  self.task['opt']))
                    if self.task['all'] is None: 
                        for element in self.task:
                        
                            if self.task[element] != None and element != 'which' and element != 'dataset' and element != 'param' and element != 'temp' and element != 'opt' and element != 'separation':
                                print('true',element, self.task[element])
                                if self.task['flux_over_thresholdscan'] is not None: #if plot flux over thrs, provide all info of dataset to the plotting function
                                    opt = temp_info # using opt to pass the dataset info at this temperature to the plotting function
                                else:
                                    opt = self.task['opt']

                                getattr(self, element)(self.task['param'],temp_label, ST, opt = opt, path = path)
            else:
                data_info = datainfo.DataInfo(self.task['dataset'], self.task['temp'], self.task['separation'])
                if self.task['all'] is None:
                    print('all is none')
                    for element in self.task:
                        #print(element)
                       
                        if self.task[element] != None and element != 'which' and element != 'dataset' and element != 'param' and element != 'temp' and element != 'separation':
                            #print('true',element, self.task[element])
                            path = data_info.FPTfolder_l #+ 'plots/'
                            ST = data_info.ST_l
                           
                            getattr(self, element)(self.task['param'],self.task['temp'], ST,  opt = self.task['opt'], path = path)
                           


        
        if self.task['which'] == 'noise':
            if self.task['analysis']:
                #self.noise_study(self.task['analysis'])
                #self.spread_study(self.task['analysis'])
                #self.map_baseline(self.task['analysis'])
                self.nsplit_study(self.task['analysis'])

        if self.task['which'] == 'search':
            key = self.task['key']
            cond = lambda a: eval(self.task['criterion'])
            id = self.task['velopix']
            self.find(key, self.find_pixel, cond = cond, id = id)

        # plot manager
        if self.task['which'] == 'idkplot':
            print('task plot')
            data_info = datainfo.DataInfo(self.task['dataset'], self.task['temp'])
            if self.task['velo'] is None: 
                if self.task['all'] is None: 
                    print(self.task)
                    for element in self.task:
                        #print(element)
                       
                        if self.task[element] != None and element != 'which' and element != 'dataset' and element != 'param' and element != 'temp':
                            print('if true',element, self.task[element])
                            path = data_info.FPTfolder_l #+ 'plots/'
                            ST = data_info.ST_l
                            print('path', path, 'ST', ST)
                            getattr(self, element)(self.task['param'],self.task['temp'], ST, path = path)
                           
                elif self.task['all'] is True:
                    print('all = True')
                    self.plot_all(self.task['plot_type'])

            elif self.task['velo'] is True:
                for mod in range(52):
                    path = 'Module' + str(mod).zfill(2) + '/' + self.task['scan_type'] + '/'
                    if os.path.exists(path):
                        if not os.path.exists(path + 'plot') and self.task["which"] == "plot":
                            os.makedirs(path + 'plot')
                        #for id in asicinfo.asic_list:
                        if self.task['all'] is None:
                            for element in self.task:
                                if self.task[element] != None and element != 'which' and element != 'scan_type' and element != 'velo':
                                    getattr(self, element)(self.task['scan_type'], path = path)
                        elif self.task['all'] is True:
                            self.plot_all(self.task['scan_type'], path = path)

        if self.task["which"] == "monitoring":
            if self.task['history']:
                self.group_history(25)

        if self.task['which'] == 'status':
            if self.task['status_of'] == 'equalisation':
                self.status_manager(self.task['dir'])
            if self.task['status_of'] == 'update':
                #self.make_overview(self.task['dir'], which = 'quick')
                #self.make_overview(self.task['dir'], which = 'normal')
                self.make_overview(self.task['dir'], which = 'control')
            if self.task['status_of'] == 'masked': 
                self.status_masked(self.task['dir'])
            if self.task['status_of'] == 'spread':
                self.status_spread(self.task['dir'])
            if self.task['status_of'] == 'threshold':
                self.status_threshold(self.task['dir'])
            if self.task['status_of'] == 'baseline':
                self.status_mean(self.task['dir'])
                #self.status_mean(self.task['dir'])
                #self.status_if_actual(self.task['dir'])
            if self.task['status_of'] == 'dacsettings':
                #self.status_dacsettings(self.task['dir'], self.task['tag']) 
                self.correlation_dacsettings_noise(self.task['dir'], self.task['tag'])
            if self.task['status_of'] == 'summary':
                self.make_calib_summary(self.task['dir'])

        if self.task['which'] == 'recipe':
            if self.task['recipe_of'] == 'threshold':
                rc = recipehandler.RecipeHandler()
                rc.build_thr_recipe(self.task['dir'])

        #if self.task['which'] == 'ivelo':
            #gi = ivelo.Ivelo()
            #gi.run(self.task['dir'])


