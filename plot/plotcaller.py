import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from plot import plot
from analyse import analyse
import pandas as pd
import numpy as np

class PlotCaller(plot.Plot, analyse.Matrices):
    def __init__(
        self
    ):
        """
        template call functions to plot functions
        canvas definition
        """

    def plot(self, fname, f, id = '', opt = None, path = '') -> None:
        output_path = path + 'plot/' + fname + '.png'
        print('Creating ' + output_path + ' ...')
        fig, ax = plt.subplots(figsize = [6,5])
        f(ax, id, opt = opt, path = path) #this is a function in plot.py which takes data using open_csv plots in the figure created here.
        fig.savefig(output_path)
        plt.close('all')

    def plot_histograms(self, f, param = None, temp = None, ST = None, opt = None, s2_opt = None, path = ''):
        if param is None:
            fname = 'ParameterHistogram_'+path.split('/')[-2]+'ST'+ST+'_'+str(temp)
        else:
            fname = param +'Histogram_'+path.split('/')[-2]+'ST'+ST+'_'+str(temp)
        
        if opt is None:
            output_path = path + 'plots/' + fname + '.png'
        else:
            self.check_and_create_path(path + 'plots/discrepancy/')
            output_path = path + 'plots/discrepancy/' + opt + 'Discrepancy_' + fname + '.png'
        print('Creating ' + output_path + ' ...')
        
        if param != None:
            params_used = []
            params_used.append(param)
        else:
            if s2_opt != None:
                params_used = ['E0', 'f', 's','s2']
            else:
                params_used = ['E0', 'f', 's']
                #params_used = ['E0', 'f', 's', 'target']
        

        fig, axs = plt.subplots(1, len(params_used), figsize=(5*len(params_used), 5), dpi=300)
        if len(params_used) == 1:
            f(axs, params_used[0], temp=temp, opt = opt, path=path)
        else:
            for i, parameter in enumerate(params_used):
             
                f(axs[i], parameter, temp=temp,opt = opt, path=path)
        
        fig.tight_layout()
        

        fig.savefig(output_path)
        plt.close('all')


    
    def plot_hmap(self, f, param,temp, ST, path = '') -> None:

        fname = param +'_Heatmap_'+ path.split('/')[-2]+'_ST_'+ST+'_'+temp
        output_path = path + 'plots/'+'heat_map_plots/'+param +'/'+ fname + '.png'

        print('Creating ' + output_path + ' ...')
        plt.figure(figsize = [7,5], dpi = 300)
        f(param,temp, path = path) #this is a function in plot.py which takes data using open_csv plots in the figure created here.
        plt.savefig(output_path)
        plt.close('all')


    def plot_summary(self, results_selected=['calibration_baseline'], results_collection_path='/data/bfys/apuicerc/N037_new/results_collection.csv', typeplot = 'gain',savepath='/data/bfys/apuicerc/N037_new/plots'):
        
        datasets_to_plot = 'Main'

        #'resultscalibration_separation/even', 'resultscalibration_separation/odd', 'resultscalibration_separation/rows'
        #'resultscalibration_separation/even', 'resultscalibration_separation/odd','resultscalibration_separationspp/even', 'resultscalibration_separationspp/odd',
        # 'calibration_baseline','resultsultimate_newbaseline', 'resultscalibrationfitfunc4withterm0'
        
        #yaxis_results_selected = ['calibration', 'even', 'odd', 'rows']
        #yaxis_results_selected = ['calibration', 'even', 'odd', 'evenspp', 'oddspp']

        ## FOR BASELINE INCLUSION VARIATION
        #['resultscalibration', 'calibration_baseline','resultsultimate_newbaseline', 'resultscalibration_newbaseline']
        #yaxis_results_selected = ['calibration', '_baseline', '_newbaseline'] #'fitfunc4withterm0'
        if datasets_to_plot == 'Main':
            results_selected=['calibration_baseline']
            yaxis_results_selected = ['calibration']
            plotting_opt = [True, True, True] #corresponding to value to be plotted, mean, ASIC, ASICgood


        ## FOR EXPOSURE TIME VARIATION resultscalibration
        if datasets_to_plot == 'Exposuretime':
            results_selected = ['calibration_baseline','resultscalibrationequalT', 'resultscalibration_nacq25','resultscalibration_nacq25to50']
            yaxis_results_selected = ['calibration', 'calibration_nacq50', 'calibration_acq0to25','calibration_acq25to50'] #'fitfunc4withterm0'
            plotting_opt = [False, False, False]
        ## FOR FITFUNCTION VARIATION
        if datasets_to_plot == 'Fitfunction':
        #'resultscalibration', 'resultscalibrationnoterm0', 'resultscalibrationfitfunc2', 'resultscalibrationfunc2noterm0', 'resultscalibrationfitfunc3', 'resultscalibrationfitfunc3withterm0','resultscalibrationfunc2sAB'
            results_selected = ['resultscalibration', 'resultscalibrationnoterm0','resultscalibrationfitfunc3withterm0','resultscalibrationfitfunc3', 'resultscalibrationfitfunc2', 'resultscalibrationfunc2noterm0','resultscalibrationfunc2sAB']
        #yaxis_results_selected = ['calibration', 'resultscalibrationnoterm0','resultscalibrationfitfunc3withterm0','resultscalibrationfitfunc3', 'resultscalibrationfitfunc2', 'resultscalibrationfunc2noterm0','resultscalibrationfunc2sAB'] #'fitfunc4withterm0'
            yaxis_results_selected = ['calibration, $F_{1,0}$', '$F_{1}$','$F_{AB,0}$','$F_{AB}$', '$F_{2s,0}$', '$F_{2s}$','$F_{2s,AB}$'] #'fitfunc4withterm0'
        ##  CHANGE THE VARIABLE BELOW TO TRUE IF YOU WANT TO REMOVE A SPECIFIC KEYWORD FROM THE LABELS
            plotting_opt = [False, False, False]

        remove_word_from_labels = False
        if remove_word_from_labels:
            # Define a function to handle the special cases and general case
            def remove_results_prefix(s, index, keyword_to_remove):
                if index == 0:
                    return 'calibration'
                elif index == 1:
                    return 'calibrationnoterm0'
                elif s.startswith(keyword_to_remove):
                    return s[len(keyword_to_remove):]
                return s

            # Remove the prefix 'results' from each string in the list
            keyword_to_remove = 'resultscalibration'

            # Use list comprehension with enumeration to apply the function to each element
            cleaned_results_for_label = [remove_results_prefix(s, i, keyword_to_remove) for i, s in enumerate(results_selected)]

        
            yaxis_results_selected = cleaned_results_for_label
        else:
            yaxis_results_selected = yaxis_results_selected
        
        data_storage = {
            'Kmean_list': {
                'warm': [],
                'cold': []
            },
            'uKmean_list': {
                'warm': [],
                'cold': []
            }
            ,
            'KASIC_list': {
                'warm': [],
                'cold': []
            }
            ,
            'uKASIC_list': {
                'warm': [],
                'cold': []
            },
            'KASICgood_list': {
                'warm': [],
                'cold': []
            }
            ,
            'uKASICgood_list': {
                'warm': [],
                'cold': []
            }
        }
        
        label_l = ['warm', 'cold']
        
        if typeplot == 'gain':
            print("Reading data from " + results_collection_path)
            data = pd.read_csv(results_collection_path)

                    
            mean_opt, ASIC_opt, ASICgood_opt = plotting_opt

            
            
            #Set the below variable to True if you want to have a cut in the x axis
            #Convinient if there is a large difference in value for the gain
            cut_x_axis_opt=False

            if cut_x_axis_opt:
                fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(15, 5), dpi=300)
                # Adjust the spacing between the subplots
                fig.subplots_adjust(wspace=0.025)

                x1_lim_min, x1_lim_max = 13, 13.5#13.3
                x2_lim_min, x2_lim_max = 14.1, 14.5
                step_ticks = 0.1

                ax1.set_xlim( x1_lim_min, x1_lim_max)
                ax2.set_xlim(x2_lim_min, x2_lim_max)
                ax1.set_xticks(np.arange(x1_lim_min, x1_lim_max, step_ticks))
                ax2.set_xticks(np.arange(x2_lim_min, x2_lim_max+step_ticks, step_ticks))
                axes = [ax1, ax2]
                
                    # Hide the right spine on the first plot and the left spine on the second plot
                ax1.spines['right'].set_visible(False)
                ax2.spines['left'].set_visible(False)
                
                # Add diagonal lines to indicate the break
                d = 0.015  # size of diagonal lines in axes coordinates
                kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
                ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)        # top-left diagonal
                ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

                kwargs.update(transform=ax2.transAxes)  # switch to the right axes
                ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-right diagonal
                ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal


            else:
                fig, ax1 = plt.subplots(1, 1, sharey=True, figsize=(15, 5), dpi=300)
                
                x_lim_min, x_lim_max = 12.25, 15
                step_ticks = 0.25

                ax1.set_xlim(x_lim_min, x_lim_max)
                ax1.set_xticks(np.arange(x_lim_min, x_lim_max+step_ticks, step_ticks))

                axes = [ax1]
            
            categories = ['even', 'odd', 'rows']
            save_fig_label = ''
            for i in range(len(results_selected)):
                
                specified_result = results_selected[i]
                
                if specified_result.split('_')[-1].split('/')[0] == 'separation':
                    save_fig_label = '_separation'

                    for cat in categories:
                        separation_specified_result = specified_result + '/' + cat
            
            for j in range(len(results_selected)):
                
                specified_result = results_selected[j]
                
                print("Plotting " + specified_result)
                result_params = data[data['Results'] == specified_result].drop('Results', axis=1)
                
                for i in range(len(result_params)):

                    # Get the parameters for the current result
                    label = result_params.iloc[i].values[0]
                    ST = result_params.iloc[i].values[1]
                    baseline = result_params.iloc[i].values[2]
                    result_param = [float(x) for x in result_params.iloc[i].tolist()[3:]]
                    
                    targetmean, utargetmean, targetASIC, utargetASIC, targetASICgood, utargetASICgood, Kmean, Kstd, uKmean, KASIC, uKASIC, KASICgood, uKASICgood = result_param
                    # now result_param contains all the floats from results_collection.csv
                    print('Plotting for dataset {} {}'.format(specified_result, label))

                    color = 'blue' if label.startswith('cold') else 'orange'

                    custom_separation_ASIC = 0.1
                    custom_separation_ASICgood = 0.2
                    elinewidth = 2
                    capsize = 5
                    
                    if specified_result == 'resultscalibration' or specified_result == 'calibration_baseline':
                        for ax in axes:
                            ax.vlines(Kmean, -1, len(results_selected), colors='grey', linestyles='--')
                            ax.fill_betweenx([-1, len(results_selected)], Kmean - uKmean, Kmean + uKmean, color=color, alpha=0.2)

                    alpha = 1.0
                    
                    if label.startswith('cold'):
                            data_storage['Kmean_list']['cold'].append(Kmean)
                            data_storage['uKmean_list']['cold'].append(uKmean)
                            data_storage['KASIC_list']['cold'].append(KASIC)
                            data_storage['uKASIC_list']['cold'].append(uKASIC)
                            data_storage['KASICgood_list']['cold'].append(KASICgood)
                            data_storage['uKASICgood_list']['cold'].append(uKASICgood)
                            alpha = 1
                    else:
                            data_storage['Kmean_list']['warm'].append(Kmean)
                            data_storage['uKmean_list']['warm'].append(uKmean)
                            data_storage['KASIC_list']['warm'].append(KASIC)
                            data_storage['uKASIC_list']['warm'].append(uKASIC)
                            data_storage['KASICgood_list']['warm'].append(KASICgood)
                            data_storage['uKASICgood_list']['warm'].append(uKASICgood)
                    

                    
                    for ax in axes:
                        if mean_opt :
                            label_mean = label + ', Mean = {:.2f}$\pm${:.2f}'.format(Kmean, uKmean)
                        elif j == 0:
                            label_mean =  label
                        else:
                            label_mean = ''
                        ax.errorbar(Kmean, j, xerr=uKmean, fmt='o', ecolor=color,
                                    capsize=capsize, elinewidth=elinewidth, color=color, alpha=alpha,
                                    label=label_mean )
                        
                        
                        if ASIC_opt:
                            ax.errorbar(KASIC, j - custom_separation_ASIC, xerr=uKASIC, fmt='o', ecolor='red',
                                        capsize=capsize, elinewidth=elinewidth, color='red',
                                        label='ASIC' if j == 0 and i == 1 else '')
                        
                        if ASICgood_opt:
                            ax.errorbar(KASICgood, j - custom_separation_ASICgood, xerr=uKASICgood, fmt='o', ecolor='purple',
                                        capsize=capsize, elinewidth=elinewidth, color='purple',
                                        label='ASICgood' if j == 0 and i == 1 else '')
            
            if len(results_selected) != 1:

                sys_uncertainty_l = []
                for label in label_l:
                    sys_uncertainty = self.systematic_uncertainty_finder(results_selected,data_storage['Kmean_list'][label], data_storage['uKmean_list'][label])
                    sys_uncertainty_l.append(sys_uncertainty)

                print('Systematic uncertainties found for warm and cold', sys_uncertainty_l)
            else:
                sys_uncertainty_l = self.systematic_uncertainty_finder(results_selected,data_storage['Kmean_list'][label], data_storage['uKmean_list'][label])


            for ax in axes:
                ax.set_xlabel('K [e-/DAC]')
                ax.set_ylim(-0.5, len(results_selected) - 0.5)
                ax.set_yticks(np.arange(len(results_selected)))
                ax.set_yticklabels(yaxis_results_selected)
                ax.grid()



            
            

            extra_plotting = False
            systematic_uncertainty_opt = True
            estimate_gain_opt = True

            custom_handles = []
            if extra_plotting == True:
                if systematic_uncertainty_opt == True:
                    dummycounter = 0

                    for temp in label_l:
                        
                        Kmean_reference = data_storage['Kmean_list'][temp][0]
                        uKmean_reference = data_storage['uKmean_list'][temp][0]
                        sys_uncertainty_label = sys_uncertainty_l[dummycounter]

                        meanvaluelegend =  ' K[e-DAC] = {:.2f}$\pm${:.2f}$\pm${:.2f}'.format(Kmean_reference, uKmean_reference,sys_uncertainty_label)
                        color = 'blue' if temp.startswith('cold') else 'orange'
                        custom_handles.append(plt.Line2D([0], [0], color=color, lw=elinewidth, label=temp + meanvaluelegend))

                        dummycounter += 1
                        for ax in axes:
                            ax.fill_betweenx([-1, len(results_selected)], Kmean_reference - uKmean_reference - sys_uncertainty_label, Kmean_reference + uKmean_reference + sys_uncertainty_label, color='grey', alpha=0.2, label = 'Systematic uncertainty')
                
                if estimate_gain_opt == True:
                    # Given values
                    DAC_mV = 0.38  # mV per DAC
                    G_paper_mV_per_ke = 24.6  # mV per ke⁻
                    gain_variation_percentage = 3.3/100
                    # Convert the gain from mV/ke⁻ to mV/e⁻
                    G_paper_mV_per_e = G_paper_mV_per_ke / 1000  # because 1 ke⁻ = 1000 e⁻
                    
                    # Calculate the number of electrons per DAC
                    K_paper = DAC_mV / G_paper_mV_per_e
                    estimate_std = K_paper * gain_variation_percentage
                    tag = 'Estimate paper'

                    ax.errorbar(K_paper, j, xerr=estimate_std, fmt='o', ecolor='green',
                                    capsize=capsize, elinewidth=elinewidth, color='green', alpha=alpha,
                                    label='Estimate K[e-/DAC] = ' + '{:.2f}$\pm${:.2f}'.format(K_paper, estimate_std) )

                    custom_handles.append(plt.Line2D([0], [0], color='green', alpha=alpha,
                                    label='Estimate K[e-/DAC] = ' + '{:.2f}$\pm${:.2f}'.format(K_paper, estimate_std)))

                    # Calculate the difference between K_paper and Kmean
                    Kmean_reference = data_storage['Kmean_list']['cold'][0]
                    uKmean_reference = data_storage['uKmean_list']['cold'][0]
                    sys_uncertainty_reference = sys_uncertainty_l[1]

                    def sigma_difference_comparer(Kmean_reference, uKmean_reference, sys_uncertainty_reference, K_compare, uK_compare,tag = None):
                        print('Comparing {} dataset to {} ...'.format(label, tag))
                        difference = K_compare - Kmean_reference
                        
                        # Here calculating the total uncertainty associated to the data point: including the error of the mean and systematic uncertainty using quadrature
                        total_error_K_data = np.sqrt(sys_uncertainty_reference**2 + uKmean_reference**2)
                        #total_error_K_data = uKmean_reference
                        # Combine uncertainties
                        combined_uncertainty = np.sqrt(uK_compare**2 + total_error_K_data**2)

                        # Calculate number of standard deviations
                        sigma_difference = difference / combined_uncertainty

                        print('Distance in standard deviations from {}: {:.2f}'.format(tag, sigma_difference))


                    sigma_difference_comparer(Kmean_reference, uKmean_reference, sys_uncertainty_reference, K_paper, estimate_std, tag)
                    
                    temp_reference = 'cold'
                    Kmean_reference = data_storage['Kmean_list'][temp_reference][0]
                    uKmean_reference = data_storage['uKmean_list'][temp_reference][0]
                    values_to_compare = ['KASIC','KASICgood']
                    for compare in values_to_compare:
                        K_compare = data_storage[compare+'_list'][temp_reference][0]
                        uK_compare = data_storage['u'+compare+'_list'][temp_reference][0]
                    
                        sigma_difference_comparer(Kmean_reference, uKmean_reference, sys_uncertainty_reference, K_compare, uK_compare, tag=compare)
            

            meanvalueslegend = False
            if meanvalueslegend:
                meanvaluelegend =  '{:.2f}$\pm${:.2f}'.format(Kmean, uKmean)
            else:
                meanvaluelegend = ''

            if extra_plotting == False:
                custom_handles = [
                plt.Line2D([0], [0], color='orange', lw=elinewidth, label='warm'),
                plt.Line2D([0], [0], color='blue', lw=elinewidth, label='cold')
            ]

            ASICvalueslegend = False
            if ASICvalueslegend:
                ASICvaluelegend =  '{:.2f}$\pm${:.2f}'.format(KASIC, uKASIC)
                ASICgoodvaluelegend =  '{:.2f}$\pm${:.2f}'.format(KASICgood, uKASICgood)
            else:
                ASICvaluelegend, ASICgoodvaluelegend = '', ''

            if ASIC_opt:
                custom_handles.append(plt.Line2D([0], [0], color='red', lw=elinewidth, label= 'ASIC ' + ASICvaluelegend))

            if ASICgood_opt:
                custom_handles.append(plt.Line2D([0], [0], color='purple', lw=elinewidth, label='ASICgood ' + ASICgoodvaluelegend))

            if cut_x_axis_opt:
                ax2.legend(handles=custom_handles)
            else:
                if extra_plotting == True:
                    ax.legend(handles=custom_handles,loc='upper right')
                    x_lim_min, x_lim_max = 12.25, 16
                    step_ticks = 0.25

                    ax.set_xlim(x_lim_min, x_lim_max)
                    ax.set_xticks(np.arange(x_lim_min, x_lim_max+step_ticks, step_ticks))
                else:
                    ax.legend()
                #ax.legend(handles = custom_handles, loc='upper right')

            

        if typeplot =='pixels':
            data_in_x_axis = False
            if data_in_x_axis:
                save_fig_label = '_'+ typeplot

                substring = '/N037_new'

                overall_path_to_datasets = results_collection_path.split(substring)[0] + substring
            
                log_file_string = 'log_file_'

                fig, ax1 = plt.subplots(1, 1, sharey=True, figsize=(10, 5), dpi=300)
                ax = [ax1]
                
                #creating lists to append value while looping
                all_pixels = []
                all_percentages = []

                for j in range(len(results_selected)):
                    
                    specified_result = results_selected[j]
                    result_path = overall_path_to_datasets + '/'+ specified_result
                    print("Plotting " + specified_result)
                    
                    log_file_path = result_path+'/'+log_file_string + specified_result
                    
                    
                    log_files_paths = self.get_csv_files(prefix=log_file_string, path = result_path)

                    for file in log_files_paths:
                        print("Reading data from " + file)
                        label = file.split('.')[0].split('_')[-1]

                        df = pd.read_csv(file, header=None)
                        
                        # Assuming the structure you provided is exact, we can locate the rows manually
                        fittype_row = df.loc[df[0] == 'good']
                        pixels = int(fittype_row[1].values[0])
                        percentage = float(fittype_row[2].values[0])

                        all_pixels.append(pixels)
                        all_percentages.append(percentage)

                        #print(f'Pixels: {pixels}')
                        #print(f'Percentage: {percentage}')

                        width=0.25
                        if label == 'warm':
                            color = 'orange'
                            alpha = 0.7
                            x_pos = j - width/2  # Shift warm bar to the left
                        elif label == 'cold':
                            color = '#add8e6'
                            alpha = 0.7
                            x_pos = j + width/2 # Shift cold bar to the right
                        
                        ax1.bar(x_pos, percentage, width=width, alpha=alpha, color=color,edgecolor='black', linewidth=0.5)
                        ax1.text(x_pos, percentage - 0.5, f'{pixels}', ha='center', va='bottom', fontsize=8, color='black')

                        if specified_result == 'resultscalibration':
                            ax1.axhline(y=percentage, color=color, alpha = 1,linestyle='--', linewidth=2)
        
                ax1.set_xlim(-0.5, len(results_selected) - 0.5)
                ax1.set_xticks(np.arange(len(results_selected)))
                ax1.set_xticklabels(yaxis_results_selected)
                ax1.set_ylim(min(all_percentages)-1,max(all_percentages)+1)
                ax1.grid()

                ax1.set_xlabel('Dataset ID')
                ax1.set_ylabel('Percentage of good pixels w.r.t. total pixels [%]')

                custom_handles = [
                    mpatches.Patch(color='orange', alpha=alpha, label='warm'),
                    mpatches.Patch(color='#add8e6', alpha=alpha, label='cold')
                ]
                ax1.legend(loc = 'upper right',handles = custom_handles)
            
            else:
                save_fig_label = '_' + typeplot
                
                substring = '/N037_new'
                overall_path_to_datasets = results_collection_path.split(substring)[0] + substring
                
                log_file_string = 'log_file_'
                
                fig, ax1 = plt.subplots(1, 1, sharex=True, figsize=(10, 5), dpi=300)
                ax = [ax1]
                
                # Creating lists to append values while looping
                all_pixels = []
                all_percentages = []
                
                for j in range(len(results_selected)):
                    specified_result = results_selected[j]
                    result_path = overall_path_to_datasets + '/' + specified_result
                    print("Plotting " + specified_result)
                    
                    log_file_path = result_path + '/' + log_file_string + specified_result
                    
                    log_files_paths = self.get_csv_files(prefix=log_file_string, path=result_path)
                
                    for file in log_files_paths:
                        print("Reading data from " + file)
                        label = file.split('.')[0].split('_')[-1]
                
                        df = pd.read_csv(file, header=None)
                        
                        # Assuming the structure you provided is exact, we can locate the rows manually
                        fittype_row = df.loc[df[0] == 'good']
                        pixels = int(fittype_row[1].values[0])
                        percentage = float(fittype_row[2].values[0])
                
                        all_pixels.append(pixels)
                        all_percentages.append(percentage)
                
                        # print(f'Pixels: {pixels}')
                        # print(f'Percentage: {percentage}')
                
                        height = 0.25
                        if label == 'warm':
                            color = 'orange'
                            alpha = 0.7
                            y_pos = j - height / 2  # Shift warm bar down
                        elif label == 'cold':
                            color = '#add8e6'
                            alpha = 0.7
                            y_pos = j + height / 2  # Shift cold bar up
                        
                        ax1.barh(y_pos, percentage, height=height, alpha=alpha, color=color, edgecolor='black', linewidth=0.5)
                        ax1.text(percentage - 1, y_pos, f'{pixels}', va='center', ha='left', fontsize=8, color='black')
                
                        if specified_result == 'resultscalibration':
                            ax1.axvline(x=percentage, color=color, alpha=1, linestyle='--', linewidth=2)
                
                ax1.set_ylim(-0.5, len(results_selected) - 0.5)
                ax1.set_yticks(np.arange(len(results_selected)))
                ax1.set_yticklabels(yaxis_results_selected)
                ax1.set_xlim(min(all_percentages) - 1, max(all_percentages) + 1.5)
                ax1.set_xticks(np.arange(int(min(all_percentages)) - 1, int(max(all_percentages)) + 1.5, 1))
                ax1.grid()
                
                ax1.set_ylabel('Dataset ID')
                ax1.set_xlabel('Percentage of good pixels w.r.t. total pixels [%]')
                
                custom_handles = [
                    mpatches.Patch(color='orange', alpha=alpha, label='warm'),
                    mpatches.Patch(color='#add8e6', alpha=alpha, label='cold')
                ]
                ax1.legend(loc='upper right', handles=custom_handles)
            

        plt.tight_layout()

        if len(results_selected) == 1:
            save_fig_label = save_fig_label + '_'+results_selected[0]

        plt.savefig(savepath + '/results_summary'+save_fig_label+'.png')
        print('Plot saved to {}'.format(savepath + '/results_summary'+save_fig_label+'.png'))

    def plot_energy_vs_target(self, results_selected=['calibration_baseline'], results_collection_path='/data/bfys/apuicerc/N037_new/results_collection.csv', typeplot = 'gain',savepath='/data/bfys/apuicerc/N037_new/plots',  E_xrays = 5900, E_ehp = 3.69, E_ehp_err = 0.11):
        
        datasets_to_plot = 'Main'
        ASIC_name = 'VP3-1'

        #K_total_uncertainty = K_matrix * np.sqrt((n_ehp_err / n_ehp) ** 2 + (utarget_matrix / target_matrix) ** 2)
        
        
        n_ehp = E_xrays / E_ehp
        n_ehp_err = n_ehp * E_ehp_err / E_ehp

        #'resultscalibration_separation/even', 'resultscalibration_separation/odd', 'resultscalibration_separation/rows'
        #'resultscalibration_separation/even', 'resultscalibration_separation/odd','resultscalibration_separationspp/even', 'resultscalibration_separationspp/odd',
        # 'calibration_baseline','resultsultimate_newbaseline', 'resultscalibrationfitfunc4withterm0'
        
        #yaxis_results_selected = ['calibration', 'even', 'odd', 'rows']
        #yaxis_results_selected = ['calibration', 'even', 'odd', 'evenspp', 'oddspp']

        ## FOR BASELINE INCLUSION VARIATION
        #['resultscalibration', 'calibration_baseline','resultsultimate_newbaseline', 'resultscalibration_newbaseline']
        #yaxis_results_selected = ['calibration', '_baseline', '_newbaseline'] #'fitfunc4withterm0'
        if datasets_to_plot == 'Main':
            results_selected=['calibration_baseline']
            yaxis_results_selected = ['calibration']
            plotting_opt = [True, False, False] #corresponding to value to be plotted, mean, ASIC, ASICgood

        remove_word_from_labels = False
        if remove_word_from_labels:
            # Define a function to handle the special cases and general case
            def remove_results_prefix(s, index, keyword_to_remove):
                if index == 0:
                    return 'calibration'
                elif index == 1:
                    return 'calibrationnoterm0'
                elif s.startswith(keyword_to_remove):
                    return s[len(keyword_to_remove):]
                return s

            # Remove the prefix 'results' from each string in the list
            keyword_to_remove = 'resultscalibration'

            # Use list comprehension with enumeration to apply the function to each element
            cleaned_results_for_label = [remove_results_prefix(s, i, keyword_to_remove) for i, s in enumerate(results_selected)]

        
            yaxis_results_selected = cleaned_results_for_label
        else:
            yaxis_results_selected = yaxis_results_selected
        
        data_storage = {
            'targetmean_list': {
                'warm': [],
                'cold': []
            },
            'utargetmean_list': {
                'warm': [],
                'cold': []
            },
            'targetASIC_list': {
                'warm': [],
                'cold': []
            },
            'utargetASIC_list': {
                'warm': [],
                'cold': []
            },
            'Kmean_list': {
                'warm': [],
                'cold': []
            },
            'uKmean_list': {
                'warm': [],
                'cold': []
            }
            ,
            'KASIC_list': {
                'warm': [],
                'cold': []
            }
            ,
            'uKASIC_list': {
                'warm': [],
                'cold': []
            },
            'KASICgood_list': {
                'warm': [],
                'cold': []
            }
            ,
            'uKASICgood_list': {
                'warm': [],
                'cold': []
            }
        }
        
        label_l = ['warm', 'cold']
        
        if typeplot == 'gain':
            print("Reading data from " + results_collection_path)
            data = pd.read_csv(results_collection_path)

                    
            mean_opt, ASIC_opt, ASICgood_opt = plotting_opt

            
            
            #Set the below variable to True if you want to have a cut in the x axis
            #Convinient if there is a large difference in value for the gain
            cut_x_axis_opt=False

            if cut_x_axis_opt:
                fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(15, 5), dpi=300)
                # Adjust the spacing between the subplots
                fig.subplots_adjust(wspace=0.025)

                x1_lim_min, x1_lim_max = 13, 13.5#13.3
                x2_lim_min, x2_lim_max = 14.1, 14.5
                step_ticks = 0.1

                ax1.set_xlim( x1_lim_min, x1_lim_max)
                ax2.set_xlim(x2_lim_min, x2_lim_max)
                ax1.set_xticks(np.arange(x1_lim_min, x1_lim_max, step_ticks))
                ax2.set_xticks(np.arange(x2_lim_min, x2_lim_max+step_ticks, step_ticks))
                axes = [ax1, ax2]
                
                    # Hide the right spine on the first plot and the left spine on the second plot
                ax1.spines['right'].set_visible(False)
                ax2.spines['left'].set_visible(False)
                
                # Add diagonal lines to indicate the break
                d = 0.015  # size of diagonal lines in axes coordinates
                kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
                ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)        # top-left diagonal
                ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal

                kwargs.update(transform=ax2.transAxes)  # switch to the right axes
                ax2.plot((-d, +d), (-d, +d), **kwargs)        # top-right diagonal
                ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal


            else:
                fig, ax1 = plt.subplots(1, 1, sharey=True, figsize=(8, 5), dpi=300)
                
                x_lim_min, x_lim_max = 80, 150
                step_ticks = 5

                ax1.set_xlim(x_lim_min, x_lim_max)
                ax1.set_xticks(np.arange(x_lim_min, x_lim_max+step_ticks, step_ticks))

                ax1.set_xlabel('Target Threshold [DAC]')
                ax1.set_ylabel('Energy [eV]')
                energy_ylimmin, energy_ylimmax = 1000, 10000
                ax1.set_ylim(energy_ylimmin, energy_ylimmax)
                #ax.set_yticks(np.arange(len(results_selected)))
                #ax.set_yticklabels(yaxis_results_selected)
                ax1.grid()

                # Create a secondary y-axis on the left
                ax2 = ax1.twinx()
                ax2.set_ylabel('Number of electrons n$_{e-}$')
                ax2.set_ylim(energy_ylimmin/E_ehp, energy_ylimmax/E_ehp)

                
                

                axes = [ax1]
            
            categories = ['even', 'odd', 'rows']
            save_fig_label = ''
            for i in range(len(results_selected)):
                
                specified_result = results_selected[i]
                
                if specified_result.split('_')[-1].split('/')[0] == 'separation':
                    save_fig_label = '_separation'

                    for cat in categories:
                        separation_specified_result = specified_result + '/' + cat
            
            for j in range(len(results_selected)):
                
                specified_result = results_selected[j]
                
                print("Plotting " + specified_result)
                result_params = data[data['Results'] == specified_result].drop('Results', axis=1)
                
                for i in range(len(result_params)):

                    # Get the parameters for the current result
                    label = result_params.iloc[i].values[0]
                    ST = result_params.iloc[i].values[1]
                    baseline = result_params.iloc[i].values[2]
                    result_param = [float(x) for x in result_params.iloc[i].tolist()[3:]]
                    
                    targetmean, utargetmean, targetASIC, utargetASIC, targetASICgood, utargetASICgood, Kmean, Kstd, uKmean, KASIC, uKASIC, KASICgood, uKASICgood = result_param
                    # now result_param contains all the floats from results_collection.csv
                    print('Plotting for dataset {} {}'.format(specified_result, label))

                    color = 'blue' if label.startswith('cold') else 'orange'

                    custom_separation_ASIC = 0.1
                    custom_separation_ASICgood = 0.2
                    elinewidth = 2
                    capsize = 5
                    
                    
                    alpha = 1.0
                    
                    if label.startswith('cold'):
                            data_storage['targetmean_list'][label].append(targetmean)
                            data_storage['utargetmean_list'][label].append(utargetmean)
                            data_storage['targetASIC_list'][label].append(targetASIC)
                            data_storage['utargetASIC_list'][label].append(utargetASIC)
                            data_storage['Kmean_list'][label].append(Kmean)
                            data_storage['uKmean_list'][label].append(uKmean)
                            data_storage['KASIC_list'][label].append(KASIC)
                            data_storage['uKASIC_list'][label].append(uKASIC)
                            data_storage['KASICgood_list'][label].append(KASICgood)
                            data_storage['uKASICgood_list'][label].append(uKASICgood)
                            alpha = 1
                    else:
                            data_storage['targetmean_list'][label].append(targetmean)
                            data_storage['utargetmean_list'][label].append(utargetmean)
                            data_storage['targetASIC_list'][label].append(targetASIC)
                            data_storage['utargetASIC_list'][label].append(utargetASIC)
                            data_storage['Kmean_list']['warm'].append(Kmean)
                            data_storage['uKmean_list']['warm'].append(uKmean)
                            data_storage['KASIC_list']['warm'].append(KASIC)
                            data_storage['uKASIC_list']['warm'].append(uKASIC)
                            data_storage['KASICgood_list']['warm'].append(KASICgood)
                            data_storage['uKASICgood_list']['warm'].append(uKASICgood)
                    

                    
                    
                    
                    for ax in axes:
                        if mean_opt :
                            label_mean = label + ', Mean = {:.2f}$\pm${:.2f}'.format(Kmean, uKmean)
                        elif j == 0:
                            label_mean =  label
                        else:
                            label_mean = ''
                        #ax.errorbar(targetmean, E_xrays, xerr=utargetmean, fmt='o', ecolor=color,capsize=capsize, elinewidth=elinewidth, color=color, alpha=alpha,label=label_mean )
                        
                        
                        if ASIC_opt:
                            ax.errorbar(targetASIC, E_xrays - custom_separation_ASIC, xerr=utargetASIC, fmt='o', ecolor='red',
                                        capsize=capsize, elinewidth=elinewidth, color='red',
                                        label='ASIC' if j == 0 and i == 1 else '')
                        
                        if ASICgood_opt:
                            ax.errorbar(KASICgood, E_xrays - custom_separation_ASICgood, xerr=uKASICgood, fmt='o', ecolor='purple',
                                        capsize=capsize, elinewidth=elinewidth, color='purple',
                                        label='ASICgood' if j == 0 and i == 1 else '')
            
            if len(results_selected) != 1:

                sys_uncertainty_l = []
                for label in label_l:
                    sys_uncertainty = self.systematic_uncertainty_finder(results_selected,data_storage['Kmean_list'][label], data_storage['uKmean_list'][label])
                    sys_uncertainty_l.append(sys_uncertainty)

               
            else:
                sys_uncertainty_l = self.systematic_uncertainty_finder(results_selected,data_storage['Kmean_list'][label], data_storage['uKmean_list'][label])
            
            print('Systematic uncertainties found for warm and cold', sys_uncertainty_l)
            

            ax.hlines(y=E_xrays, xmin=80, xmax=150, color='black', linestyle='--', linewidth=elinewidth, alpha=0.5, 
                label='Fe55, $E_\gamma$ = 5900 eV')
                
            # Adding text label near the line (if needed)
            ax.text(134, E_xrays-50, r'Fe55, $E_\gamma$ = 5900 eV', color='black', verticalalignment='top')


            extra_plotting = True
            systematic_uncertainty_opt = True
            estimate_gain_opt = True
            K_total_uncertainty_l= []
            target_total_uncertainty_l = []

            custom_handles = []
            if extra_plotting == True:
                if systematic_uncertainty_opt == True:
                    dummycounter = 0

                    for temp in label_l:
                        
                        targetmean_reference = data_storage['targetmean_list'][temp][0]
                        utargetmean_reference = data_storage['utargetmean_list'][temp][0]

                        Kmean_reference = data_storage['Kmean_list'][temp][0]
                        uKmean_reference = data_storage['uKmean_list'][temp][0]
                        sys_uncertainty_label = sys_uncertainty_l[dummycounter]

                    
                        K_total_uncertainty = np.sqrt(sys_uncertainty_label ** 2 + uKmean_reference ** 2)
                        K_total_uncertainty_l.append(K_total_uncertainty)

                        target_total_uncertainty =  targetmean_reference * np.sqrt((n_ehp_err / n_ehp) ** 2 + (K_total_uncertainty / Kmean_reference) ** 2)
                        target_total_uncertainty_l.append(target_total_uncertainty)

                        meanvaluelegend =  ', K[e-DAC] = {:.2f}$\pm${:.2f}$\pm${:.2f}'.format(Kmean_reference, uKmean_reference,sys_uncertainty_label)
                        color = 'blue' if temp.startswith('cold') else 'orange'
                        custom_handles.append(plt.Line2D([0], [0], color=color, lw=elinewidth, label=temp +', '+ASIC_name + meanvaluelegend))
                        
                        ax.errorbar(targetmean_reference, E_xrays, xerr=target_total_uncertainty, fmt='o', ecolor=color,
                                    capsize=capsize, elinewidth=elinewidth, color=color, alpha=alpha,
                                    label=label_mean )

                        dummycounter += 1
                        #for ax in axes:
                            #ax.fill_betweenx([E_xrays, len(results_selected)], targetmean_reference - utargetmean_reference - target_total_uncertainty, targetmean_reference + utargetmean_reference + target_total_uncertainty, color='grey', alpha=0.2, label = 'Systematic uncertainty')
                

                

                if estimate_gain_opt == True:
                    # Given values
                    DAC_mV = 0.38  # mV per DAC
                    G_paper_mV_per_ke = 24.6  # mV per ke⁻
                    gain_variation_percentage = 3.3/100
                    # Convert the gain from mV/ke⁻ to mV/e⁻
                    G_paper_mV_per_e = G_paper_mV_per_ke / 1000  # because 1 ke⁻ = 1000 e⁻
                    
                    # Calculate the number of electrons per DAC
                    K_paper = DAC_mV / G_paper_mV_per_e
                    target_paper = n_ehp/K_paper

                    estimate_std = K_paper * gain_variation_percentage
                    estimate_target_uncertainty =  target_paper * np.sqrt((n_ehp_err / n_ehp) ** 2 + (estimate_std / K_paper) ** 2)
                    tag = 'Estimate paper'

                    ax.axvline(x=target_paper, ymin=0, ymax=15000, color='green', linestyle='-', linewidth=elinewidth, alpha=alpha, label='Estimate K[e-/DAC] = {:.2f}$\pm${:.2f}'.format(K_paper, estimate_std))
                    ax.fill_betweenx([0, 15000], target_paper - estimate_target_uncertainty, target_paper + estimate_target_uncertainty, color='green', alpha=0.2)

                    custom_handles.append(plt.Line2D([0], [0], color='green', alpha=alpha,
                                    label='Estimate K[e-/DAC] = ' + '{:.2f}$\pm${:.2f}'.format(K_paper, estimate_std)))

                    # Calculate the difference between K_paper and Kmean
                    Kmean_reference = data_storage['Kmean_list']['cold'][0]
                    uKmean_reference = data_storage['uKmean_list']['cold'][0]
                    sys_uncertainty_reference = sys_uncertainty_l[1]

                    def sigma_difference_comparer(Kmean_reference, uKmean_reference, sys_uncertainty_reference, K_compare, uK_compare,tag = None):
                        print('Comparing {} dataset to {} ...'.format(label, tag))
                        difference = K_compare - Kmean_reference
                        
                        # Here calculating the total uncertainty associated to the data point: including the error of the mean and systematic uncertainty using quadrature
                        total_error_K_data = np.sqrt(sys_uncertainty_reference**2 + uKmean_reference**2)
                        #total_error_K_data = uKmean_reference
                        # Combine uncertainties
                        combined_uncertainty = np.sqrt(uK_compare**2 + total_error_K_data**2)

                        # Calculate number of standard deviations
                        sigma_difference = difference / combined_uncertainty

                        print('Distance in standard deviations from {}: {:.2f}'.format(tag, sigma_difference))


                    sigma_difference_comparer(Kmean_reference, uKmean_reference, sys_uncertainty_reference, K_paper, estimate_std, tag)
                    
                    temp_reference = 'cold'
                    Kmean_reference = data_storage['Kmean_list'][temp_reference][0]
                    uKmean_reference = data_storage['uKmean_list'][temp_reference][0]
                    values_to_compare = ['KASIC','KASICgood']
                    for compare in values_to_compare:
                        K_compare = data_storage[compare+'_list'][temp_reference][0]
                        uK_compare = data_storage['u'+compare+'_list'][temp_reference][0]
                    
                        sigma_difference_comparer(Kmean_reference, uKmean_reference, sys_uncertainty_reference, K_compare, uK_compare, tag=compare)
            

            meanvalueslegend = False
            if meanvalueslegend:
                meanvaluelegend =  '{:.2f}$\pm${:.2f}'.format(Kmean, uKmean)
            else:
                meanvaluelegend = ''

            if extra_plotting == False:
                custom_handles = [
                plt.Line2D([0], [0], color='orange', lw=elinewidth, label='warm'),
                plt.Line2D([0], [0], color='blue', lw=elinewidth, label='cold')
            ]

            ASICvalueslegend = False
            if ASICvalueslegend:
                ASICvaluelegend =  '{:.2f}$\pm${:.2f}'.format(KASIC, uKASIC)
                ASICgoodvaluelegend =  '{:.2f}$\pm${:.2f}'.format(KASICgood, uKASICgood)
            else:
                ASICvaluelegend, ASICgoodvaluelegend = '', ''

            if ASIC_opt:
                custom_handles.append(plt.Line2D([0], [0], color='red', lw=elinewidth, label= 'ASIC ' + ASICvaluelegend))

            if ASICgood_opt:
                custom_handles.append(plt.Line2D([0], [0], color='purple', lw=elinewidth, label='ASICgood ' + ASICgoodvaluelegend))

            if cut_x_axis_opt:
                ax2.legend(handles=custom_handles)
            else:
                if extra_plotting == True:
                    ax.legend(handles=custom_handles,loc='upper right', fontsize = 'medium')
               
                else:
                    ax.legend(fontsize = 'medium')
                #ax.legend(handles = custom_handles, loc='upper right')

            

        ax.set_title('Energy of incident particle vs Target')
        

        plt.tight_layout()

        if len(results_selected) == 1:
            save_fig_label = save_fig_label + '_'+results_selected[0]

        plt.savefig(savepath + '/results_summary_energy_vs_target'+save_fig_label+'.png')
        print('Plot saved to {}'.format(savepath + '/results_summary_energy_vs_target'+save_fig_label+'.png'))
        



    def plot_collection(self, specified_result, results_collection_path='/data/bfys/apuicerc/N037_new/results_collection.csv', specifics=[True, True, True], opt = None, temp = None):

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5), dpi=300)

        axs = [ax1, ax2]
        
        ## set the below variable to True for performing count difference between datasets in separation of data
        ## as it is coded below it normalises the counts per bin by dividing by the total amount of good pixels
        ## then the counts per bin can be subtracted, in  this case with even and odd
        ## inorder for the code to work properly and give plot with normalised difference in coutn, -t cold or warm needs to be specifed
        difference_counts_opt = True

        if difference_counts_opt == True:
            if temp == None:
                print('Error found trying to plot difference in counts between even and odd columns...')
                exit('Need to specify temp (-t cold, -t warm) for the code plot to work')
            categories = ['even', 'odd']
        else:
            categories = ['even', 'odd', 'rows']


        broken_specified_result = specified_result.split('_')[-1]
        print('broken_specified_result in plot_collection plotcaller', broken_specified_result)
        if 'separation' in broken_specified_result:
            ## THE LINE SOF COD EBELOW ARE FOR STORAGE AND COMPARSION BETWEEN EVEN AND ODD
            even_data = {'target': None, 'edac': None}
            odd_data = {'target': None, 'edac': None}
            
            for cat in categories:
                

                cat_specified_result = specified_result + '/'+cat 
                print('Plotting collection histograms...')
                print('specified_result', cat_specified_result, 'opt', opt, 'results_collection_path', results_collection_path, 'specifics', specifics)
                
                output = self.results_collection(axs ,specified_result=cat_specified_result, results_collection_path=results_collection_path, specifics=[True, False, False], opt = opt, cat = cat, temp = temp, difference_counts_opt = difference_counts_opt)

                if output is not None:
                    target_matrix  = output[0]
                    edac_matrix = output[1]
                    if cat == 'even':
                        even_data['target'] = target_matrix
                        even_data['edac'] = edac_matrix
                    if cat == 'odd':
                        odd_data['target'] = target_matrix
                        odd_data['edac'] = edac_matrix
            
            if even_data['target'] is not None and odd_data['target'] is not None:
                n = 51
                bin_edges_target = np.linspace(70, 150, n)
                bin_edges_K = np.linspace(10, 20, n) 
                
                # Compute histograms for even and odd categories
                even_target_counts, _ = np.histogram(even_data['target'], bins=bin_edges_target)
                odd_target_counts, _ = np.histogram(odd_data['target'], bins=bin_edges_target)

                even_target_counts = even_target_counts/len(even_data['target'])
                odd_target_counts = odd_target_counts/len(odd_data['target'])

                target_count_diff = even_target_counts - odd_target_counts
                

                even_edac_counts, _ = np.histogram(even_data['edac'], bins=bin_edges_K)
                odd_edac_counts, _ = np.histogram(odd_data['edac'], bins=bin_edges_K)
                
                even_edac_counts = even_edac_counts/len(even_data['edac'])
                odd_edac_counts = odd_edac_counts/len(odd_data['edac'])

                edac_count_diff = even_edac_counts - odd_edac_counts

                print('sum', sum(even_data['target'])- sum(odd_data['target']))

                def above_zero_check(counts, bins):
                    values = []
                    
                    for i in len(range(counts)):
                        value = counts[i]*bins[i]
                        values.append(value)

                    return sum(values)


                def max_value_comparision(counts, extra_offset = True):
                    # small fucntion to get the maximum or minimim counts diff value to rescale the y axis limits
                    counts_max = abs(max(counts))
                    counts_min =  abs(min(counts))
                    final_counts_lim = 0
                    if counts_max >= counts_min:
                        final_counts_lim = counts_max
                    else:
                        final_counts_lim = counts_min
                    
                    if extra_offset != None:
                        #addign extra off set so that the y axis doesnt exatly match the the height of the bins
                        final_counts_lim += final_counts_lim*0.1

                    return final_counts_lim


                # Plot the differences
                

                bin_centers_target = (bin_edges_target[:-1] + bin_edges_target[1:]) / 2
                ax1.bar(bin_centers_target, target_count_diff, width=bin_edges_target[1] - bin_edges_target[0], color='purple', edgecolor='black', linewidth=0.5)
                #ax1.plot([], [], ' ', label='{}'.format(above_zero_check(target_count_diff, bin_centers_target)))
                ax1.set_title('Difference in Target Counts (Even - Odd)')
                ax1.set_xlabel('Target [DAC]')
                ax1.set_ylabel('Count Difference')
                ax1.set_ylim(-max_value_comparision(target_count_diff), max_value_comparision(target_count_diff) )
                ax1.legend()

                bin_centers_edac = (bin_edges_K[:-1] + bin_edges_K[1:]) / 2
                ax2.bar(bin_centers_edac, edac_count_diff, width=bin_edges_K[1] - bin_edges_K[0], color='purple', edgecolor='black', linewidth=0.5)
                #ax2.plot([], [], ' ', label='{}'.format(above_zero_check(edac_count_diff, bin_centers_edac)))
                ax2.set_title('Difference in eDAC Counts (Even - Odd)')
                ax2.set_xlabel('K [-e/DAC]')
                ax2.set_ylabel('Count Difference')
                ax2.set_ylim(-max_value_comparision(edac_count_diff), max_value_comparision(edac_count_diff) )
                ax2.legend()

        else:
            print('Plotting collection histograms...')
            print('specified_result ?no separation?', specified_result, 'opt', opt, 'results_collection_path', results_collection_path, 'specifics', specifics, 'opt', opt, 'temp', temp)
        
            self.results_collection(axs,specified_result=specified_result, results_collection_path=results_collection_path, specifics=specifics, opt = opt, temp = temp)


        categories = ['even', 'odd', 'rows']
        specified_result_fixed = specified_result.split('/')[-1]
        print('specified_result_fixed', specified_result_fixed)
            
            
        if specified_result_fixed in categories:
            specified_result = specified_result.split('/')[-2] +specified_result.split('/')[-1]
        else:
            specified_result = specified_result.split('/')[-1]

        plt.tight_layout()
        if opt == 'discrepancy':
            fig.suptitle('Target and eDAC Discrepancy distributions for {} w.r.t {} \n'.format(specified_result, opt))
            
            plt.savefig('/data/bfys/apuicerc/N037_new/plots/' + specified_result +'_'+opt+ '_discrepancy_target_eDAC_histograms_warm_cold.png')
            print('Figure saved at {}'.format('/data/bfys/apuicerc/N037_new/plots/' + specified_result +'_'+opt+'_discrepancy_target_eDAC_histograms_warm_cold.png'))
        elif opt == 'QQ':
            #fig.suptitle('Target and eDAC Discrepancy distributions for {} w.r.t {} \n'.format(specified_result, opt))
            
            plt.savefig('/data/bfys/apuicerc/N037_new/plots/' + specified_result +'_'+opt+ '_target_eDAC_histograms_warm_cold.png')
            print('Figure saved at {}'.format('/data/bfys/apuicerc/N037_new/plots/' + specified_result +'_'+opt+'_target_eDAC_histograms_warm_cold.png'))
        else:
            
            if difference_counts_opt == True:
                fig.suptitle('{} {}'.format(specified_result, temp))
                extra_label = 'countsdiff_'
            else:
                fig.suptitle('Distributions for {}\n'.format(specified_result))
                extra_label = ''

            if temp != None:
                end_label = extra_label + temp #here just adding the temp that was plotted if specified when executing the code.
            else:
                # here adding both warm and cold since the code will plot and go over them both when plotting the data since both of them will be saved in results_colletion.csv file.
                end_label = extra_label +'warm_cold'

            plt.savefig('/data/bfys/apuicerc/N037_new/plots/' + specified_result + '_target_eDAC_histograms_'+end_label+'.png')
            print('Figure saved at {}'.format('/data/bfys/apuicerc/N037_new/plots/' + specified_result + '_target_eDAC_histograms_'+end_label+'.png'))
    
    def ind_pixel_plotter(self,pixel_row, pixel_column, radii, baseline_opt = None,FPTfolder_l = None, ST_l = None,label_l = None, p0 = None, ASIC_information = None, ASICgood_information = None, together ='False'):
    
        print('Running script for pixels around pixel {0} x {1} in {2} radii'.format(pixel_row, pixel_column, radii, baseline_opt = baseline_opt))
        if pixel_row or pixel_column is not None:

            rows, columns = self.fit_coincidence_locator(FPTfolder_l, ST_l,label_l,pixel_row, pixel_column, radii, baseline_opt = baseline_opt)
            print('Rows of pixels with good fit in area considered {}'.format(rows))
            print('Columns of pixels with good fit in area considered {}'.format(columns))
        else:
            print('No individual pixels will be plotted')
            rows = pixel_row
            columns = pixel_column

        if ASIC_information == None:
            print('No ASIC information given')
        else:
            print('ASIC information given')

        if ASICgood_information == None:
            print('No ASICgood information given')
        else:
            print('ASICgood information given')


        fname = 'flux_over_thresholdscan_'+FPTfolder_l.split('/')[-2]+'ST'+ST_l+'_'+str(label_l)
     
        output_path = FPTfolder_l + 'plots/ind_pixel_plots/' + fname + '.png'
        print('Creating ' + output_path + ' ...')

        
        residuals_pulls_question = input("Do you want to plot residuals, pulls, both or neither? [r/p/y/n]: ")
        residuals = None
        pulls = None

        if residuals_pulls_question.lower() == 'r':
            residuals = True
            pulls = None
            number_subplots = 2 #one for flux, another for residuals and another for pulls
        elif residuals_pulls_question.lower() == 'p':
            residuals = None
            pulls = True
            number_subplots = 2 #one for flux, another for residuals and another for pulls
        elif residuals_pulls_question.lower() == 'y':
            residuals = True
            pulls = True
            number_subplots = 3 #one for flux, another for residuals and another for pulls
        else:
            residuals = None
            pulls = None
            number_subplots = 1
        
        if number_subplots == 1:
            fig, axs = plt.subplots(number_subplots, 1, figsize=(7, 5*number_subplots), dpi=300)
            self.pixelscanplot_baseline(axs,rows, columns, ASIC_information= ASIC_information,ASICgood_information = ASICgood_information,FPTfolder = FPTfolder_l, ST = ST_l, p0 = p0, label = label_l,residuals = residuals, pulls = pulls)
            
        else:
            fig, axs = plt.subplots(number_subplots, 1, figsize=(15, 5*number_subplots), dpi=300)
            self.pixelscanplot_baseline(axs,rows, columns, ASIC_information= ASIC_information,ASICgood_information = ASICgood_information, FPTfolder = FPTfolder_l, ST = ST_l, p0 = p0, label = label_l,residuals = residuals, pulls = pulls)
            

        #plt.figure(figsize = (15,5))
        #self.pixelscanplot_baseline(axs,rows, columns, ASIC_information= ASIC_information,FPTfolder = FPTfolder_l, ST = ST_l, p0 = p0, label = label_l,save_fig='True')

        fig.suptitle('ST={} from {} {}'.format(ST_l, str(FPTfolder_l.split('/')[-2]), str(label_l)))
        
        
        plt.tight_layout()
        fig.savefig(output_path, format='png')
     
        plt.close('all')
        print('Figure saved in {}'.format(output_path))

        #if len(rows) > 1:
          
            # commenting out the rest of the function
            #one run separate if you want a plot with the all of the pixels together
            #if together == 'True':
            #    print('Creating plots with pixels together')
            #    for file in range(0,len(FPTfolder_l)):
            #        #self.pixelscanplot_baseline(FPTfolder_l[file], FPTfilename_l[file], ST_l[file], fromThr_l[file], toThr_l[file], stepThr_l[file],label_l[file], p0[file], rows_l, columns_l, baselineparams=baselineparams, baseline_opt=baseline_opt, alpha = 0.05, save_fig='True')
            #else:
            #    #self.pixelscanplot_baseline(FPTfolder_l, ST_l, fromThr_l, toThr_l, stepThr_l,label_l, p0, rows_l, columns_l, baselineparams=baselineparams, baseline_opt=baseline_opt, alpha = 0.05, save_fig='True')

        #if type_data == 'acqbaseline':
            #one run separate if you want a plot with the all of the pixels together
            #if together == 'True':
            #    print('Creating plots with pixels together')
            #    for file in range(2,len(FPTfolder_l)):
            #        self.pixelscanplot_baseline(FPTfolder_l[file], FPTfilename_l[file], ST_l[file], fromThr_l[file], toThr_l[file], stepThr_l[file],label_l[file], p0[file], rows_l, columns_l, baselineparams=baselineparams, baseline_opt=baseline_opt, alpha = 0.05, save_fig='True')


        if together == False:    
            #saves pixels plots individually
            print('together = False')
            for pixel in range(len(rows)):
                    row = [rows[pixel]]
                    column = [columns[pixel]]
                    
                    self.pixelscanplot_baseline(FPTfolder_l[file], FPTfilename_l[file], ST_l[file], fromThr_l[file], toThr_l[file], stepThr_l[file],label_l[file], p0[file], row, column, baselineparams=baselineparams, baseline_opt=baseline_opt, alpha = 0.05, save_fig='True')
                
            print('Common pixel between data sets plots saved. All-in-one plot = {0}'.format(together))

