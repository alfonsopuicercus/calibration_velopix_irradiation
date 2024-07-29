import glob as glob
import os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import pandas as pd
from scipy.stats import skew, kurtosis, kurtosistest, probplot, norm

from scipy.optimize import curve_fit

from utils import organizer


class Plot(organizer.FileOrganizer):
    def __init__(
        self
    ):
        """
        a class for plotting functions
        """
    def param_hist(self,ax, param,temp,opt = None, path = ''):
        data = self.open_csv(prefix=param, label = temp, path = path)[0]
        udata = self.open_csv(prefix='u'+param, label = temp, path = path)[0]

        data=np.ma.masked_where(data<=0, data)
        udata=np.ma.masked_where(udata<=0, udata)

        mean_data = np.mean(data)
        std_data = np.std(data)

        calibration_plots_share_x_and_y_axis = True
        if calibration_plots_share_x_and_y_axis:
            if param == 'E0':
                minx, maxx = 1500, 1580
                maxy = 600
                data=np.ma.masked_where(data>=maxx, data)
            if param == 'target':
                minx, maxx = 70, 160
                maxy = 600
                data=np.ma.masked_where(data>=maxx, data)
            if param == 'f':
                
                minx, maxx = 0, 0.01
                maxy = 1200
                data=np.ma.masked_where(data>=maxx, data)

            if param == 's': 
                minx, maxx = 0, 15
                maxy = 2000
                data=np.ma.masked_where(data>=maxx, data)

           
        
        if opt == None:
            number_of_pixels = len(data[data>0])
            
            values_flat = data.flatten()
    
            ax.hist(values_flat, bins=50, edgecolor='black', linewidth=0.5)
            decimal_places = 2 if param != 'f' else 5
            ax.axvline(x=np.mean(values_flat), color='orange', 
            label='Mean = {:.{dp}f} +{:.{dp}f}'.format(np.mean(values_flat), np.std(values_flat)/np.sqrt(number_of_pixels), dp=decimal_places))
    
            
            ASIC_data = pd.read_csv(os.path.join(path, 'ASIC_params'+str(temp)+'.csv'), delimiter=',', names=['poptASIC', 'Value', 'Uncertainty'])
            ASICgood_data = pd.read_csv(os.path.join(path, 'ASICgood_params'+str(temp)+'.csv'), delimiter=',', names=['poptASICgood', 'Value', 'Uncertainty'])

            parameter_rows = ASIC_data[ASIC_data['poptASIC']==param]
            if not parameter_rows.empty:
                ASIC_value_to_plot = float(parameter_rows['Value'].values[0])
                ASIC_uncertainty_to_plot = float(parameter_rows['Uncertainty'].values[0])
                ax.axvline(x=ASIC_value_to_plot, color='red', linestyle='--', label = 'ASIC = {:.{dp}f} +{:.{dp}f}'.format(ASIC_value_to_plot,ASIC_uncertainty_to_plot, dp=decimal_places))

            parameter_rows = ASICgood_data[ASICgood_data['poptASICgood']==param]
            if not parameter_rows.empty:
                ASICgood_value_to_plot = float(parameter_rows['Value'].values[0])
                ASICgood_uncertainty_to_plot = float(parameter_rows['Uncertainty'].values[0])
                ax.axvline(x=ASICgood_value_to_plot, color='purple', linestyle='--', label = 'ASICgood = {:.{dp}f} +{:.{dp}f}'.format(ASICgood_value_to_plot,ASICgood_uncertainty_to_plot, dp=decimal_places))

            ax.set_title('{} Distribution'.format(param))
            ax.set_xlabel('{} Value'.format(param))
            ax.set_ylabel('Frequency')

            legend_loc = 'upper right' if param not in ['E0','target'] else 'upper right'
            ax.legend(loc=legend_loc, prop={'size': 9})
            ax.set_xlim(minx, maxx)
            ax.set_ylim(0, maxy)

        if opt != None:
            data = data.flatten()
            udata = udata.flatten()

            # Remove masked values from data
            mask = np.ma.getmaskarray(data)
            data = data[~mask]
            udata = udata[~mask]

            ASIC_data = pd.read_csv(os.path.join(path, 'ASIC_params'+str(temp)+'.csv'), delimiter=',', names=['poptASIC', 'Value', 'Uncertainty'])
            ASICgood_data = pd.read_csv(os.path.join(path, 'ASICgood_params'+str(temp)+'.csv'), delimiter=',', names=['poptASICgood', 'Value', 'Uncertainty'])

            
            if opt == 'ASIC':
                parameter_rows = ASIC_data[ASIC_data['poptASIC']==param]
            else:
                parameter_rows = ASICgood_data[ASICgood_data['poptASICgood']==param]

            if not parameter_rows.empty:
                ASIC_value_to_plot = float(parameter_rows['Value'].values[0])
                ASIC_uncertainty_to_plot = float(parameter_rows['Uncertainty'].values[0])

            std_from_compared_value = []
            for i in range(len(data.flatten())):
                std_from_compared_value.append((data[i]-ASIC_value_to_plot)/np.sqrt(udata[i]**2  + ASIC_uncertainty_to_plot**2))

            counts, bins = np.histogram(std_from_compared_value, bins =20, range=(-10,10))
            total3sigma = 0
            for i in range(10-3,10+3): 
                total3sigma+=counts[i]

            bins_to_color = [-3,-2, -1, 0, 1, 2,3]

            ax.hist(std_from_compared_value, bins = 20, edgecolor='black', linewidth=0.5, range = (-10,10))
            ax.hist(std_from_compared_value, bins=bins_to_color, edgecolor='black', linewidth=0.5, color='green', alpha=0.7)
            ax.set_xlabel('Standard deviation $\sigma$')
            ax.set_ylabel('Frequency')
    
            ax.set_title('Discrepancy of pixels from {}'.format(opt))
            #ax.set_title("Pixels within 3$\sigma$ of {} = {}. {:.2f}% w.r.t goodpixels, {:.2f}% of total".format(opt, total3sigma, total3sigma/len(data.flatten())*100, total3sigma/(256*256)*100), y = 1.005)
            
            # Add the legend
            percent_goodpixels = total3sigma / len(data.flatten()) * 100
            percent_total = total3sigma / (256 * 256) * 100
            legend_texts = ['% of pixels within 3$\sigma$','{:.2f}% of goodpixels'.format(percent_goodpixels), '{:.2f}% of total'.format(percent_total)]

            for legend_text in legend_texts:
                ax.plot([], [], ' ', label=legend_text)
                
            ax.legend(loc = 'upper left', fontsize = 'small')
            
    def heat_map(self, param,temp, path = ''):

        data = self.open_csv(prefix=param, label = temp, path = path)[0]
      
        udata = self.open_csv(prefix='u'+param, label = temp, path = path)[0]
        
        data=np.ma.masked_where(data<=0, data)
        udata=np.ma.masked_where(data<=0, udata)
        
        mean_data = np.mean(data)
        std_data = np.std(data)
      
        if param == 'E0':
            flatten_data = data.flatten()
            flatten_data = flatten_data[flatten_data>0]
            print('{} Number of goodpixels = {} -> {:.2f} % w.r.t total, with {} mean = {}'.format(temp, len(flatten_data), len(flatten_data)/(256**2)*100, param, mean_data))

       
        #Calculating weighted average of parameters and the uncertainty of the average  #A: recheck this is the part in the paper where he goes into weights
        swdata = 0
        swudata = 0
        for row in range(256):
            for column in range(256):

                #Ignoring the masked elements (otherwise will return -- instead of a number)
                if np.ma.is_masked(data[row][column]):
                    continue
                
                #Calculating weights
                #if statements are made to eliminate uncertainties too close to 0. Otherwise, for the fourth data set,
                #average and uncertainty of the parameters will be NaN or inf. If statements may be removed for first to third data sets
                #A: recheck, quick note here how are these if statements removed for the rest of the data sets to that it can be computed
                if udata[row][column]>10**(-15):
                    wudata=1/udata[row][column]**2
                else:
                    wudata=0
                

                #Calculating sum of weighted parameter 
                swdata+=data[row][column]*wudata
              
                
                #Calculating sum of weights 
                swudata+=wudata
            
        
        #Finally, calculating weighted average parameter and its uncertainty, A. reckec i am guess
        dataavg=swdata/swudata
        udataavg=swudata**(-0.5)
        param = str(param)
        themes = ['autumn', 'summer', 'winter']
        if param == 'E0':
            theme = themes[0]
            mincoeff, maxcoeff = 0.99, 1.01
        if param == 'target':
            theme = themes[0]
            mincoeff, maxcoeff = 0.95, 1.05
        if param == 'f': 
            theme = themes[1]
            mincoeff, maxcoeff = 0.4, 1.6
        if param == 's': 
            theme = themes[2]
            mincoeff, maxcoeff = 0.7, 1.3
    
        
        CMapp=plt.get_cmap(theme)
        CMapp.set_bad([0.94, 0.94, 0.94])

        vmin=mincoeff*dataavg
        vmax=maxcoeff*dataavg

        if data is not None:
            plt.pcolormesh(data, cmap=CMapp, vmin=mincoeff*dataavg, vmax=maxcoeff*dataavg)

            ticks=np.linspace(0, 256, 5)
            plt.xticks(ticks)
            plt.yticks(ticks)
            plt.title('{} parameter, mean={:.4f}$\pm${:.4f}'.format(param, mean_data, std_data), fontsize=13)
            plt.ylabel('Pixel rows', fontsize=10)
            plt.xlabel('Pixel columns', fontsize=10)
            plt.colorbar()
    
    def fittype_heat_map(self, param ,temp, path = ''):

        data = self.open_csv(prefix=param, label = temp, path = path)[0]

        def fittype_counter(matrix):
            good_pixels = len(data[data==1])
            percentage_good_pixels = round(good_pixels/(256**2)*100, 2)

            bad_pixels = len(data[data==0])
            percentage_bad_pixels = round(bad_pixels/(256**2)*100, 2)
            fitnotfound_pixels = len(data[data==-1])
            percentage_fitnotfound_pixels = round(fitnotfound_pixels/(256**2)*100, 2)
            cut_pixels = len(data[data==-2])
            percentage_cut_pixels = round(cut_pixels/(256**2)*100 ,2)

            print('Amount of pixels per category:')
            print('good_pixels: {} or {}\%'.format(good_pixels, percentage_good_pixels))
            print('bad_pixels: {} or {}\%'.format(bad_pixels, percentage_bad_pixels))
            print('fitnotfound_pixels: {} or {}\%'.format(fitnotfound_pixels, percentage_fitnotfound_pixels))
            print('cut_pixels: {} or {}\%'.format(cut_pixels, percentage_cut_pixels))

            counts_pixels = [good_pixels,bad_pixels,fitnotfound_pixels, cut_pixels ]
            percentages = [percentage_good_pixels, percentage_bad_pixels, percentage_fitnotfound_pixels, percentage_cut_pixels]

            return counts_pixels, percentages

        counts_pixels, percentages=  fittype_counter(data)
        percentage_good_pixels, percentage_bad_pixels, percentage_fitnotfound_pixels, percentage_cut_pixels = percentages

        Cmap=mpl.colors.ListedColormap([[0.94, 0.94, 0.94], [0.5, 0.5, 0.5], [0.9, 0.9, 0.7], [0.8, 0.4, 0.4]])

        param = str(param)

        if data is not None:
            plt.pcolormesh(data, cmap=Cmap, vmin=-2, vmax=1)

            ticks=np.linspace(0, 256, 5)
            plt.xticks(ticks)
            plt.yticks(ticks)
            plt.title('{} parameter'.format(param), fontsize=13)
            plt.ylabel('Pixel rows', fontsize=10)
            plt.xlabel('Pixel columns', fontsize=10)
            cbar = plt.colorbar(ticks=(-1.64, -0.88, -0.13, 0.62))
            tick_labels = [
                f'Cut data, {percentage_cut_pixels}%',
                f'Fit not found, {percentage_fitnotfound_pixels}%',
                f'Bad fit, {percentage_bad_pixels}%',
                f'Good fit, {percentage_good_pixels}%'
            ]
            cbar.set_ticklabels(tick_labels)
            # Set the font size for the tick labels
            for label in cbar.ax.get_yticklabels():
                label.set_fontsize(8)  # Adjust the font size as needed

    def flux_heat_map(self, param ,temp, path = ''):
        #data = np.loadtxt(path+'Fluxperpixel_Module0_VP3-1_ECS_data_ST_2s621ms'+'_'+str(temp)+'_THR_'+str(1530)+'.csv', dtype=float, delimiter=',')
        data = self.open_csv(prefix=param, label=temp, path=path)[0]
        data=np.ma.masked_where(data<0, data)
        data_masked = np.ma.masked_where(data<=75, data)
        mean_data_mask = np.mean(data_masked)
        data_flat = data.flatten()
        max_hits, min_hits = max(data_flat), min(data_flat)
        
        param = str(param)

        if data is not None:
            # Normalize the data
            norm = mcolors.Normalize(vmin=min_hits, vmax=mean_data_mask)

            plt.pcolormesh(data, norm=norm, cmap='viridis')  # Using 'viridis' for better contrast

            ticks = np.linspace(0, data.shape[0], 5)
            plt.xticks(ticks)
            plt.yticks(ticks)
            plt.title('VP3-1 ASIC Flux from Fe55 X-rays, {} data'.format(temp),fontsize=13)
            plt.ylabel('Pixel rows', fontsize=10)
            plt.xlabel('Pixel columns', fontsize=10)
            cbar = plt.colorbar()
            cbar.set_label('Hits/s')
            



    def results_collection(self, axs,specified_result, results_collection_path = '', specifics = [True, True, True], opt = None, cat = None, temp = None, difference_counts_opt = None):
        """
        This function plots histograms for target and eDAC distributions based on the specified result and specifics.
        
        Args:
            results_collection_path (str): The path to the results_collection.csv file.
            specified_result (str): The specified result folder with the processed data to plot histograms for.
            specifics (list): A list of booleans indicating whether to plot specific values in the distributions. Contains 3 booleans for  mean, ASIC, and ASICgood.
        
        Returns:
            None
            Saves the plots in the specified N037_new/plots folder.
        """
        # Read the data from the results_collection.csv file into a DataFrame
        specfics_mean, specfics_ASIC, specfics_ASICgood = specifics
        print('specifics', specifics)
        categories = ['even', 'odd', 'rows']
        print("Reading data from " + results_collection_path)
        data = pd.read_csv(results_collection_path)
        print('specified_result', specified_result)

        # Get all the parameters for the specified result
        result_params = data[data['Results'] == specified_result].drop('Results', axis=1)

        line_colors = ['red', 'blue']
        color_counter = 0
        ax1, ax2 = axs

        for i in range(len(result_params)):
            

            # Get the parameters for the current result
            label = result_params.iloc[i][0]

            #here performing label check just in case we want to have a specific temperature plotted.

            if temp != None and label != temp:

                continue


            ST = result_params.iloc[i][1]
            baseline = result_params.iloc[i][2]
            result_param = [float(value) for value in result_params.iloc[i].tolist()[3:]]
            
            targetmean, utargetmean, targetASIC, utargetASIC, targetASICgood, utargetASICgood, Kmean, Kstd, uKmean, KASIC, uKASIC,KASICgood, uKASICgood = result_param
            # now result_param contains all the floats from results_collection.csv
            print('Plotting for dataset {} label {} baseline {} opt {}'.format(specified_result, label, baseline, opt))

            target_matrix_path = os.path.join('/data/bfys/apuicerc/N037_new/', specified_result, 'targetmatrix_ST_{}_{}.csv'.format(ST, label))
            target_matrix = np.loadtxt(target_matrix_path,  dtype=float, delimiter=',')
            mask = (target_matrix > 0) & (target_matrix < 1000)
            target_matrix = target_matrix[mask].flatten()
            utarget_matrix_path = os.path.join('/data/bfys/apuicerc/N037_new/', specified_result, 'utargetmatrix_ST_{}_{}.csv'.format(ST, label))
            utarget_matrix = np.loadtxt(utarget_matrix_path,  dtype=float, delimiter=',')
            utarget_matrix = utarget_matrix[mask].flatten()


            edac_matrix_path = os.path.join('/data/bfys/apuicerc/N037_new/', specified_result, 'K_matrix_ST_{}_{}.csv'.format(ST, label))
            edac_matrix = np.loadtxt(edac_matrix_path,  dtype=float, delimiter=',')
                   
            mask = (edac_matrix > 0) & (edac_matrix < 100)
            edac_matrix = edac_matrix[mask].flatten()
            uedac_matrix_path = os.path.join('/data/bfys/apuicerc/N037_new/', specified_result, 'uK_matrix_ST_{}_{}.csv'.format(ST, label))
            uedac_matrix = np.loadtxt(uedac_matrix_path,  dtype=float, delimiter=',')
            uedac_matrix = uedac_matrix[mask].flatten()

            n_pixels = len(target_matrix[target_matrix>0])

            if label == 'warm':
                color = 'orange'
                alpha = 0.7
            elif label == 'cold':
                color = '#add8e6'
                alpha = 0.7
         
            if opt == None:
                if cat not in categories:
                    opt_label = ''
                else: 
                    ## HERE JUST ADDING THE LINE BELOW TO SKIP OVER ROWS SINC EI WANT A GRAOH TO COMPARE EVEN AND ODD
                    #if cat =='rows':
                    #    continue


                    opt_label = cat
                    if temp != None:
                        #here adding some plotting specifics to plot neatly the different groups considered if only a specific temperature is plotted.

                        alpha = 0.5
                        if cat == 'even':
                            color = 'green'
                            
                        elif cat == 'odd':
                            color = 'red'
        

                        elif cat == 'rows':
                            color = 'blue'

                skewness_target = skew(target_matrix)
                kurtosis_target = kurtosis(target_matrix)
                skewness_K = skew(edac_matrix)
                kurtosis_K = kurtosis(edac_matrix)
                
                print('number of pixels {} {:.2f}%'.format(len(target_matrix), len(target_matrix)/(256**2)*100))

                n = 51
                bin_edges_target = np.linspace(70, 150, n)
                bin_edges_K = np.linspace(10, 20, n) 

                

                if difference_counts_opt == True and cat != 'rows':
                    output = [target_matrix, edac_matrix]
              
                    return output           
                    

                                           
                else:
                    ax1.hist(target_matrix, bins = bin_edges_target, alpha = alpha, color = color, edgecolor='black', linewidth=0.5, label = '{} {}'.format(label,opt_label))
                    ax2.hist(edac_matrix, bins = bin_edges_K, alpha = alpha ,color = color, edgecolor='black', linewidth=0.5, label = '{} {}'.format(label,opt_label))
                    
                line_styles = ['-', '--', ':']
             

                if specfics_mean and targetmean != 0:
                    ax1.axvline(x = targetmean, color = line_colors[color_counter], linestyle = line_styles[0], label = 'Mean = {:.2f} +- {:.2f} '.format(float(targetmean), float(utargetmean)/np.sqrt(n_pixels)))
              
                if specfics_ASIC and targetASIC != 0:
                    ax1.axvline(x = targetASIC, color = line_colors[color_counter], linestyle = line_styles[1], label = 'ASIC = {:.2f} +- {:.2f}'.format(float(targetASIC), float(utargetASIC)))  ##np.mean(baseline_width[baseline_width>0])
       
                if specfics_ASICgood and targetASICgood != 0 :
                    ax1.axvline(x = targetASICgood, color = line_colors[color_counter], linestyle = line_styles[2], label = 'ASICgood = {:.2f} +- {:.2f}'.format(float(targetASICgood), float(utargetASICgood)))

                
            
                if specfics_mean and Kmean != 0:
                    ax2.axvline(x = Kmean, color = line_colors[color_counter], linestyle = line_styles[0], label = 'Mean = {:.2f} +- {:.2f} '.format(float(Kmean), float(uKmean)))
      
                if specfics_ASIC and KASIC != 0:
                    ax2.axvline(x = KASIC, color = line_colors[color_counter], linestyle = line_styles[1], label = 'ASIC = {:.2f} +- {:.2f}'.format(float(KASIC), float(uKASIC)))  ##np.mean(baseline_width[baseline_width>0])
        
                if specfics_ASICgood and KASICgood != 0 :
                    ax2.axvline(x = KASICgood, color = line_colors[color_counter], linestyle = line_styles[2], label = 'ASICgood = {:.2f} +- {:.2f}'.format(float(KASICgood), float(uKASICgood)))
         

                ax1.set_title('Target')
                ax2.set_title('eDAC')
                ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
                ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
                ax1.set_xlabel('Target [DAC]')
                ax2.set_xlabel('K [-e/DAC]')
                ax1.set_ylabel('Frequency')
        
                color_counter += 1
            
                

           
            elif opt == 'ASIC' or opt == 'ASICgood':
                
                axs = [ax1, ax2]
                axs_correspondance = ['Target', 'K[e-/DAC]']
                axs_data_correspondance = [target_matrix, edac_matrix] # correspondance
                axs_udata_correspondance = [utarget_matrix, uedac_matrix] # correspondance
                axs_ASIC_correspondance = [[targetASIC, targetASICgood], [KASIC, KASICgood]]
                for ax in axs:
                    print('label', label, 'axs.index(ax)', axs.index(ax), axs_correspondance[axs.index(ax)])
                    data = axs_data_correspondance[axs.index(ax)]
                    udata = axs_udata_correspondance[axs.index(ax)]
                    if opt == 'ASIC':
                        ASIC_value_to_plot = axs_ASIC_correspondance[axs.index(ax)][0]
                    else:
                        ASIC_value_to_plot = axs_ASIC_correspondance[axs.index(ax)][1]

                    std_from_compared_value = []
                    for i in range(len(data)):

                        std_from_compared_value.append((data[i] - ASIC_value_to_plot) / udata[i])

                    counts, bins = np.histogram(std_from_compared_value, bins=20, range=(-10, 10))
                    total3sigma = sum(counts[10 - 3:10 + 3])

                    bins_to_color = [-3, -2, -1, 0, 1, 2, 3]

                    ax.hist(std_from_compared_value, bins=20, range=(-10, 10), color = color, alpha = 0.5, label = '{}, pixels within 3$\sigma$ = {}. {:.2f}% w.r.t goodpixels, {:.2f}% of total'
                                .format(label, total3sigma, total3sigma / len(data.flatten()) * 100, total3sigma / (256 * 256) * 100))
                    ax.hist(std_from_compared_value, bins=bins_to_color,range=(-10, 10),  edgecolor='black', color = color, alpha = 0.7, linewidth=0.5)
                    #ax.hist(std_from_compared_value, bins=bins_to_color, edgecolor='black', linewidth=0.5, color='green', alpha=0.7)
                    ax.set_xlabel('Standard deviation $\sigma$')
                    ax.set_ylabel('Frequency')
                    ax.set_xlim(-10, 10)
                    ax.set_xticks(np.arange(-10, 11, 1))
                    
                    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.13),
                              )
                    #ax.set_title('Discrepancy of {} distribution'.format(param))
                    ax.set_title(axs_correspondance[axs.index(ax)])

            elif opt == 'QQ':
                
                if cat not in categories:
                    opt_label = ''
                else: opt_label = cat

                ax1, ax2 = axs
                
                def QQ( ax, data, data_uncertainties, label, unc = None, title = '', centered = None):
                    if centered != None:
                        data = data - np.mean(data)
                    res = probplot(data, dist="norm")
                    
                    if unc != None:
                        ax.errorbar(res[0][0], res[0][1], yerr=data_uncertainties, color=color, alpha=0.1)
                    
                    ax.scatter(res[0][0], res[0][1], color=color,  s=10, alpha=0.5, label = label)
                    # Plot the red lines
                    ax.plot(res[0][0], res[1][0]*res[0][0] + res[1][1], color=line_colors[i], linestyle='--', label = 'Expected from {}'.format(label))
                    
                    # Add labels and title
                    ax.set_xlabel('Theoretical Quantiles')
                    ax.set_ylabel('Sample Quantiles')
            
                    ax.set_title('QQ Plot {}'.format(title))
                    ax.legend()

                if cat != None:
                    totallabel = label + str(cat)
                    print('cat', cat)
                else:
                    totallabel = label
                QQ(ax1,target_matrix, utarget_matrix, totallabel, unc = None, title = 'target', centered = None)
                QQ(ax2,edac_matrix, uedac_matrix, totallabel, unc= None, title = 'K')
                #QQ(ax2,edac_matrix, uedac_matrix, label, unc = None, title = 'K')

        
    def fit_coincidence_locator(self,FPTfolder,ST, label, x,y, r,baseline_opt = False, cond='False'):
    
        rows_full_l = []
        columns_full_l = []
        lengths_full_l = []

        rows = []
        columns = []
        
        if r == 0:
            n_pixels_plot = 1
        else:
            n_pixels_plot = 5 #numebr of pixel desired to plot, have to add this t the input of the function so that it can be chosen by the user when executing the code.
        
        

        if type(FPTfolder) == str:
            rows_file, columns_file = self.good_fit_locator(FPTfolder,ST, label, x,y, r, baseline_opt, cond=True)
            
            total_pixels = len(rows_file)
            indices = [int(i * total_pixels / n_pixels_plot) for i in range(n_pixels_plot)]
            selected_pixels = [(rows_file[idx], columns_file[idx]) for idx in indices]
            for idx in indices:
                rows.append(rows_file[idx])
                columns.append(columns_file[idx])
            return rows, columns
            
        else:
            for i in range(0,len(FPTfolder)):
                rows_file, columns_file = self.good_fit_locator(FPTfolder[i],ST[i], label[i], x,y, r, baseline_opt, cond=True)
                rows_full_l.append(rows_file)
                columns_full_l.append(columns_file)
                lengths_full_l.append(len(rows_file))
            
            

        # this is from an earlier part of the code in which this was used to find pixels with good fit tha where both on cold and warm datasets
        between_warm_cold_both_with_good_fit = False
        if between_warm_cold_both_with_good_fit == True:
            common_pixel_row = []
            common_pixel_column = []
            for j in range(0,len(lengths_full_l)):
                if len(rows_full_l[j]) == min(lengths_full_l):
                    index = j
                
            for entry in range(0,len(rows_full_l[index])):
                #here going every pixel in the smallest pixel with good fit data set

                for data_set in range(len(lengths_full_l)):
                    if data_set == index:
                        continue

                    for pixel in range(0,lengths_full_l[data_set]):
                        if rows_full_l[index][entry] == rows_full_l[data_set][pixel] and columns_full_l[index][entry] == columns_full_l[data_set][pixel]:
                            common_pixel_row.append(rows_full_l[index][entry])
                            common_pixel_column.append(columns_full_l[index][entry])

            return common_pixel_row, common_pixel_column

        

    def good_fit_locator(self,FPTfolder,ST, label, x,y, r,baseline_opt = False, cond='False'):
        rows_l = []
        columns_l = []
        print('goodfitlocator FPTfolder', FPTfolder)
        if baseline_opt == True:
            fittypematrix = np.loadtxt(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'_baseline.csv', dtype = float, delimiter=',')
        else:
            fittypematrix = np.loadtxt(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv', dtype = float, delimiter=',')

        print('\nFinding pixels with good fit {}'.format(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv'))
        if r != 0:
            for row in range(x-r, x+r):
                for column in range(y-r,y+r):
                    if fittypematrix[row][column]== 1:
                        rows_l.append(row)
                        columns_l.append(column)
        else:
            rows_l.append(x)
            columns_l.append(y)

        #print('Rows = ', len(rows_l))
        #print('Columns = ', len(columns_l))
        if cond == 'True':
            print(rows_l)
            print(columns_l)
        
        return rows_l, columns_l


    def pixelscanplot_baseline(self, ax,rows, columns, FPTfolder, ST, label, p0 , ASIC_information,ASICgood_information,residuals = None, pulls = None):
        
        #print('nan_cut = {}, thrsASICrange= ({},{})'.format(nan_cut, thrsASIC[0], thrsASIC[-1]))
        
        if ASIC_information != None:
            Thrs_ASIC, fluxASIC,unc_fluxASIC, poptASIC, poptuncASIC = ASIC_information
        else:
            Thrs_ASIC = np.arange(1480,1600,1)

        if ASICgood_information != None:
            Thrs_ASICgood, fluxASICgood,unc_fluxASICgood, poptASICgood, poptuncASICgood = ASICgood_information 

        x=np.linspace(Thrs_ASIC[0],Thrs_ASIC[-1], 2000) #For plotting the fits

        FPT=[]
        unc_FPT=[]

        if rows is not None:
            for thr in Thrs_ASIC:
                f=np.loadtxt(FPTfolder+'Fluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
                u=np.loadtxt(FPTfolder+'UncertaintyFluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
                FPT.append(f)
                unc_FPT.append(u)

            
            #Plotting flux of each pixel
            for pixel in range(len(rows)):
                ci=float(pixel)/float(len(rows))

                row = int(rows[pixel])
                column = int(columns[pixel])

                flux_pixel = [i[row][column] for i in FPT]
                unc_flux_pixel = [i[row][column] for i in unc_FPT]

                trimcut=np.nonzero(flux_pixel)
                thrs=Thrs_ASIC[trimcut[0][0]:trimcut[0][-1]+1]
                flux_pixel=np.trim_zeros(flux_pixel)
                unc_flux_pixel=unc_flux_pixel[trimcut[0][0]:trimcut[0][-1]+1]
            
                popt, pcov=curve_fit(self.fitfunction, thrs, flux_pixel, sigma=unc_flux_pixel, p0=p0)
                pcov=np.sqrt(np.diag(pcov))
        
                if type(ax) is not np.ndarray:
                    ax.errorbar(thrs, flux_pixel, yerr = unc_flux_pixel, marker='o', linestyle='none', color=(0.8-0.8*ci, 0.4, 0.4-0.4*ci))
                    ax.plot(x, self.fitfunction(x, *popt), linestyle='-' , linewidth=1.7, color=(1-ci, 0.4, 0.6-0.6*ci), label='Flux and fit on pixel {}x{} \n E0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(rows[pixel], columns[pixel], popt[-1], pcov[-1], popt[0], pcov[0],popt[2], pcov[2]))
                    
                    
                else:
                    ax[0].errorbar(thrs, flux_pixel, yerr = unc_flux_pixel, marker='o', linestyle='none', color=(0.8-0.8*ci, 0.4, 0.4-0.4*ci))
                    ax[0].plot(x, self.fitfunction(x, *popt), linestyle='-' , linewidth=1.7, color=(1-ci, 0.4, 0.6-0.6*ci), label='Flux and fit on pixel {}x{} \n E0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(rows[pixel], columns[pixel], popt[-1], pcov[-1], popt[0], pcov[0],popt[2], pcov[2]))
                    
                    if residuals is not None or pulls is not None:
                        residuals_pixel= flux_pixel - self.fitfunction(np.array(thrs), *popt)
                        if residuals is not None:
                            ax[1].plot(thrs, residuals_pixel, marker='o', linestyle='-', color=(1-ci, 0.4, 0.6-0.6*ci), label = 'Pixel {}x{}'.format(rows[pixel], columns[pixel]))
                    if pulls is not None:
                        if unc_flux_pixel is not None:
                            pulls_pixel = [residuals_pixel[entry]/unc_flux_pixel[entry] for entry in range(len(residuals_pixel))]
                            
                            if residuals is None:
                                ax[1].plot(thrs, pulls_pixel, marker='o', linestyle='-', color=(1-ci, 0.4, 0.6-0.6*ci), label = 'Pixel {}x{}'.format(rows[pixel], columns[pixel]))
                            else:
                                ax[2].plot(thrs, pulls_pixel, marker='o', linestyle='-', color=(1-ci, 0.4, 0.6-0.6*ci), label = 'Pixel {}x{}'.format(rows[pixel], columns[pixel]))

            if residuals  and pulls is None:
                if ASIC_information != None:
                    ax.errorbar(Thrs_ASIC, fluxASIC, marker='o',color=[0.3, 0.4, 0.3])
                    ax.plot(x, self.fitfunction(x, *poptASIC), color=[0.3, 0.4, 0.3], label='Fit on average pixel, ASIC. \nE0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(poptASIC[-1], poptuncASIC[-1], poptASIC[0], poptuncASIC[0], poptASIC[2], poptuncASIC[2])) #Additions are to round up uncertainties. Hard coded
                ax.set_xlabel('Threshold (DAC)')
                ax.set_ylabel('Flux (hits/s)')
                ax.set_xlim(Thrs_ASIC[0], Thrs_ASIC[-1])
                ax.set_ylim(0, 1.75)
                ax.legend(fontsize=8, loc='upper right', bbox_to_anchor=(1,1))
        
        if residuals and pulls is not None:
            if ASIC_information != None:
                ax[0].errorbar(Thrs_ASIC, fluxASIC, yerr=unc_fluxASIC, fmt='o', color='blue', label='Flux on average pixel, ASIC, {}'.format(label))
                ax[0].plot(Thrs_ASIC, self.fitfunction(Thrs_ASIC, *poptASIC), color='blue', label='Fit on average pixel ASIC. \nE0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(poptASIC[-1], poptuncASIC[-1], poptASIC[0], poptuncASIC[0], poptASIC[2], poptuncASIC[2]))

                # Plot residuals on ax2
                residuals_ASIC = fluxASIC - self.fitfunction(np.array(Thrs_ASIC), *poptASIC)
                ax[1].plot(Thrs_ASIC, residuals_ASIC, marker='o', linestyle='-', color='blue', label = 'ASIC')
                ax[1].axhline(y=0, color='black', linestyle='--')
                
                # Plot residuals on ax3
                if unc_fluxASIC is not None:
                    pulls = [residuals_ASIC[entry]/unc_fluxASIC[entry] for entry in range(len(residuals_ASIC))]
                    ax[2].plot(Thrs_ASIC, pulls, marker='o', linestyle='-', color='blue', label = 'ASIC')
                    ax[2].axhline(y=0, color='black', linestyle='--')
                else:
                    print('unc_fluxASIC is not array:', unc_fluxASIC)
            
            if ASICgood_information is not None:
                ax[0].errorbar(Thrs_ASICgood, fluxASICgood, yerr=unc_fluxASICgood, fmt='o', label='Flux on average pixel, ASICgood, {}'.format(label))
                ax[0].plot(Thrs_ASICgood, self.fitfunction(Thrs_ASICgood, *poptASICgood),color = 'blue', label='Fit on average pixel, ASICgood. \nE0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(poptASICgood[-1], poptuncASICgood[-1], poptASICgood[0], poptuncASICgood[0], poptASICgood[2], poptuncASICgood[2]))
                
                residuals_ASICgood = fluxASICgood - self.fitfunction(np.array(Thrs_ASICgood), *poptASICgood)
                ax[1].plot(Thrs_ASICgood, residuals_ASICgood, marker='o', linestyle='-', label = 'ASICgood')

                pulls_ASICgood = [residuals_ASICgood[entry]/unc_fluxASICgood[entry] for entry in range(len(residuals_ASICgood))]
                ax[2].plot(Thrs_ASICgood, pulls_ASICgood, marker='o', linestyle='-', label = 'ASICgood')
            
        

        elif residuals is not None and pulls is None:
            if ASIC_information != None:
                ax[0].errorbar(Thrs_ASIC, fluxASIC, yerr=unc_fluxASIC, fmt='o', color='blue', label='Flux on average pixel, ASIC, {}'.format(label))
                ax[0].plot(Thrs_ASIC, self.fitfunction(Thrs_ASIC, *poptASIC), color='blue', label='Fit on average pixel ASIC. \nE0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(poptASIC[-1], poptuncASIC[-1], poptASIC[0], poptuncASIC[0], poptASIC[2], poptuncASIC[2]))

                # Plot residuals on ax2
                residuals_ASIC = fluxASIC - self.fitfunction(np.array(Thrs_ASIC), *poptASIC)
                ax[1].plot(Thrs_ASIC, residuals_ASIC, marker='o', linestyle='-', color='blue', label = 'ASIC')
                ax[1].axhline(y=0, color='black', linestyle='--')
            
            if ASICgood_information is not None:
                ax[0].errorbar(Thrs_ASICgood, fluxASICgood, yerr=unc_fluxASICgood, fmt='o', label='Flux on average pixel, ASICgood, {}'.format(label))
                ax[0].plot(Thrs_ASICgood, self.fitfunction(Thrs_ASICgood, *poptASICgood),color = 'blue', label='Fit on average pixel, ASICgood. \nE0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(poptASICgood[-1], poptuncASICgood[-1], poptASICgood[0], poptuncASICgood[0], poptASICgood[2], poptuncASICgood[2]))
                
                residuals_ASICgood = fluxASICgood - self.fitfunction(np.array(Thrs_ASICgood), *poptASICgood)
                ax[1].plot(Thrs_ASICgood, residuals_ASICgood, marker='o', linestyle='-', label = 'ASICgood')


        elif residuals is None and pulls is not None:
            if ASIC_information != None:
                ax[0].errorbar(Thrs_ASIC, fluxASIC, yerr=unc_fluxASIC, fmt='o', color='blue', label='Flux on average pixel, ASIC, {}'.format(label))
                ax[0].plot(Thrs_ASIC, self.fitfunction(Thrs_ASIC, *poptASIC), color='blue', label='Fit on average pixel ASIC. \nE0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(poptASIC[-1], poptuncASIC[-1], poptASIC[0], poptuncASIC[0], poptASIC[2], poptuncASIC[2]))

                # Plot residuals on ax2
                residuals_ASIC = fluxASIC - self.fitfunction(np.array(Thrs_ASIC), *poptASIC)

                if unc_fluxASIC is not None:
                    pulls = [residuals_ASIC[entry]/unc_fluxASIC[entry] for entry in range(len(residuals_ASIC))]
                    ax[1].plot(Thrs_ASIC, pulls, marker='o', linestyle='-', color='blue', label = 'ASIC')
                    ax[1].axhline(y=0, color='black', linestyle='--')
                else:
                    print('unc_fluxASIC is not array:', unc_fluxASIC)
            
            if ASICgood_information is not None:
                ax[0].errorbar(Thrs_ASICgood, fluxASICgood, yerr=unc_fluxASICgood, fmt='o', label='Flux on average pixel, ASICgood, {}'.format(label))
                ax[0].plot(Thrs_ASICgood, self.fitfunction(Thrs_ASICgood, *poptASICgood),color = 'blue', label='Fit on average pixel, ASICgood. \nE0={:.1f}$\pm${:.1f}, f={:.5f}$\pm${:.5f}, s={:.2f}$\pm${:.2f}'.format(poptASICgood[-1], poptuncASICgood[-1], poptASICgood[0], poptuncASICgood[0], poptASICgood[2], poptuncASICgood[2]))
                
                residuals_ASICgood = fluxASICgood - self.fitfunction(np.array(Thrs_ASICgood), *poptASICgood)
                

                pulls_ASICgood = [residuals_ASICgood[entry]/unc_fluxASICgood[entry] for entry in range(len(residuals_ASICgood))]
                ax[1].plot(Thrs_ASICgood, pulls_ASICgood, marker='o', linestyle='-', label = 'ASICgood')
        

            
            # Print residuals and uncertainties
            # limit = 30 
            # print(len(Thrs_ASIC), len(residuals_ASIC), len(unc_fluxASIC), len(pulls))
            # print(Thrs_ASIC[:limit])
            # print(residuals_ASIC[:limit])
            # print(unc_fluxASIC[:limit])
            # print(pulls[:limit])
            # Adjust layout to prevent overlap
            
        #Plotting asthetics
        ax[0].set_xlabel('Threshold (DAC)')
        ax[0].set_ylabel('Flux (hits/s)')
        ax[0].set_ylim(0,2)
        
        ax[0].legend(fontsize=8, loc='upper right', bbox_to_anchor=(1,1))

        if residuals is not None:
            ax[1].set_xlabel('Threshold (DAC)')
            ax[1].set_ylabel('Residuals')
            ax[1].set_title('Residuals of Flux')
            ax[1].legend(fontsize=8, loc='upper right', bbox_to_anchor=(1,1))
            ax[1].grid(True)
            ylimres = 0.3
            #ax[1].set_ylim(-ylimres,ylimres)
        

        if residuals and pulls is not None:
            ax[2].set_xlabel('Threshold (DAC)')
            ax[2].set_ylabel('Pulls')
            ax[2].set_title('Pulls of Flux')
            ax[2].legend(fontsize=8, loc='upper right', bbox_to_anchor=(1,1))
            ax[2].grid(True)
            ylimpulls = 10
            ax[2].set_ylim(-ylimpulls,ylimpulls)

        if residuals is None and pulls is not None:
            ax[1].set_xlabel('Threshold (DAC)')
            ax[1].set_ylabel('Pulls')
            ax[1].set_title('Pulls of ASIC Curve')
            ax[1].legend(fontsize=8, loc='upper right', bbox_to_anchor=(1,1))
            ax[1].grid(True)
            ylimpulls = 10
            ax[1].set_ylim(-ylimpulls,ylimpulls)

        plt.tight_layout()
            
    def trim_hist(self, ax, vp, opt = '', path = ''):
        data = self.open_csv(vp, nt.name['trim'], skiprows = 0, ravel = True, path = path)

        mask = self.build_mask_pattern(ravel = True, build4 = True)
        bins = np.arange(-0.5, 16.5, 1)
        if data is not None:
            ax.hist(data[mask == 2], bins = bins, histtype = 'step', linewidth = 1.3,
                color = 'k', label = 'even col.', alpha = 0.9)
            ax.hist(data[mask == 1], bins = bins, histtype = 'step', linewidth = 1.3,
                color = 'r', label = 'odd col.', alpha = 0.9)
            ax.hist(data[mask == 4], bins = bins, histtype = 'step', linewidth = 1.3,
                color = 'k', label = 'even col. 16th rows', alpha = 0.9, ls = 'dashed')   
            ax.hist(data[mask == 3], bins = bins, histtype = 'step', linewidth = 1.3,
                color = 'r', label = 'odd col. 16th rows', alpha = 0.9, ls = 'dashed')   

        if ax.has_data():
            ax.set_title(vp)
            ax.set_xlim(-0.5, 15.5)
            ax.set_ylim(0, 7000)
            ax.set_axisbelow(True)
            ax.grid()
            ax.legend()
            ax.set_xlabel('Trim', fontsize = 11)
            ax.set_ylabel('Number of counts', fontsize = 11)
            ax.tick_params(axis='both', which='major', labelsize=10)

    def n_corr(self, ax, vp, opt = 'mean', path = ''):
        if opt == 'mean':
            self.n_scatter(ax, vp, ['n_mean_t0', 'n_width_t0'], path = path)
            self.n_scatter(ax, vp, ['n_mean_t15', 'n_width_t15'], path = path)
            self.n_scatter(ax, vp, ['n_mean_t16', 'n_width_t16'], path = path)
        elif opt == 'rate':
            self.n_scatter(ax, vp, ['n_rate_t0', 'n_width_t0'], path = path)
            self.n_scatter(ax, vp, ['n_rate_t15', 'n_width_t15'], path = path)
            self.n_scatter(ax, vp, ['n_rate_t16', 'n_width_t16'], path = path)

        if ax.has_data():
            ax.set_title(vp)
            ax.set_xlim(1300, 1900)
            ax.set_ylim(0, 40)
            ax.grid()
            k_patch = mpatches.Patch(color = 'gray', label = 'All, even col.')
            r_patch = mpatches.Patch(color = 'r', label = 'All, odd col.')
            b_patch = mpatches.Patch(color = 'b', label = '16th, even col.')
            ax.legend(handles = [k_patch, r_patch, b_patch], loc = 'upper right', fontsize = 10)
            ax.set_xlabel('DAC code', fontsize = 11)
            ax.set_ylabel('Sigma', fontsize = 11)
            ax.tick_params(axis='both', which='major', labelsize=10) 
    
    def plot_target_edac_histograms( specified_result, results_collection_path = '/data/bfys/apuicerc/N037_new/results_collection.csv',specifics = [True, True, True]):
        """
        This function plots histograms for target and eDAC distributions based on the specified result and specifics.
        
        Args:
            results_collection_path (str): The path to the results_collection.csv file.
            specified_result (str): The specified result folder with the processed data to plot histograms for.
            specifics (list): A list of booleans indicating whether to plot specific values in the distributions. Contains 3 booleans for  mean, ASIC, and ASICgood.
        
        Returns:
            None
            Saves the plots in the specified N037_new/plots folder.
        """
        # Read the data from the results_collection.csv file into a DataFrame
        specfics_mean, specfics_ASIC, specfics_ASICgood = specifics

        print("Reading data from " + results_collection_path)
        data = pd.read_csv(results_collection_path)

        # Get all the parameters for the specified result
        result_params = data[data['Results'] == specified_result].drop('Results', axis=1)
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5), dpi=300)
        
        line_colors = ['red', 'blue']
        color_counter = 0

    
        for i in range(len(result_params)):
            
            # Get the parameters for the current result
            label = result_params.iloc[i][0]
            ST = result_params.iloc[i][1]
            baseline = result_params.iloc[i][2]
            result_param = result_params.iloc[i].tolist()[3:]
            
            targetmean, utargetmean, targetASIC, utargetASIC, targetASICgood, utargetASICgood, Kmean, uKmean, KASIC, uKASIC,KASICgood, uKASICgood = result_param
            # now result_param contains all the floats from results_collection.csv
            print('Plotting for dataset {} {}'.format(specified_result, label))
            
            target_matrix_path = os.path.join('/data/bfys/apuicerc/N037_new/', specified_result, 'target_dac_matrix_ST_{}_{}.csv'.format(ST, label))
            target_matrix = np.loadtxt(target_matrix_path,  dtype=float, delimiter=',')
            target_matrix = target_matrix[target_matrix>0]
            target_matrix = target_matrix[target_matrix<1000]

            edac_matrix_path = os.path.join('/data/bfys/apuicerc/N037_new/', specified_result, 'eDAC_matrix_ST_{}_{}.csv'.format(ST, label))
            edac_matrix = np.loadtxt(edac_matrix_path,  dtype=float, delimiter=',')

            
            if label == 'warm':
                color = 'orange'
            elif label == 'cold':
                color = '#add8e6'
            ax1.hist(target_matrix, bins = 50, alpha = 0.7, color = color, label = '{} {}'.format(label, 'target'))
            ax2.hist(edac_matrix, bins = 50, alpha = 0.7, color = color, label = '{} {}'.format(label, 'eDAC'))

    
            line_styles = ['-', '--', ':']
            line_idx = 0

            if specfics_mean and targetmean != 0:
                ax1.axvline(x = targetmean, color = line_colors[color_counter], linestyle = line_styles[line_idx], label = 'Mean = {:.2f} +- {:.2f} '.format(float(targetmean), float(utargetmean)))
                line_idx = (line_idx+1)%len(line_styles)
            if specfics_ASIC and targetASIC != 0:
                ax1.axvline(x = targetASIC, color = line_colors[color_counter], linestyle = line_styles[line_idx], label = 'ASIC = {:.2f} +- {:.2f}'.format(float(targetASIC), float(utargetASIC)))  ##np.mean(baseline_width[baseline_width>0])
                line_idx = (line_idx+1)%len(line_styles)
            if specfics_ASICgood and targetASICgood != 0 :
                ax1.axvline(x = targetASICgood, color = line_colors[color_counter], linestyle = line_styles[line_idx], label = 'ASICgood = {:.2f} +- {:.2f}'.format(float(targetASICgood), float(utargetASICgood)))
                line_idx = (line_idx+1)%len(line_styles)
            
            line_idx = 0
            if specfics_mean and Kmean != 0:
                ax2.axvline(x = Kmean, color = line_colors[color_counter], linestyle = line_styles[line_idx], label = 'Mean = {:.2f} +- {:.2f} '.format(float(Kmean), float(uKmean)))
                line_idx = (line_idx+1)%len(line_styles)
            if specfics_ASIC and KASIC != 0:
                ax2.axvline(x = KASIC, color = line_colors[color_counter], linestyle = line_styles[line_idx], label = 'ASIC = {:.2f} +- {:.2f}'.format(float(KASIC), float(uKASIC)))  ##np.mean(baseline_width[baseline_width>0])
                line_idx = (line_idx+1)%len(line_styles)
            if specfics_ASICgood and KASICgood != 0 :
                ax2.axvline(x = KASICgood, color = line_colors[color_counter], linestyle = line_styles[line_idx], label = 'ASICgood = {:.2f} +- {:.2f}'.format(float(KASICgood), float(uKASICgood)))
                line_idx = (line_idx+1)%len(line_styles)


            ax1.set_title('Target')
            ax2.set_title('eDAC')
            fig.suptitle('Target and eDAC distributions for {}\n'.format(specified_result))
            ax1.legend(bbox_to_anchor=(1.02, 1), loc='upper left')
            ax2.legend(bbox_to_anchor=(1.02, 1), loc='upper left')

            color_counter += 1    
    


        plt.tight_layout()
        plt.savefig('/data/bfys/apuicerc/N037_new/plots/' + specified_result + '_target_eDAC_histograms_warm_cold.png')
        print('Figure saved at {}'.format('/data/bfys/apuicerc/N037_new/plots/' + specified_result + '_target_eDAC_histograms_warm_cold.png'))
    

            


if __name__ == "__main__":
    folderpath = '/data/bfys/apuicerc/N037_new/resultscheck/'
    filename = 'E0matrix'
    shutter_time = '2s621ms'  # in milliseconds
    from_threshold = 1480  # in ADU
    to_threshold = 1600  # in ADU
    step_threshold = 1  # in ADU
    label = ''
    p0 = [0.0009, 8,12, 1551] #new stuff  # initial guess of parameters
    nohits_cut = 0.8  # mask the pixel if it has too many zero-valued points
    alpha = 0.05  # level of significance

    results_check_path = '/data/bfys/apuicerc/N037_new/resultscheck/'

    plot = Plot()
    #plot.open_csv(key = filename,label=label, skiprows = 1, ravel = True,path = folderpath)
