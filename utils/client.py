import argparse

parser = argparse.ArgumentParser()
subparser = parser.add_subparsers(help = 'Available options')

# a scheme for a parser for further development
dummy_parser = subparser.add_parser('...', help = '...')
dummy_parser.add_argument('dummy', action = None, choices = None, default = None,
                         nargs = None, metavar = None, help = None)


# analysis parser
analysis_parser = subparser.add_parser('analysis', help = 'handle analysis of files. Computes flux files, fits data and finds matrices with further parameters.')
analysis_parser.set_defaults(which = 'analysis')
analysis_parser.add_argument('generator', action = 'store', choices = ['flux', 'fitter' ,'matrices'],
                                help = 'select flux to find flux files, fitter to fit data and matrices to find target or e/DAC matrices')
analysis_parser.add_argument('--dataset', '-dt', action = 'store', const = 'check', default = 'check', nargs = '?',
                                help = 'select the dataset, default - check')
analysis_parser.add_argument('--temp', '-t', action = 'store', const = None , default = None, nargs = '?',
                                help = 'temperature, default - None. Analyses all temperatures within the dataset')
analysis_parser.add_argument('--s2_opt', '-s2', action = 'store', const = True , default = None, nargs = '?',
                                help = 'option for using fitfunction with 2 s parameters, default - None.')
analysis_parser.add_argument('--separation', '-sep', action = 'store', const = None , default = None,choices = ['even', 'odd', 'rows'], nargs = '?',
                                help = 'option for separation in even, odd columns 1and 16th rows, default - None.')



# noise analysis parser
noise_parser = subparser.add_parser('noise', help = 'noise analysis of a quick/normal/control scan')
noise_parser.set_defaults(which = 'noise')
noise_parser.add_argument('analysis', action = 'store', choices = ['groups'], 
                          help = 'analyse & visualise odd/even column and 16th rows pixel noise')
noise_parser.add_argument('scan_type', action = 'store', choices = ['quick', 'normal', 'control'],
                          help = 'select the type of the scan to analyse')




# plot parser
plot_parser = subparser.add_parser('plot', help = 'plots of all kinds for calibration data')
plot_parser.set_defaults(which = 'plot')

plot_parser.add_argument('--dataset', '-dt', action = 'store', const = None, default = None, nargs = '?',
                                help = 'select the dataset, default - None')
plot_parser.add_argument('--temp', '-t', action = 'store', const = None , default = None, nargs = '?',
                                help = 'temperature, default - None. Analyses all temperatures within the dataset')

plot_parser.add_argument('-plot_heatmap', '-hm', action = 'store_const', const = True, 
                         help = 'make a 2D heatmap of parameter for pixels')

plot_parser.add_argument('-plot_results_collection', '-rc', action = 'store_const', const = True, 
                         help = 'make a plot showing target and K[e/DAC] histograms for warm and cold scans.')

plot_parser.add_argument('-distribution', '-d', action = 'store_const', const = True, default = None,
                         help = 'make a histogram distribution of parameter. Can be specified with --param option.')

plot_parser.add_argument('-flux_over_thresholdscan', '-flux_vs_thr', action = 'store_const', const = True, default = None,
                         help = 'make a flux [Hits/s] vs threshold scan [DAC] plot (for individual pixels and ASIC).\n Can calculate residuals, later requested while code is executed.')
                         
plot_parser.add_argument('-opt', '-o', action = 'store', const = True, default = None,choices = [None,'discrepancy', 'QQ', 'summary'],nargs = '?',
                         help = 'addition option to plot, default - None. Discrepancy compares in number of std to a specific value. Summary creates plot showing results for K[e-/DAC].')                        

plot_parser.add_argument('--param', '-p', action = 'store', const = None , default = None,choices = ['E0', 'f', 's','target','Fittype', 'totalflux', 'K'], nargs = '?',
                                help = 'select parameter to plot, default - None. Analyses all parameters within the dataset.\n Can be used with: -plot_heatmap, -distribution')
                    
plot_parser.add_argument('--all', '-a', action = 'store_const', const = True,
                         help = 'make all available plots for certain type of scan. this doesnt work taken from other code')
plot_parser.add_argument('--separation', '-sep', action = 'store', const = None , default = None,choices = ['even', 'odd', 'rows'], nargs = '?',
                                help = 'option for separation in even, odd columns and 16th rows, default - None.')


#recipe parser
recipe_parser = subparser.add_parser('recipe', help = 'first recipe')
recipe_parser.set_defaults(which = 'recipe')
recipe_parser.add_argument('recipe_of', action = 'store', choices = ['threshold'], help =  '---')
recipe_parser.add_argument('--dir', '-d', action = 'store', required = False, default = '',
                         help = 'path to the equalisation directory')
