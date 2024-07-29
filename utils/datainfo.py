import numpy as np
import scipy
from scipy.optimize import curve_fit
from scipy.special import erfc


### FILE CONTAINING ALL DATA SPECIFICS



class FitFunctions:

    def fitfunction(self,E, f, A, s, E0):
        term0 = f*A /2 * s/np.sqrt(np.pi)*np.exp(-((E-E0)/s)**2)
        term1 = f * A / 2 * (E0 - E) * erfc((E - E0) / s)
        term2 = (1 - f) * A / (2 * np.sqrt(2 * np.pi) * s) * erfc((E - E0) / (np.sqrt(2) * s))
        
        return (term0 + term1+ term2)

    def fitfunctionnoterm0(self,E, f, A, s, E0):
        #term0 = f*A /2 * s/np.sqrt(np.pi)*np.exp(-((E-E0)/s)**2)
        term1 = f * A / 2 * (E0 - E) * erfc((E - E0) / s)
        term2 = (1 - f) * A / (2 * np.sqrt(2 * np.pi) * s) * erfc((E - E0) / (np.sqrt(2) * s))
        
        return (term1+ term2)
   
    
    def fitfunction2s(self,E, f, A, s,s2, E0):
        term0 = f*A /2 * s/np.sqrt(np.pi)*np.exp(-((E-E0)/s)**2)
        term1 = f * A / 2 * (E0 - E) * erfc((E - E0) / s)
        #term2 = (1 - f) * A / (2 * np.sqrt(2 * np.pi) * s) * erfc((E - E0) / (np.sqrt(2) * s))
        term2 = (1 - f) * A / (2 * np.sqrt(2 * np.pi) * s2) * erfc((E - E0) / (np.sqrt(2) * s2))
        
        return (term0+ term1 + term2) 

    def fitfunction2sAB(self,E, B, A, s,s2, E0):
        #tterm0 = A* s/np.sqrt(np.pi)*np.exp(-((E-E0)/s)**2)
        term1 = A* (E0 - E) * erfc((E - E0) / s)
        #term2 = (1 - f) * A / (2 * np.sqrt(2 * np.pi) * s) * erfc((E - E0) / (np.sqrt(2) * s))
        term2 = B/( np.sqrt(2 * np.pi) * s2) * erfc((E - E0) / (np.sqrt(2) * s2))
        
        return (term1 + term2) 

    def fitfunctionAB(self,E, B, A, s, E0):
        term0 = A* s/np.sqrt(np.pi)*np.exp(-((E-E0)/s)**2)
        term1 = A* (E0 - E) * erfc((E - E0) / s)
        term2 = B / ( np.sqrt(2 * np.pi) * s) * erfc((E - E0) / (np.sqrt(2) * s))
        
        return (term0 + term1 + term2)
    
    def fitfunctionABnoterm0(self,E, B, A, s, E0):
        #term0 = A* s/np.sqrt(np.pi)*np.exp(-((E-E0)/s)**2)
        term1 = A* (E0 - E) * erfc((E - E0) / s)
        term2 = B / ( np.sqrt(2 * np.pi) * s) * erfc((E - E0) / (np.sqrt(2) * s))
        
        return (term1 + term2)
    
    def fitfunction4(self,E, B, A, s, E0):
        #FITFUNCTION4WITHTERM0
        #term0 = A* s*np.sqrt(2/np.pi)*np.exp(-(1/2)*((E-E0)/s)**2)
        #term1 = A* (E0 - E) * erfc((E - E0) /(s*np.sqrt(2)))
        #term2 = B /2 * erfc((E - E0) / (np.sqrt(2) * s))

        term0 = A* s*np.sqrt(2/np.pi)*np.exp(-(1/2)*((E-E0)/s)**2)
        term1 = A* (E0 - E) * erfc((E - E0) /(s*np.sqrt(2)))
        term2 = B *np.sqrt(np.pi/2)*s* erfc((E - E0) / (np.sqrt(2) * s))
        return (term0 + term1 + term2) 

    def fitfunctionpaper(self,E, B, A, s, E0):
        
        term0 = A* np.exp(-(1/2)*((E-E0)/s)**2)
        term1 = B* erfc((E - E0) / (s*np.sqrt(2)))
        #term2 = B / ( np.sqrt(2 * np.pi) * s) * erfc((E - E0) / (np.sqrt(2) * s))
        
        return (term0 + term1)

    

class DataInfo:
    def __init__(self, type_data=None, label_to_consider=None, separation = None):
        #self.fit_functions = FitFunctions()
        info_list = self.information(type_data)
        self.assemble_new(info_list, label_to_consider,separation)

    
    def __dictflux__(self):
        """
        Returns a dictionary containing data specifics for fluxperpixel analysis.
        """
        dictflux = {}
        for i in range(len(self.datafolder_l)):
            dictflux[i] = {
                'datafolder': self.datafolder_l[i],
                'FPTfolder': self.FPTfolder_l[i],
                'ST': self.ST_l[i],
                'nacq': self.nacq_l[i],
                'fromThr': self.fromThr_l[i],
                'toThr': self.toThr_l[i],
                'stepThr': self.stepThr_l[i],
                'label': self.label_l[i],
                'bad_cut': self.bad_cut_l[i],
                'condition': self.condition_l[i]
            }
        return dictflux
    
    def information(self,type_data):

        # example of addition of a dataset to be analysed
        if type_data == 'new_dataset':
            datafolder_l = ['/path/to/data']*2
            FPTfolder_l = ['/path/to/results/folder']*2
            FPTfilename_l = ['FPT_filename']*2
            ST_l = ['ST_value']*2
            nacq_l = [100]*2
            fromThr_l = [1480]*2
            toThr_l = [1600]*2
            stepThr_l = [1]*2
            label_l = ['warm']*2
            bad_cut_l = [18.34765625]*2
            condition_l = ['False']*2
            p0 = [[0.0009, 8, 12, 1551]]*2
            p0_baseline = [[0.0009, 8, 12, 125]]*2
            baseline_opt = False
            separation_opt = False
            predict_opt_l = [False]*2
            output_opt_l = [False]*2


        if type_data == 'baseline':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultsultimate_baseline/'] #resultsultimate_baseline
            FPTfilename_l=2*['Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = True
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = self.fit_functions.fitfunction


        if type_data == 'newbaseline': #remeber to change nan_cut in pixelscanplot_basleine funtion for this dataset
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultsultimate_newbaseline/'] 
            FPTfilename_l=2*['Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 42, 44] 
            toThr_l=[ 193, 198] 
            stepThr_l=[1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction

        
        if type_data == 'acqbaseline':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultsultimate_acqbaseline/'] 
            FPTfilename_l=2*['Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 42, 44]
            toThr_l=[ 193, 198] 
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False'] 
            p0=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False  #set true if the data to be analysed gas alreasdy gone through fitter_baseline() and thus E0 is target
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction
        
        if type_data == 'baseline_information':

            ### IMPORTANT INFORMATION
            ### HERE NOISE BASELINE MATRIX FOR WARM AND COLD ARE GIVEN RESPECTIVELY, THIS MIGHT AFFECT SOME OF THE ANALYSIS CODE
            ### CHECK baseline_importer() in organiser to make sure the correct baseline is being used when temperature is specified

            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230626_warm/', '/data/bfys/apuicerc/N037_new/N037/20230622_cold/']
            ## importing the new info for the noise baseline, dont know if it is the same as above
            #baseline information required to import and process the files.
            FPTfolder_l=[ '/data/bfys/apuicerc/N037_new/N037/20230626_warm/', '/data/bfys/apuicerc/N037_new/N037/20230622_cold/']
            #FPTfolder_l= 4*['/data/bfys/apuicerc/N037_new/results4baseline/']
            FPTfilename_l=2*['Module0_VP3-1_TrimBest_Noise_Predict']
            ST_l=['1s310ms' , '2s621ms']
            nacq_l=[ 20, 40]
            fromThr_l=[ 1370, 1370]
            toThr_l=[ 1600, 1445] 
            stepThr_l=[5, 5]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            predict_opt_l = 2*[True] #True for obtaining baseline from equalisation process
            baseline_opt = False
            output_opt_l = 2*[True]
            separation_opt = False
            fit_function = None

        

        if type_data == 'calibration':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibration/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction

        if type_data == 'calibrationnoterm0':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationnoterm0/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunctionnoterm0
        
        if type_data == 'calibration_baseline':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/calibration_baseline/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = True
            separation_opt = None
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction

        if type_data == 'calibrationfitfunc2':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationfitfunc2/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12,8, 1551], [0.00032, 14.5, 17,8, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
        
        if type_data == 'calibrationfunc2noterm0':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationfunc2noterm0/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12,8, 1551], [0.00032, 14.5, 17,8, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction2s

        if type_data == 'calibration_fitfunc3':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationfitfunc3/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            #p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            #p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
    
            p0=[[(1-0.0009)*8/2, (0.0009*8)/2,12, 1551], [(1-0.00032)*14.5/2, (0.00032*14.5)/2, 17, 1535]] #new stuff
            p0_baseline=[[(1-0.0009)*8/2, (0.0009*8)/2,12, 125], [(1-0.00032)*14.5/2, (0.00032*14.5)/2, 17, 115]] #new stuff
             
            
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunctionABnoterm0

        if type_data == 'calibration_fitfunc3withterm0':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationfitfunc3withterm0/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            #p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            #p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
    
            p0=[[(1-0.0009)*8/2, (0.0009*8)/2,12, 1551], [(1-0.00032)*14.5/2, (0.00032*14.5)/2, 17, 1535]] #new stuff
            p0_baseline=[[(1-0.0009)*8/2, (0.0009*8)/2,12, 125], [(1-0.00032)*14.5/2, (0.00032*14.5)/2, 17, 115]] #new stuff
             
            
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunctionAB

        if type_data == 'calibration_fitfunc4withterm0':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationfitfunc4withterm0/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            #p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            #p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
    
            p0=[[(1-0.0009)*8/2, (0.0009*8)/2,12, 1551], [(1-0.00032)*14.5/2, (0.00032*14.5)/2, 17, 1535]] #new stuff
            p0_baseline=[[(1-0.0009)*8/2, (0.0009*8)/2,12, 125], [(1-0.00032)*14.5/2, (0.00032*14.5)/2, 17, 115]] #new stuff
             
            
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction4

        

        if type_data == 'calibrationfunc2sAB':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationfunc2sAB/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[(1-0.0009)*8/2, (0.0009*8)/2,12,8, 1551], [(1-0.00032)*14.5/2, (0.00032*14.5)/2, 17,8, 1535]] #new stuff
            p0_baseline=[[(1-0.0009)*8/2, (0.0009*8)/2, 12,8, 125], [(1-0.00032)*14.5/2, (0.00032*14.5)/2, 17,8, 115]]

            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction2sAB

        if type_data == 'calibration_fitfunctimepix':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationfitfunctimepix/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            #p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            #p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
    
            p0=[[0.7,0.1,2.6, 1551], [0.7,0.1,2.6, 1535]] #new stuff
            p0_baseline=[[0.7,0.1,2.6, 125], [0.7,0.1,2.6, 115]] #new stuff
             
            
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunctionpaper


        if type_data == 'calibration_separation':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibration_separation/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = True
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction

        if type_data == 'calibration_separationspp':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibration_separationspp/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = True
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction

        if type_data == 'calibration_equalT':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationequalT/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 50, 50] #[ 100, 50] HERE I HAVE CHANGE THE NACQ FOR THE FIRST DATA SET TO MAKE THEM HAVE THE SAME EXPOSURE TIME
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction
        
        if type_data == 'calibration_nacq25':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibration_nacq25/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 25, 25] #[ 100, 50] HERE I HAVE CHANGE THE NACQ FOR THE FIRST DATA SET TO MAKE THEM HAVE THE SAME EXPOSURE TIME
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction

        if type_data == 'calibration_nacq25to50':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibration_nacq25to50/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 50, 50] #[ 100, 50] HERE I HAVE CHANGE THE NACQ FOR THE FIRST DATA SET TO MAKE THEM HAVE THE SAME EXPOSURE TIME
            #for this dataset in which i take into account acquistions from 25 to 50 the flux code has to be modified such that the for loop iterating over the acqs go from range (24,nacq)
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction

        if type_data == 'calibration_newbaseline': #remeber to change nan_cut in pixelscanplot_basleine funtion for this dataset
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibration_newbaseline/'] 
            FPTfilename_l=2*['Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 42, 44] 
            toThr_l=[ 193, 198]
            #fromThr_l=[ 1480, 1480]
            #toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            
            stepThr_l=[1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction
        
        if type_data == 'calibration_debug':
            datafolder_l =[ '/data/bfys/apuicerc/N037_new/N037/20230803_warmFe55/', '/data/bfys/apuicerc/N037_new/N037/20230823_coldFe55/']
            
            FPTfolder_l= 2*['/data/bfys/apuicerc/N037_new/resultscalibrationdebug/'] 
            FPTfilename_l=2*['Fluxperpixel_Module0_VP3-1_ECS_data_ST_']
            ST_l=['2s621ms' , '2s621ms']
            nacq_l=[ 100, 50]
            fromThr_l=[ 1480, 1480]
            toThr_l=[ 1600, 1600] #A: recheck, all these values in the lists are obtained from the file names of data.
            stepThr_l=[ 1, 1]
            label_l = [ 'warm', 'cold']
            bad_cut_l = [18.34765625, 17.6941176471]
            condition_l = [ '(thr == 1494 and acq == 18) or (thr == 1527 and acq == 78)', 'False']
            p0=[[0.0009, 8,12, 1551], [0.00032, 14.5, 17, 1535]] #new stuff
            p0_baseline=[[0.0009, 8, 12, 125], [0.00032, 14.5, 17, 115]]
            baseline_opt = False
            separation_opt = False
            predict_opt_l = 2*[False] #True for obtaining baseline from equalisation process
            output_opt_l = 2*[False]
            fit_function = FitFunctions().fitfunction2s

        

        return [type_data ,datafolder_l, FPTfolder_l, FPTfilename_l, ST_l, nacq_l, fromThr_l, toThr_l, stepThr_l, label_l,bad_cut_l, condition_l, p0, p0_baseline, baseline_opt, separation_opt, predict_opt_l, output_opt_l, fit_function]

    def assemble_new(self, info_list, label_to_consider=None, separation = None):
        type_data, datafolder_l, FPTfolder_l, FPTfilename_l, ST_l, nacq_l, fromThr_l, toThr_l, stepThr_l, label_l,bad_cut_l, condition_l, p0, p0_baseline, baseline_opt, separation_opt, predict_opt_l, output_opt, fit_function = info_list
        
        if label_to_consider is None:
            
            self.datafolder_l = datafolder_l
            newFPTfolder_l = []
            for index in range(len(FPTfolder_l)):
                if separation != None:
                    newFPTfolder_l.append(FPTfolder_l[index]+str(separation)+'/')
                else:
                    newFPTfolder_l.append(FPTfolder_l[index])

            self.FPTfolder_l = newFPTfolder_l
            self.FPTfilename_l = FPTfilename_l
            self.ST_l = ST_l
            self.nacq_l = nacq_l
            self.fromThr_l = fromThr_l
            self.toThr_l = toThr_l
            self.stepThr_l = stepThr_l
            self.label_l = label_l
            self.bad_cut_l = bad_cut_l
            self.condition_l = condition_l
            self.p0 = p0
            self.p0_baseline = p0_baseline
            self.baseline_opt = baseline_opt
            self.separation_opt = separation_opt
            self.predict_opt_l = predict_opt_l
            self.output_opt = output_opt
            self.fit_function = fit_function
        else:
            index = [i for i, label in enumerate(label_l) if label == label_to_consider][0]
            
            self.datafolder_l = datafolder_l[index]
            if separation != None:
                self.FPTfolder_l = FPTfolder_l[index]+str(separation)+'/'
            else:
                self.FPTfolder_l = FPTfolder_l[index]

            self.FPTfilename_l = FPTfilename_l[index]
            self.ST_l = ST_l[index]
            self.nacq_l = nacq_l[index]
            self.fromThr_l = fromThr_l[index]
            self.toThr_l = toThr_l[index]
            self.stepThr_l = stepThr_l[index]
            self.label_l = label_l[index]
            self.bad_cut_l = bad_cut_l[index]
            self.condition_l = condition_l[index]
            self.p0 = p0[index]
            self.p0_baseline = p0_baseline[index]
            self.baseline_opt = baseline_opt
            self.separation_opt = separation_opt
            self.predict_opt_l = predict_opt_l[index]
            self.output_opt = output_opt[index]
            self.fit_function = fit_function



    def get_info(self):
        return (self.type_data,self.datafolder_l,self.FPTfolder_l, self.FPTfilename_l, self.ST_l, self.nacq_l, self.fromThr_l,
                self.toThr_l, self.stepThr_l, self.label_l,self.bad_cut_l, self.condition_l, self.p0,
                self.p0_baseline, self.baseline_opt, self.separation_opt, self.predict_opt_l,
                self.output_opt, self.fit_function)
    




#example_data_info = DataInfo()
#example_data_info.create_dataset()
#print(example_data_info.get_info())
#print(example_data_info)

#print(example_data_info.get_info())
