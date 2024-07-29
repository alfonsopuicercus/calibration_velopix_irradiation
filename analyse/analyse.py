import time
import numpy as np
import glob as glob
import scipy
from scipy.optimize import curve_fit
from scipy.special import erfc
import scipy.stats as stats
from utils import organizer
import os
import inspect #used to find the number of parameters a function takes



class AnalyseFlux():
    def __init__(
            self,
            datafolder_l=None,
            FPTfolder_l=None,
            FPTfilename_l=None,
            ST_l=None,
            nacq_l=None,
            fromThr_l=None,
            toThr_l=None,
            stepThr_l=None,
            label_l=None,
            bad_cut_l=None,
            condition_l=None,
            p0=None,
            p0_baseline=None,
            baseline_opt=None,
            separation_opt=None,
            predict_opt_l=None,
            output_opt=None,
            fit_function = None
    ):
        """
        """

        # Checks if there is a separation in groups of pixels which shoudl mena there are subfolders for the groups
        separation_category =  str(FPTfolder_l.split('/')[-2]) #shoudl be even odd or rows if group separation exisitent
        # Mapping categories to their corresponding values
        category_map = {
            'even': 2,
            'odd': 3,
            'rows': 4
        }

        # Get the corresponding cat value or None if not found
        cat = category_map.get(separation_category, None)
        

        if FPTfolder_l.split('/')[-2] == 'resultscalibration_newbaseline':
            print('using fluxpixelthr_ultimate_baseline for newbaseline approach')
            self.fluxpixelthr_ultimate_baseline(
                        filepath=datafolder_l,
                        savepath=FPTfolder_l,
                        ST=ST_l,
                        nacq=nacq_l,
                        fromThr=fromThr_l,
                        toThr=toThr_l,
                        stepThr=stepThr_l,
                        label=label_l,
                        bad_cut=bad_cut_l,
                        condition=condition_l)


        else:
            
            self.fluxpixelthr_ultimate(
                filepath=datafolder_l,
                savepath=FPTfolder_l,
                ST=ST_l,
                nacq=nacq_l,
                fromThr=fromThr_l,
                toThr=toThr_l,
                stepThr=stepThr_l,
                label=label_l,
                bad_cut=bad_cut_l,
                condition=condition_l, separation_opt = separation_opt, cat = cat)
        

    def fluxpixelthr_ultimate(self,filepath,savepath, ST, nacq, fromThr, toThr, stepThr, label, bad_cut, condition,separation_opt=None, cat = None,CUT_BAD=True, redo = False, good_only = False, ignore_firstrow = False):
        '''
        This one is used for the ultimate dataset. based function for all the others

        Creates a list of flux per pixel for each threshold, the uncertainty of the flux per pixel for each threshold,
        the total flux on the ASIC per threshold and the uncertainty of the total flux on the ASIC per threshold.
        The analyzed data files are named in the form <filepath+fileprefix+str(ST)+'_THR_'+str(thr)+'_'+str(acq)+'.csv'>
        (e.g. /data/bfys/LucasvdH/VELOdata/coolsource/Module0_VP0-1_ECS_data_ST_1s310ms_THR_1320_3.csv).
        The saved files will be named in the form <filepath+'<datatype>'+fileprefix+str(ST)+'_THR_'+str(thr)+'.csv'>, 
        where datatype is Fluxperpixel, FluxperASIC, UncertaintyFluxperpixel or UncertaintyFluxperASIC.

        Arguments:
        - filepath (str): Filepath (folder) to the data (e.g. '/data/bfys/LucasvdH/VELOdata/coolsource/').
        - fileprefix (str): Prefix of the data files (e.g. 'Module0_VP0-1_ECS_data_ST_').
        - ST (str): shutter time (e.g. '1s310ms').
        - nacq (int): Number of acquisitions per threshold (e.g. 100).
        - fromThr (int): Starting threshold of the scan (e.g. 1500).
        - toThr (int): Ending threshold of the scan (e.g. 1600).
        - stepThr (int): Increment of the threshold between scans (e.g. 5).
        - label (str): Specifies cooling condition for better file recognition and post-management (e.g. 'cool').
        - CUT_BAD (bool): Default True. Skip data acquisitions that contain bad data. 
        - bad_cut (float, int): Default 2. Used to identify bad data acquisitions. 
        If flux per pixel on the last columns exceeds bad_cut, it is faulty. 
        '''
        Thrs=[]
        Flux=[]
        F_ASIC=np.zeros((int((toThr-fromThr)/stepThr)+1, 1))
        unc_F_ASIC=np.zeros((int((toThr-fromThr)/stepThr)+1, 1))
        fileprefix = 'Module0_VP3-1_ECS_data_ST_'
        print('separation_opt', separation_opt)
        try:
            mask=AnalyseFlux.get_mask(filepath, fileprefix) # Uses function to obtain the desired mask, in this case mask matrix from the equalisation process
            MASK = True
            print('MASK', MASK)

        except:
            mask = np.zeros((256,256))
            MASK = False
        
        # Check for separation option, if not None then create the mask for separation into groups of pixels.
        if separation_opt != None and separation_opt != False:
            MASK_separation = True
            print('MASK_separation',MASK_separation)
            if cat == 4:
                mask_separation = self.create_mask_pattern(build4=True) 
            else:
                mask_separation = self.create_mask_pattern(build4=False) # Uses function to create the desired mask
        else:
            MASK_separation = None
            mask_separation = np.zeros((256,256))

        # Obtains the acquisition time in seconds
        st=ST.split('ms')
        st=st[0].split('s')
        sec=float(st[0])
        milisec=float(st[1])/1000
        acqtime=sec+milisec

        
        print('\n\nANALYZING for {}{}{}.csv'.format(filepath, fileprefix, str(ST))) # Iterates over all thresholds
        for thr in range(fromThr, toThr+1, stepThr):
            Thrs.append(thr)
            if thr % 5 == 0:
                print('Analyzing for thr={}'.format(thr)) # Prints the current threshold every 5th threshold
            singlehits = []
            #t_hits: Intermediary variable to sum the total hits, thereafter converted to flux
        
            t_hits=np.zeros((256, 256)) #A: here creating a grid to store the number of hit per pixel
            goodacq=0
            emptyscan = 0
            badscan = 0
            bad_cut_scan = 0

            for acq in range(nacq): # Iterates over all acquisitions
                if eval(condition):  # Checks for conditions to skip acquisition (due to bad data takign or issues in the data)
                    continue

                if ignore_firstrow == 'True': # if True, ignore the first row of the data
                    hits=np.loadtxt(filepath+fileprefix+str(ST)+'_THR_'+str(thr)+'_'+str(acq)+'.csv', dtype=int, delimiter=',', skip_header = 1)
                else:
                    hits=np.loadtxt(filepath+fileprefix+str(ST)+'_THR_'+str(thr)+'_'+str(acq)+'.csv', dtype=int, delimiter=',')


                if (np.isclose(np.mean(hits), 0.0, rtol=1e-09, atol=1e-09)): #Check whether the scan is empty and if so it continues to the next acquisition.
                    emptyscan+=1 #Add 1 to the number of empty scans
                    continue

            
                if (MASK==True): # if True, applies the mask
                    hits = np.ma.masked_where(mask>0, hits)  # Masks pixels based on mask matrix from equalisation

                    copy_hits_mask = hits.copy() # Copies the hits matrix to not disturb the original data for bad scan check
      
                    if MASK_separation == True: # if True, applies the mask for separation of pixels in groups
                        
                        hits_scan =np.ma.masked_where(mask_separation != cat, hits) # Masks pixels based on mask matrix from equalisation
                      
                        scan = np.sum(hits_scan.compressed()) # Counts the number of hits
             
                        hits = np.ma.filled(hits_scan, fill_value=0)
                    else:
                        scan = np.sum(hits.compressed()) # Counts the number of hits
                else:
                    scan = np.sum(hits)

                mean_lastcol, BAD=AnalyseFlux.check_badscan(copy_hits_mask, bad_cut)  # Checks whether the scan is bad
                #print(mean_scan, BAD) 
                
                if (BAD == True): # If the scan is bad
                    badscan +=1 #Add 1 to the number of bad scans
                    if (CUT_BAD == True): # If True, cut the scan
                        bad_cut_scan += 1 #Add 1 to the number of bad cut scans
                        continue

                singlehits.append(scan) # Adds the number of hits in this scan to the list
                t_hits+=hits # Adds the number of hits in this scan to the hits grid
                goodacq+=1 #if data acquisition no defective, considered good for analysis
            
            debugging_scans = False # Set to True to print debugging information on the scans
            if debugging_scans == True:
                print(np.mean(hits),goodacq/nacq, emptyscan/nacq, badscan/nacq, bad_cut_scan/nacq)
    

            totalacqtime=goodacq*acqtime #Total acquisition time in seconds

            #Saves the flux on each pixel and, seperately, its uncertainty, as a 256x256 array
            np.savetxt(savepath+'Fluxperpixel_'+fileprefix+str(ST)+'_'+str(label)+'_THR_'+str(thr)+'.csv', t_hits/totalacqtime, delimiter=',')
            np.savetxt(savepath+'UncertaintyFluxperpixel_'+fileprefix+str(ST)+'_'+str(label)+'_THR_'+str(thr)+'.csv', np.sqrt(t_hits)/totalacqtime, delimiter=',')


            mean = sum(singlehits)/float(goodacq)
            squared_diff_sum = np.sum((x - mean)**2 for x in singlehits)

            
            F_ASIC[int((thr-fromThr)/stepThr)][0]=mean/(acqtime*256*256) #Calculates ASIC average flux
            unc_F_ASIC[int((thr-fromThr)/stepThr)][0]=np.sqrt(squared_diff_sum/goodacq)/(acqtime*256*256) #Calculates ASIC uncertainty
            
      

        #Saves the total flux on the ASIC per threshold and, separately, its uncertainty, as a numberofthrs x 1 array
        np.savetxt(savepath+'FluxperASIC_'+fileprefix+str(ST)+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', F_ASIC, delimiter=',')
        np.savetxt(savepath+'UncertaintyFluxperASIC_'+fileprefix+str(ST)+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', unc_F_ASIC, delimiter=',')

        print('All Flux files saved for {} {}.'.format(savepath, label))



    
    def fluxpixelthr_ultimate_baseline(self,filepath,savepath, ST, nacq, fromThr, toThr, stepThr, label, bad_cut, condition,separation_opt=None, cat = None,CUT_BAD=True, redo = False, good_only = False, ignore_firstrow = False):
        '''
        This one is used for the newbaseline dataset. Same as fluxpixelthr() but with an additional baseline correction before calculating the flux.


        Creates a list of flux per pixel for each threshold, the uncertainty of the flux per pixel for each threshold,
        the total flux on the ASIC per threshold and the uncertainty of the total flux on the ASIC per threshold.
        The analyzed data files are named in the form <filepath+fileprefix+str(ST)+'_THR_'+str(thr)+'_'+str(acq)+'.csv'>
        (e.g. /data/bfys/LucasvdH/VELOdata/coolsource/Module0_VP0-1_ECS_data_ST_1s310ms_THR_1320_3.csv).
        The saved files will be named in the form <filepath+'<datatype>'+fileprefix+str(ST)+'_THR_'+str(thr)+'.csv'>, 
        where datatype is Fluxperpixel, FluxperASIC, UncertaintyFluxperpixel or UncertaintyFluxperASIC.

        Arguments:
        - filepath (str): Filepath (folder) to the data (e.g. '/data/bfys/LucasvdH/VELOdata/coolsource/').
        - fileprefix (str): Prefix of the data files (e.g. 'Module0_VP0-1_ECS_data_ST_').
        - ST (str): shutter time (e.g. '1s310ms').
        - nacq (int): Number of acquisitions per threshold (e.g. 100).
        - fromThr (int): Starting threshold of the scan (e.g. 1500).
        - toThr (int): Ending threshold of the scan (e.g. 1600).
        - stepThr (int): Increment of the threshold between scans (e.g. 5).
        - label (str): Specifies cooling condition for better file recognition and post-management (e.g. 'cool').
        - CUT_BAD (bool): Default True. Skip data acquisitions that contain bad data. 
        - bad_cut (float, int): Default 2. Used to identify bad data acquisitions. 
        If flux per pixel on the last columns exceeds bad_cut, it is faulty. 
        '''
        
        Thrs=[]
        Flux=[]
        fileprefix = 'Module0_VP3-1_ECS_data_ST_'
        
        
        try:
            mask=AnalyseFlux.get_mask(filepath, fileprefix)
            MASK = True
            print('MASK', MASK)

        except:
            mask = np.zeros((256,256))
            MASK = False

        MASK_separation = None
        mask_separation = np.zeros((256,256))

        st=ST.split('ms')
        st=st[0].split('s')
        sec=float(st[0])
        milisec=float(st[1])/1000
        acqtime=sec+milisec
        
        baseline_matrix = Fitter().baseline_importer(temp = label)
        baseline_ASIC = np.mean(baseline_matrix)
        baseline_matrix_max = int(baseline_matrix[baseline_matrix>0].max())
        baseline_matrix_min = int(baseline_matrix[baseline_matrix>0].min())

        print('baseline_matrix_min, baseline_matrix_max, empty', baseline_matrix_min, baseline_matrix_max, 256*256-len(baseline_matrix[baseline_matrix>0]))
        

        newthr_max = toThr-baseline_matrix_min
        newthr_min = fromThr-baseline_matrix_max
        newthr_range = np.arange(newthr_min, newthr_max+1,1)
        
        #saving in the array below all all the goodacq in the thr range for later use of computation of ASIC error
        goodacq_l = []

        print('thr_min, thr_max, len(thr)', fromThr, toThr, len(np.arange(fromThr, toThr+1, stepThr)))
        print('newthr_min, newthr_max, len(newthr_range)', newthr_min, newthr_max, len(newthr_range))
        print('int((newthr_max-newthr_min)/stepThr+1)', int((newthr_max-newthr_min)/stepThr+1))
        array_of_hitmatrices = np.zeros((len(newthr_range), 256, 256))
        array_of_errhitmatrices = np.zeros((len(newthr_range), 256, 256))

        array_of_singlehits = np.zeros((len(newthr_range), nacq))

        F_ASIC=np.zeros((int((newthr_max-newthr_min)/stepThr+1), 1))
        unc_F_ASIC=np.zeros((int((newthr_max-newthr_min)/stepThr+1), 1))


        print('\n\nANALYZING for {}{}{}.csv'.format(filepath, fileprefix, str(ST)))
        for thr in range(fromThr, toThr+1, stepThr):
            dummy_thr = 0
            Thrs.append(thr)
            print('Analyzing for thr={}'.format(thr))
            singlehits = []
            #t_hits: Intermediary variable to sum the total hits, thereafter converted to flux
            #t_hits=np.ma.masked_where(mask>0, np.zeros((256, 256), int))
            t_hits=np.zeros((256, 256)) #A: here creating a grid to store the number of hit per pixel
            goodacq=0
            emptyscan_MD = 0
            emptyscan_old = 0
            badscan = 0
            bad_cut_scan = 0

            for acq in range(nacq):
                if eval(condition):
                    continue

                if ignore_firstrow == 'True':
                    hits=np.loadtxt(filepath+fileprefix+str(ST)+'_THR_'+str(thr)+'_'+str(acq)+'.csv', dtype=int, delimiter=',', skip_header = 1)
                else:
                    hits=np.loadtxt(filepath+fileprefix+str(ST)+'_THR_'+str(thr)+'_'+str(acq)+'.csv', dtype=int, delimiter=',')


                if (np.isclose(np.mean(hits), 0.0, rtol=1e-09, atol=1e-09)):
                    emptyscan_MD+=1
                    continue

                if (np.mean(hits)==0.0): #A: recheck, this is added as a suggestion from MD, i think it check whether the scan is empty and if so it continues to the next acquisition.
                    emptyscan_old += 1
                    continue

                if (MASK==True):
                    hits = np.ma.masked_where(mask>0, hits)

                    copy_hits_mask = hits.copy()
      
                    if MASK_separation == True:
                        
                        hits_scan =np.ma.masked_where(mask_separation != cat, hits)
                      
                        scan = np.sum(hits_scan.compressed())
             
                        hits = np.ma.filled(hits_scan, fill_value=0)
                    else:
                        scan = np.sum(hits.compressed())
                else:
                    scan = np.sum(hits)

                mean_lastcol, BAD=AnalyseFlux.check_badscan(copy_hits_mask, bad_cut)
                #print(mean_scan, BAD)
                
                if (BAD == True):
                    badscan +=1
                    if (CUT_BAD == True):
                        bad_cut_scan += 1
                        continue

                singlehits.append(scan)
                array_of_singlehits[dummy_thr][acq] = scan
                t_hits+=hits
                goodacq+=1 #if data acquisition no defective, considered good for analysis
            print(nacq)
            print('mean(hits),goodacq/nacq, emptyscan_MD/nacq,emptyscan_old/nacq, badscan/nacq, bad_cut_scan/nacq')
            print(np.mean(hits),goodacq/nacq, emptyscan_MD/nacq,emptyscan_old/nacq, badscan/nacq, bad_cut_scan/nacq)
            
            goodacq_l.append(goodacq)
            
            totalacqtime=goodacq*acqtime
            for row in range(0,256):
                for column in range(0,256):

                    if t_hits[row][column] <= 0:
                        continue
                        
                    if baseline_matrix[row][column] <= 0:
                        continue

                    new_thr = int(thr - baseline_matrix[row][column])
                    for i in range(len(newthr_range)):
                        if new_thr == newthr_range[i]:
                            index = i

                    array_of_hitmatrices[index][row][column] = t_hits[row][column]/goodacq
                    array_of_errhitmatrices[index][row][column] = np.sqrt(t_hits[row][column])/goodacq
    
            
        for i in range(len(newthr_range)):
            thr_shift = newthr_range[i]
            #Saves the flux on each pixel and, seperately, its uncertainty, as a 256x256 array
            np.savetxt(savepath+'Fluxperpixel_'+fileprefix+str(ST)+'_'+str(label)+'_THR_'+str(int(thr_shift))+'.csv', array_of_hitmatrices[i]/acqtime, delimiter=',')
            np.savetxt(savepath+'UncertaintyFluxperpixel_'+fileprefix+str(ST)+'_'+str(label)+'_THR_'+str(int(thr_shift))+'.csv', array_of_errhitmatrices[i]/acqtime, delimiter=',')
                        

            mean = np.sum(array_of_hitmatrices[i])
            squared_diff_sum = 0

            print('int((thr_shift-newthr_min)/stepThr), mean/(acqtime*256*256)')
            print(int((thr_shift-newthr_min)/stepThr), mean/(acqtime*256*256))
            #Saves the total flux on the ASIC per threshold and, separately, its uncertainty, as a numberofthrs x 1 array
            F_ASIC[int((thr_shift-newthr_min)/stepThr)][0]=mean/(acqtime*256*256)  #A: recheck nansum chnages the nan values for zero so that it can perform the sum
            unc_F_ASIC[int((thr_shift-newthr_min)/stepThr)][0]=np.sqrt(squared_diff_sum/goodacq)/(acqtime*256*256)
            #print('t_hits/totalacqtime')
            #print(t_hits/totalacqtime,)
            print('FASIC mean/(acqtime*256*256)')
            print( mean/(acqtime*256*256))

        np.savetxt(savepath+'FluxperASIC_'+fileprefix+str(ST)+'_'+str(label)+'_THR_'+str(int(newthr_min))+'-'+str(int(newthr_max))+'-step'+str(stepThr)+'.csv', F_ASIC, delimiter=',')
        np.savetxt(savepath+'UncertaintyFluxperASIC_'+fileprefix+str(ST)+'_'+str(label)+'_THR_'+str(int(newthr_min))+'-'+str(int(newthr_max))+'-step'+str(stepThr)+'.csv', unc_F_ASIC, delimiter=',')


    def create_mask_pattern(self,ravel = False, build4 = False, size = 256, separation_type = None):
        '''
        Creates mask pattern to categorise pixels in even and odd columns as well as every 16th row.
        This was done due to a biased patter observed in the equalisation while following this type of categorisation.
        Sets:
            - even = 2
            - odd = 3
            - 16th rows = 4

        Args:
            ravel (bool, optional): If True, transforms the matrix into an flat/1D array. Defaults to False.
            build4 (bool, optional): If True, includes 16th rows categorisation. Defaults to False.
            size (int, optional): Sets size of desired square matrix. Defaults to 256.

        Returns:
            matrix (array): mask matrix with the categorisation
        '''
        
        matrix = np.zeros((size,size))  # Initialize a 256x256 matrix with zeros
        separation_type = 'evenodd16'
        #separation_type = 'superpixelmask'
    
        print('Using pixel separation type:', separation_type) 

        if separation_type == 'evenodd16' or separation_type == None: #Separates the pixels in even and odd columns
            for row in range(size):
                for col in range(size):
                    if col % 2 == 0:
                        matrix[row][col] = 2  # Assign 2 to even columns
                    else:
                        matrix[row][col] = 3  # Assign 3 to odd columns
                    
                    if row % 15 == 0 and row != 0 and build4 == True:
                        matrix[row][col] = 4  # Assign 4 to the 16th row
        
        elif separation_type == 'superpixelmask': # Separates the pixels in super pixel structure (pairs of columns)
            for row in range(size):
                for col in range(size):
                    if (col // 2) % 2 == 0:
                        matrix[row][col] = 2  # Assign 2 to pairs of columns (0,1,4,5,8,9,...)
                    else:
                        matrix[row][col] = 3  # Assign 3 to pairs of columns (2,3,6,7,10,11,...)
                    
                    if row % 15 == 0 and row != 0 and build4 == True:
                        matrix[row][col] = 4  # Assign 4 to the 16th row
        
        if ravel == True: # Transforms the matrix into an flat/1D array
            return matrix.ravel()
        else:
            return matrix # Returns the matrix



    def check_badscan(matrix, bad_cut):
        '''Checks for badscan based on certain conditions. Used in fluxperthr generation of files.
        The check is performed ont he right-most column since during some scans it is possible to see
        biases in the measurements where this columns appears as all pixels firing up.

        Args:
            matrix (array): Hit matrix from data.
            bad_cut (float): Condition that if met, scan considered as bad scan. Specified in the data specifics together with other parameters.

        Returns:
            mean (float): Mean hits of the right-most column
            BAD (bool): True or False. Sets the input scan to be considered as good or bad (the does not consider it for flux computation) respectively. 
        '''

        BAD = False
        #divide in columns
        """
        column1 = 0
        column2= 7
        means = []
        sigmas = []
        while (column2 < 256):
        """
        slice=matrix[:, -1:]
        #print('slice', np.size(slice))
        mean=slice.mean()
        """
        sigmas.append(slice.std())
        column1+=8
        column2+=8
        """
        #means = np.array(means,dtype='f')
        #sigmas = np.array(sigmas,dtype='f')
        if (mean>bad_cut-0.01 and mean< bad_cut+0.01):
            BAD = True
        return mean, BAD # Returns the mean of the right-most column and the boolean BAD (see above)

    def get_mask(filepath, filename):
        '''
        Finds mask matrix within the data folder. Used in fluxpixelthr function.

        Args:
            filepath (string): path where the mask matrix is contained.
            filename (string): name of the files to be analysed. This contains information on Module and ASIC which is used to find the mask.

        Returns:
            mask (array): mask matrix containing the location of the pixels that need to be masked based on the result from the equalisation process.
        '''
        module = filename.split("_")[0]
        velopix = filename.split("_")[1]
        mask = np.loadtxt(filepath+module+"_"+velopix+"_Matrix_Mask.csv", dtype=int, delimiter = ',')
        return mask # Returns the mask



class Fitter():

    def __init__(self,
            FPTfolder_baseline_l = None,
            s2_opt = None,
            datafolder_l=None,
            FPTfolder_l=None,
            FPTfilename_l=None,
            ST_l=None,
            nacq_l=None,
            fromThr_l=None,
            toThr_l=None,
            stepThr_l=None,
            label_l=None,
            bad_cut_l=None,
            condition_l=None,
            p0=None,
            p0_baseline=None,
            baseline_opt=None,
            separation_opt=None,
            predict_opt_l=None,
            output_opt=None,
            fit_function = None):

        if any(arg is not None for arg in (FPTfolder_baseline_l, datafolder_l, FPTfolder_l, FPTfilename_l,
                                            ST_l, nacq_l, fromThr_l, toThr_l, stepThr_l, label_l,
                                            bad_cut_l, condition_l, p0, p0_baseline, baseline_opt,
                                            separation_opt, predict_opt_l, output_opt, fit_function)):
            # Call the function if any argument is not None

            self.fitfunction = fit_function #importing the desired fitfunction to use by the dataset
            
            if s2_opt is not None:
                print('Fitting using fitfunc with 2 s parameter')
                #self.fit_function = self.fitfunction2

                if FPTfolder_baseline_l is not None:
                    
                    self.fitter2(FPTfolder_l, FPTfolder_baseline_l, ST_l, fromThr_l, toThr_l, stepThr_l, label_l, p0_baseline)
                else:
                    self.fitter2(FPTfolder_l, FPTfolder_baseline_l, ST_l, fromThr_l, toThr_l, stepThr_l, label_l, p0)

            else:
                print('Fitting using fitfunc with 1 s parameter')

                if FPTfolder_l.split('/')[-2] == 'resultscalibrationfitfunc4withterm0':
                    print('Using fitter4 function...')
                    self.fitfunction = self.fitfunction4
                    if FPTfolder_baseline_l is not None:
                        self.fitter4(FPTfolder_l, FPTfolder_baseline_l, ST_l, fromThr_l, toThr_l, stepThr_l, label_l, p0_baseline)
                    else:
                        self.fitter4(FPTfolder_l, FPTfolder_baseline_l, ST_l, fromThr_l, toThr_l, stepThr_l, label_l, p0)


                else:
                    print('Using default fitter function...')
                
                    signature = inspect.signature(self.fitfunction)
                    params_l = list(signature.parameters)[1:] 
                
                    n_parameters = len(params_l)

                    print('Using {} fitter function with {} parameters {} ...'.format(self.fitfunction.__name__, n_parameters, params_l ))
                    if FPTfolder_baseline_l is not None:
                        self.fitter(FPTfolder_l, FPTfolder_baseline_l, ST_l, fromThr_l, toThr_l, stepThr_l, label_l, p0_baseline, params_l = params_l)
                    else:
                        self.fitter(FPTfolder_l, FPTfolder_baseline_l, ST_l, fromThr_l, toThr_l, stepThr_l, label_l, p0, params_l = params_l)

            
    def chisquare(self,O, E, Oerr):
        '''
        Returns Chi square. If the uncertainty is close to 0 (raising error) or X2 is NaN, it assigns 0 to X2.
    
        Arguments:
        - O (array): observed Flux. 
        - E (array): estimated Flux from fit.
        - Oerr (array): uncertainty of the observed Flux.
        '''
        O=O.flatten() # Flatten the array
        E=E.flatten()
        Oerr=Oerr.flatten()
        X2=0
        for i in range(len(O)):
            try:
                X2+=(O[i]-E[i])**2/Oerr[i]**2
            except ZeroDivisionError:
                X2=0
                continue
        if np.isnan(X2):
            X2=0
    
        return X2 
    
    

    def baseline_importer(self, FPTfolder_baseline = None ,width = False, ASIC_baseline =False, temp = None):
        '''
        Imports the baseline matrix previously found with the equalisation process.
    
        Args:
            baselineparams (array): array containing information to locate the baseline matrix file.
            ASIC_baseline (bool, optional): Set to True to return ASIC_baseline which is assumed to be the mean of all baseline matrix entries. Defaults to False.
    
        Returns:
            baseline (array): Matrix containing the noise baseline position in DAC units for every pixel in the ASIC.
        '''

        FPTfolder_baseline_l = ['/data/bfys/apuicerc/N037_new/N037/20230626_warm/', '/data/bfys/apuicerc/N037_new/N037/20230622_cold/']
        
        # selecting the folder where noise baseline matrix is stored. For the two datasets used for the analysis defaults apply, else, it takes the input from FPT_baseline
        if temp == 'warm':
            FPTfolder_baseline = FPTfolder_baseline_l[0]
        elif temp == 'cold':
            FPTfolder_baseline = FPTfolder_baseline_l[1]
        else:
            FPTfolder_baseline = FPTfolder_baseline # input from FPT_baseline


            
        #importing the files
        baseline=np.loadtxt(FPTfolder_baseline+'Module0_VP3-1_TrimBest_Noise_Predict'+'.csv', dtype=float, delimiter=',') #this is the noise baseline
        baseline_width = np.loadtxt(FPTfolder_baseline+'Module0_VP3-1_Trim0_Noise_Width'+'.csv', dtype=float, delimiter=',') #this is the noise baseline width

        baseline_mean = int(np.mean(baseline[baseline>0])) #this is the mean of the noise baseline
        
        if ASIC_baseline == True:
            return baseline_mean #outputs the mean of the noise baseline, which is considered and used for ASIC flux analysis
        if width == True:
            #print('Importing noise baseline width...')
            return baseline_width #outputs the width of the noise baseline
        else:
            #print('Importing noise baseline...')
            return baseline #outputs the noise baseline matrix
    
    

    def fitter(self,FPTfolder, FPTfolder_baseline, ST, fromThr, toThr, stepThr, label, p0,params_l, nohits_cut=0.8, alpha=0.05):
        '''
        This one is fitter. should be similar to that of Lucas but  most likely containing corrections or additions to the code.

        Finds best fit to the flux per pixel and total flux on pixel grid to fitfunction (defined below). Saves eight 256x256 arrays 
        containing the fitted E0, f and s parameter, their uncertainties and the fit type (1 for good fit; 0 for failed chi square test or unphysical 
        values of the fitted parameters, which is the case when E0, f or s is negative; -1 if the fit could not converge; and -2
        if the data is cut before fitting). For the 0, -1 and -2 cases (so, if not a good fit), the value of the parameters is 0.
        Zeros and NaN values are trimmed. Uses files named in the form <FPTfolder+FPTfilename+ST+'_THR_'+str(thr)+'.csv'>. 
        Saves files named in the form <FPTfolder+'<datatype>_ST_'+ST+'.csv'>, where datatype is E0matrix, fmatrix, smatrix, fittypematrix
        or fixedfittypematrix.

        Arguments:
        - FPTfolder (str): Folder path to the Flux per Pixel per Threshold arrays files
        (e.g. '/data/bfys/LucasvdH/VELOdata/coolsource/').
        - FPTfilename (str): Path to the Flux per Pixel per Threshold arrays made by fluxpixelthr (needs importing)
        (e.g. 'Fluxperpixel_Module0_VP0-1_ECS_data_ST_').
        - ST (str): Shutter time (e.g. '1s310ms').
        - fromThr (int): Starting threshold of the scan (e.g. 1500).
        - toThr (int): Ending threshold of the scan (e.g. 1600).
        - stepThr (int): Increment of the threshold between scans (e.g. 5).
        - p0 (list of length 4): initial guess of the parameters to fit [f, A, s, E0].  #A: recheck you do need some initial guess for these params, what is it set to and how do we find them??
        - nohits_cut (float between 0 and 1): Default=0.8. Mask the pixel if it has too many zero-valued points 
        (e.g. nohits_cut=0.8 means that the pixel will be masked if the it got no hits for 80% or more of the thresholds).
        - alpha (float): Default=0.05. Level of significance.
        '''

        #Loading FluxperASIC first to find indices of bad thresholds
        fluxASIC=np.loadtxt(FPTfolder+'FluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', dtype=float, delimiter=',')
        
        #If for some thresholds, all acquisitions were bad (which returns nan) set a index to slice all pertinent data to leave the data after the bad thresholds
        #If all thresholds were good, nan_cut is set to 0, which doesn't slice the data
    
        if any(np.isnan(fluxASIC)):  # finding where we have a value nan in the flux and getting the position value in the array so that later we can cut that and start ignoring those
            nan_cut=np.where(np.isnan(fluxASIC))[0][-1]+1
           
        else:
            nan_cut=0 #if all thresholds were good, nan_cut is set to 0
        
        #Loading Fluxperpixel and UncertaintyFluxperpixel files. NaN cut is applied by cutting Thrs
        Thrs=list(range(fromThr, toThr+1, stepThr))[nan_cut:]
        FPT=[]
        unc_FPT=[]

        totalflux_matrix = np.zeros((256,256)) #this is the total amount of flux per pixel summed over all thresholds
        for thr in Thrs:
            f=np.loadtxt(FPTfolder+'Fluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
            u=np.loadtxt(FPTfolder+'UncertaintyFluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
            FPT.append(f)
            unc_FPT.append(u)

            totalflux_matrix += f #summes up the flux per pixel over all thresholds
        
        # Here matrices are created for th evalue and its uncertainties from the parameters of the chose flux equation
        param_matrices = {param + 'matrix': np.zeros((256, 256)) for param in params_l}
        uparam_matrices = {'u'+uparam+'matrix': np.zeros((256, 256)) for uparam in params_l}
        print('param_matrices', len(param_matrices))
        
        #Assigning a 256x256 matrix to the fitted parameters. The 0s are to be replaced by the calculated values. If no or bad fits, the element is unchanged
        #E0matrix, fmatrix, smatrix= np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256))
        #Assigning a 256x256 matrix to the uncertainties in the fitted parameters
        #uE0matrix, ufmatrix, usmatrix= np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256))
        #Assigning a 256x256 matrix to fit type. These will express the type (lack of hits (-2), unable to find (-1), bad (0) or good (1) fits) of fit performed on the pixel. One for free fit and one for E0, f and s fixed fit
        fittypematrix = np.zeros((256, 256)) # first one is for the free fit and the second one is for the fixed fit.
        

        if FPTfolder_baseline != None: #if baseline folder is not none, import the baseline matrix
            baseline_matrix = self.baseline_importer(FPTfolder_baseline)
            ubaseline_matrix = self.baseline_importer(FPTfolder_baseline, width = True)
            baseline_ASIC = int(np.mean(baseline_matrix))
        else:
            baseline_matrix = np.zeros((256,256)) #if baseline folder is none, set baseline matrix to 0
            ubaseline_matrix = np.zeros((256,256))
            baseline_ASIC = int(0)


        
        fluxASIC=fluxASIC[nan_cut:] 
        unc_fluxASIC=np.loadtxt(FPTfolder+'UncertaintyFluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', dtype=float, delimiter=',')[nan_cut:]
        
        # Subtract the constant baseline value from every term in the list
        Thrs_ASIC = [thr - baseline_ASIC for thr in Thrs]
        
        poptASIC, pcovASIC=curve_fit(self.fitfunction, Thrs_ASIC, fluxASIC, sigma=unc_fluxASIC, p0=p0) #Performing fit for flux on ASIC
        poptASIC[1]=poptASIC[1]/65352 #To rescale to the flux on a single pixel. Used to calculate X2ASIC chi squared test
        E0ASIC=poptASIC[-1] #Renaming for reading ease

        analysistime=[] #to make an estimation of total analysis time
        
        #Looping over all the pixels
        good_pixel_counter = 0
        for row in range(256):
            st=time.time()
            print('\nPerfoming the fit for the pixels on row {}'.format(row))
            for column in range(256): # Below here, a pixel with a row and column is used.

                fluxpixel=[i[row][column] for i in FPT] #obtaining the individual pixel flux values
                unc_fluxpixel=[i[row][column] for i in unc_FPT]  #obtaining the individual pixel flux uncertainty values

                #If fluxpixel contains too many zeros, don't analyze it <-----
                if (len(fluxpixel)-np.count_nonzero(fluxpixel))/len(fluxpixel)>=nohits_cut:
                    #Assign -2 for type of fit if cut for lack of hits
                    fittypematrix[row][column]=-2
                    
                    continue
                
                Thrs_pixel = [thr - baseline_matrix[row][column] for thr in Thrs] #subtracting the baseline from the pixel, will only happen if the baseline is not containing zeros.
                               
                trimcut=np.nonzero(fluxpixel) #returns an array/matrix locating the indices of nonzero values, cuts the hanging zeros on the left and right
          
                thrs=Thrs_pixel[trimcut[0][0]:trimcut[0][-1]+1]  #using the location of non zeros values to get the Thrs at which we have an nonzero values
                fluxpixel=np.trim_zeros(fluxpixel) #same as above
                unc_fluxpixel=unc_fluxpixel[trimcut[0][0]:trimcut[0][-1]+1] #same as above
            

                ###Performing free fit
                try:  #changed here self.fitfunction to fit_function
                    popt, pcov=curve_fit(self.fitfunction, thrs, fluxpixel, sigma=unc_fluxpixel, p0=p0) #here taking into account sigma = uncertainty of fluxpixel and p0 = initial guess of the parameters
                except RuntimeError: 
                    #If fit could not converge, assign -1 to fit type
                    fittypematrix[row][column]=-1
                    continue

                except TypeError: 
                    #If fit could not converge, assign -1 to fit type
                    fittypematrix[row][column]=-1 #this was a -3 to identfiy those pixels which actually gave this specific error but changed to -1 since now are just fit not found
       
                    continue

                
                #Setting expected values and calculating X2
                E=self.fitfunction(thrs, *popt) #expected values from the fit
                X2=self.chisquare(np.array(fluxpixel), E, np.array(unc_fluxpixel)) #this is the chisquared defined within this class
                
                
                #Degrees of freedom for critical X2 (fitting four parameters, so -4)
                dof=len(fluxpixel)-4  #length of the fluxpixel datapoints minus the number of parameters (4 for the fit function or 5 for other fitfunctions)
                criticalX2=scipy.stats.chi2.ppf(1-alpha, dof) #this is the critical X2 value, which depends explicitly on the confidence parameter alpha

                #Mask pixel if X2>criticalX2, X2==0, or E0, f, s are negative (by leaving the zero in the matrix)
                #If calculating X2 was not possible (division by zero, X2=0), then consider the fit bad
                if X2<criticalX2 and X2!=0 and popt[-1]>0 and popt[0]>0 and popt[1]>0 and popt[2]>0:
                    #If the fit was good, change the 0 in the 256x256 matrices by the fitted value and its uncertainty
                    
                    for parameter in range(len(params_l)): # saving the values found by the fitting the individual pixels
                        key = params_l[parameter]
                       
                        param_matrices[key+'matrix'][row][column]=popt[parameter]
                        uparam_matrices['u'+key+'matrix'][row][column]=np.sqrt(np.diag(pcov))[parameter]
                    

                    #Assign 1 to fit type if it is a good fit (passed chi square test and parameters are positive)
                    fittypematrix[row][column]=1

                    good_pixel_counter += 1
                
            
            #Printing the time that analyzing the row took and estimated time (based on average time per row) to finish
            et=time.time()
            analysistime.append(et-st)
            
            print('Analyzing this row took {:.2f} seconds,\nEstimated time to finish analyzing the data set is {:.2f} minutes'.format(et-st, sum(analysistime)/len(analysistime)*(255-row)/60))
        
        print('good_pixel_counter', good_pixel_counter)

        

        for prefix, matrix in param_matrices.items():
            if FPTfolder_baseline != None and prefix == 'E0matrix':

                prefix = 'targetmatrix'
                np.savetxt(FPTfolder + f'{prefix}_ST_{ST}_{label}.csv', matrix, delimiter=',')
                print(f'\n{prefix} has been saved in {FPTfolder}{prefix}_ST_{ST}_{label}.csv')
            else:
                np.savetxt(FPTfolder + f'{prefix}_ST_{ST}_{label}.csv', matrix, delimiter=',')
                print(f'\n{prefix} has been saved in {FPTfolder}{prefix}_ST_{ST}_{label}.csv')

        for prefix, matrix in uparam_matrices.items():
            if FPTfolder_baseline != None and prefix == 'uE0matrix':

                prefix = 'utargetmatrix'
                np.savetxt(FPTfolder + f'{prefix}_ST_{ST}_{label}.csv', matrix, delimiter=',')
                print(f'\n{prefix} has been saved in {FPTfolder}{prefix}_ST_{ST}_{label}.csv')
            else:
                np.savetxt(FPTfolder + f'{prefix}_ST_{ST}_{label}.csv', matrix, delimiter=',')
                print(f'\n{prefix} has been saved in {FPTfolder}{prefix}_ST_{ST}_{label}.csv')

        np.savetxt(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv', fittypematrix, delimiter=',')
        print('\nFittypematrix has been saved in {}'.format(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv'))
       

        print('\nAll matrices saved for {} {}'.format(FPTfolder, label))

    
    def fitter2(self,FPTfolder, FPTfolder_baseline, ST, fromThr, toThr, stepThr, label, p0,nohits_cut=0.8, alpha=0.05):
        '''
        This one is fitter. should be similar to that of Lucas but  most likely containing corrections or additions to the code.

        Finds best fit to the flux per pixel and total flux on pixel grid to fitfunction (defined below). Saves eight 256x256 arrays 
        containing the fitted E0, f and s parameter, their uncertainties and the fit type (1 for good fit; 0 for failed chi square test or unphysical 
        values of the fitted parameters, which is the case when E0, f or s is negative; -1 if the fit could not converge; and -2
        if the data is cut before fitting). For the 0, -1 and -2 cases (so, if not a good fit), the value of the parameters is 0.
        Zeros and NaN values are trimmed. Uses files named in the form <FPTfolder+FPTfilename+ST+'_THR_'+str(thr)+'.csv'>. 
        Saves files named in the form <FPTfolder+'<datatype>_ST_'+ST+'.csv'>, where datatype is E0matrix, fmatrix, smatrix, fittypematrix
        or fixedfittypematrix.

        Arguments:
        - FPTfolder (str): Folder path to the Flux per Pixel per Threshold arrays files
        (e.g. '/data/bfys/LucasvdH/VELOdata/coolsource/').
        - FPTfilename (str): Path to the Flux per Pixel per Threshold arrays made by fluxpixelthr (needs importing)
        (e.g. 'Fluxperpixel_Module0_VP0-1_ECS_data_ST_').
        - ST (str): Shutter time (e.g. '1s310ms').
        - fromThr (int): Starting threshold of the scan (e.g. 1500).
        - toThr (int): Ending threshold of the scan (e.g. 1600).
        - stepThr (int): Increment of the threshold between scans (e.g. 5).
        - p0 (list of length 4): initial guess of the parameters to fit [f, A, s, E0].  #A: recheck you do need some initial guess for these params, what is it set to and how do we find them??
        - nohits_cut (float between 0 and 1): Default=0.8. Mask the pixel if it has too many zero-valued points 
        (e.g. nohits_cut=0.8 means that the pixel will be masked if the it got no hits for 80% or more of the thresholds).
        - alpha (float): Default=0.05. Level of significance.
        '''

        #Loading FluxperASIC first to find indices of bad thresholds
        #print(FPTfolder+'FluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv')
        fluxASIC=np.loadtxt(FPTfolder+'FluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', dtype=float, delimiter=',')
        
        #If for some thresholds, all acquisitions were bad (which returns nan) set a index to slice all pertinent data to leave the data after the bad thresholds
        #If all thresholds were good, nan_cut is set to 0, which doesn't slice the data
    
        if any(np.isnan(fluxASIC)):  #A: recheck this is just finding where we have a value nan in the flux and getting the position value in the array so that later we can cut that and start ignoring those
            nan_cut=np.where(np.isnan(fluxASIC))[0][-1]+1
           
        else:
            nan_cut=0
        
        #Loading Fluxperpixel and UncertaintyFluxperpixel files. NaN cut is applied by cutting Thrs
        Thrs=list(range(fromThr, toThr+1, stepThr))[nan_cut:]
        FPT=[]
        unc_FPT=[]
        for thr in Thrs:
            f=np.loadtxt(FPTfolder+'Fluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
            u=np.loadtxt(FPTfolder+'UncertaintyFluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
            FPT.append(f)
            unc_FPT.append(u)
        

        #Assigning a 256x256 matrix to the fitted parameters. The 0s are to be replaced by the calculated values. If no or bad fits, the element is unchanged
        E0matrix, fmatrix, smatrix, s2matrix= np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256))
        #Assigning a 256x256 matrix to the uncertainties in the fitted parameters
        uE0matrix, ufmatrix, usmatrix, us2matrix= np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256))
        #Assigning a 256x256 matrix to fit type. These will express the type (lack of hits (-2), unable to find (-1), bad (0) or good (1) fits) of fit performed on the pixel. One for free fit and one for E0, f and s fixed fit
        fittypematrix, fixedfittypematrix= np.zeros((256, 256)), np.zeros((256, 256)) #A: recheck, first one is for the free fit and the second one is for the fixed fit.
        
        if FPTfolder_baseline != None:
            baseline_matrix = self.baseline_importer(FPTfolder_baseline)
            ubaseline_matrix = self.baseline_importer(FPTfolder_baseline, width = True)
            baseline_ASIC = int(np.mean(baseline_matrix))
        else:
            baseline_matrix = np.zeros((256,256))
            ubaseline_matrix = np.zeros((256,256))
            baseline_ASIC = int(0)


        #Performing fit for flux on ASIC, A: recheck, this being the fit on the whole pixel grid.
        fluxASIC=fluxASIC[nan_cut:] 
        unc_fluxASIC=np.loadtxt(FPTfolder+'UncertaintyFluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', dtype=float, delimiter=',')[nan_cut:]
        
        # Subtract the constant baseline value from every term in the list
        Thrs_ASIC = [thr - baseline_ASIC for thr in Thrs]
        
        poptASIC, pcovASIC=curve_fit(self.fitfunction2AB, Thrs_ASIC, fluxASIC, sigma=unc_fluxASIC, p0=p0)
        poptASIC[1]=poptASIC[1]/65352 #To rescale to the flux on a single pixel. Used to calculate X2ASIC  #A: recheck X2 refers to chi sqaured test so not to be confused
        E0ASIC=poptASIC[-1] #Renaming for reading ease

        analysistime=[] #Just to make an estimation of total analysis time
        #A: recheck here starting with the individual pixel analysis to perform the fit individually
        #Selecting a single pixel
        for row in range(256):
            st=time.time()
            print('\nPerfoming the fit for the pixels on row {}'.format(row))
            for column in range(256):
                fluxpixel=[i[row][column] for i in FPT] #obtaining the individual values from the FPT  #A: recheck is it making an array since it is inside []? try it out
                unc_fluxpixel=[i[row][column] for i in unc_FPT] #A: recheck here it is obtaining all the fluxes of the same pixel for all Thrs adn adding it to an array to further analyse it.

                #If fluxpixel contains too many zeros, don't analyze it <-----
                if (len(fluxpixel)-np.count_nonzero(fluxpixel))/len(fluxpixel)>=nohits_cut:
                    #Assign -2 for type of fit if cut for lack of hits
                    fittypematrix[row][column]=-2
                    fixedfittypematrix[row][column]=-2
                    continue
                
                Thrs_pixel = [thr - baseline_matrix[row][column] for thr in Thrs]
                #Cut the hanging zeros on the left and right
                trimcut=np.nonzero(fluxpixel)  #A: recheck nonzero?? returns an array/matrix locating the indices of nonzero values
                thrs=Thrs_pixel[trimcut[0][0]:trimcut[0][-1]+1] #A: recheck, using the locastion of non zeros values to get the Thrs at which we have an nonzero values
                fluxpixel=np.trim_zeros(fluxpixel) #A: recheck, pretty much doing the same thing but this time just trimming the zzeros not needed to use the location as above
                unc_fluxpixel=unc_fluxpixel[trimcut[0][0]:trimcut[0][-1]+1]
            

                ###Performing fit  #A: recheck i think this is free fit?
                try:  
                    popt, pcov=curve_fit(self.fitfunction2AB, thrs, fluxpixel, sigma=unc_fluxpixel, p0=p0) #here taking into account sigma = 
                except RuntimeError: 
                    #If fit could not converge, assign -1 to fit type
                    fittypematrix[row][column]=-1
                    continue

                except TypeError: 
                    #If fit could not converge, assign -1 to fit type
                    fittypematrix[row][column]=-1 #this was a -3 to identfiy thos epixels which actually gave this specific error but changed to -1 since now are just fit not found
                    #print('TypeError')
                    continue

                
                #Setting expected values and calculating X2
                E=self.fitfunction2AB(thrs, *popt)
                X2=self.chisquare(np.array(fluxpixel), E, np.array(unc_fluxpixel)) #A: recheck, this is the chisquared defined by us as a function which is th eone described in the paper, nothing imported from np or scipy
                
                
                #Degrees of freedom for critical X2 (fitting four parameters, so -4)
                dof=len(fluxpixel)-5  #A: recheck why si the len of fluxpixel the number of dof??? RECHECK THS AGAIN I SAW SMTH ON A PAPER TALKING BOUT THIS
                criticalX2=scipy.stats.chi2.ppf(1-alpha, dof) #A: recheck, as said in the paper, this depends explicitly on the confidence parameter alpha

                #Mask pixel if X2>criticalX2, X2==0, or E0, f, s are negative (by leaving the zero in the matrix)
                #If calculating X2 was not possible (division by zero, X2=0), then consider the fit bad
                if X2<criticalX2 and X2!=0 and popt[-1]>0 and popt[0]>0 and popt[1]>0 and popt[2]>0: #A: recheck, criticalX2 is the chi squared value calculated with th enp function to be compared with the calculaated with our analysis
                    #If the fit was good, change the 0 in the 256x256 matrices by the fitted value and its uncertainty
                    #A: recheck, here only cdhanging the values in the corresponding matrices whenever it passes the chi2 test and the values have phsyical significance ie, positive
                    E0matrix[row][column]=popt[-1]
                    uE0matrix[row][column]=np.sqrt(np.diag(pcov))[-1]
                    fmatrix[row][column]=popt[0]
                    ufmatrix[row][column]=np.sqrt(np.diag(pcov))[0]
                    smatrix[row][column]=popt[2]
                    usmatrix[row][column]=np.sqrt(np.diag(pcov))[2]

                    s2matrix[row][column]=popt[3]
                    us2matrix[row][column]=np.sqrt(np.diag(pcov))[3]


                    #Assign 1 to fit type if it is a good fit (passed chi square test and parameters are positive)
                    fittypematrix[row][column]=1
                    

                
                ###Performing fit fixing E0, f, s  #A: recheck FIXED FIT
                #Fixing E0, f, s parameters in Fitter().fitfunction  
                def E0fsfixedfitfunction(x, A):
                    return self.fitfunction2AB(x, poptASIC[0], A, poptASIC[2],poptASIC[3], E0ASIC)

                try:
                    poptfixed, pcovfixed=curve_fit(E0fsfixedfitfunction, thrs, fluxpixel, sigma=unc_fluxpixel, p0=p0[1])
                except RuntimeError:
                    #If fit could not converge, assign -1 to fit type
                    fixedfittypematrix[row][column]=-1
                    continue

                except TypeError: 
                    #If fit could not converge, assign -1 to fit type
                    fittypematrix[row][column]=-1 #this was a -3 to identfiy thos epixels which actually gave this specific error but changed to -1 since now are just fit not found -3
                    print('TypeError')
                    continue
                
                #Setting expected values and calculating X2
                Efixed=E0fsfixedfitfunction(thrs, *poptfixed)
                X2fixed=self.chisquare(np.array(fluxpixel), Efixed, np.array(unc_fluxpixel))

                #Degrees of freedom for critical X2 (fitting one parameter, so -1)
                doffixed=len(fluxpixel)-1
                criticalX2fixed=scipy.stats.chi2.ppf(1-alpha, doffixed)
                
                if X2fixed<criticalX2fixed and X2fixed!=0:
                    #Assign 1 to fit type if it is a good fit (passed chi square test)  #A: recheck other bad fit conditions rewrite the matrix entries above
                    fixedfittypematrix[row][column]=1
                
            
            #Printing the time that analyzing the row took and estimated time (based on average time per row) to finish
            et=time.time()
            analysistime.append(et-st)
            print('Analyzing this row took {:.2f} seconds\nEstimated time to finish analyzing the data set is {:.2f} minutes'.format(et-st, sum(analysistime)/len(analysistime)*(255-row)/60))
        

        #Saving the eight matrices
        if FPTfolder_baseline != None:
            np.savetxt(FPTfolder+'targetmatrix_ST_'+ST+'_'+str(label)+'_fitter.csv', E0matrix, delimiter=',')
            np.savetxt(FPTfolder+'utargetmatrix_ST_'+ST+'_'+str(label)+'_fitter.csv', uE0matrix, delimiter=',')
        else:
            np.savetxt(FPTfolder+'E0matrix_ST_'+ST+'_'+str(label)+'.csv', E0matrix, delimiter=',')
            np.savetxt(FPTfolder+'uE0matrix_ST_'+ST+'_'+str(label)+'.csv', uE0matrix, delimiter=',')
        
        np.savetxt(FPTfolder+'fmatrix_ST_'+ST+'_'+str(label)+'.csv', fmatrix, delimiter=',')
        np.savetxt(FPTfolder+'ufmatrix_ST_'+ST+'_'+str(label)+'.csv', ufmatrix, delimiter=',')
        np.savetxt(FPTfolder+'smatrix_ST_'+ST+'_'+str(label)+'.csv', smatrix, delimiter=',')
        np.savetxt(FPTfolder+'usmatrix_ST_'+ST+'_'+str(label)+'.csv', usmatrix, delimiter=',')
        np.savetxt(FPTfolder+'s2matrix_ST_'+ST+'_'+str(label)+'.csv', s2matrix, delimiter=',')
        np.savetxt(FPTfolder+'us2matrix_ST_'+ST+'_'+str(label)+'.csv', us2matrix, delimiter=',')
        np.savetxt(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv', fittypematrix, delimiter=',')
        np.savetxt(FPTfolder+'Fixedfittypematrix_ST_'+ST+'_'+str(label)+'.csv', fixedfittypematrix, delimiter=',')

        #Printing where the matrices have been saved to make it easier to find them
        if FPTfolder_baseline != None:
            print('\ntargetmatrix has been saved in {}'.format(FPTfolder+'targetmatrix_ST_'+ST+'_'+str(label)+'_fitter.csv'))
            print('\nutargetmatrix has been saved in {}'.format(FPTfolder+'utargetmatrix_ST_'+ST+'_'+str(label)+'_fitter.csv'))
        else:
            print('\nE0matrix has been saved in {}'.format(FPTfolder+'E0matrix_ST_'+ST+'_'+str(label)+'.csv'))
            print('\nuE0matrix has been saved in {}'.format(FPTfolder+'uE0matrix_ST_'+ST+'_'+str(label)+'.csv'))
    
        print('\nfmatrix has been saved in {}'.format(FPTfolder+'fmatrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nufmatrix has been saved in {}'.format(FPTfolder+'ufmatrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nsmatrix has been saved in {}'.format(FPTfolder+'smatrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nusmatrix has been saved in {}'.format(FPTfolder+'usmatrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nFittypematrix has been saved in {}'.format(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nFixedfittypematrix has been saved in {}'.format(FPTfolder+'Fixedfittypematrix_ST_'+ST+'_'+str(label)+'.csv'))

        print('All matrices saved for {} {}'.format(FPTfolder, label))
    

    def fitter4(self,FPTfolder, FPTfolder_baseline, ST, fromThr, toThr, stepThr, label, p0,nohits_cut=0.8, alpha=0.05):
        '''
        This one is fitter. should be similar to that of Lucas but  most likely containing corrections or additions to the code.

        Finds best fit to the flux per pixel and total flux on pixel grid to fitfunction (defined below). Saves eight 256x256 arrays 
        containing the fitted E0, f and s parameter, their uncertainties and the fit type (1 for good fit; 0 for failed chi square test or unphysical 
        values of the fitted parameters, which is the case when E0, f or s is negative; -1 if the fit could not converge; and -2
        if the data is cut before fitting). For the 0, -1 and -2 cases (so, if not a good fit), the value of the parameters is 0.
        Zeros and NaN values are trimmed. Uses files named in the form <FPTfolder+FPTfilename+ST+'_THR_'+str(thr)+'.csv'>. 
        Saves files named in the form <FPTfolder+'<datatype>_ST_'+ST+'.csv'>, where datatype is E0matrix, fmatrix, smatrix, fittypematrix
        or fixedfittypematrix.

        Arguments:
        - FPTfolder (str): Folder path to the Flux per Pixel per Threshold arrays files
        (e.g. '/data/bfys/LucasvdH/VELOdata/coolsource/').
        - FPTfilename (str): Path to the Flux per Pixel per Threshold arrays made by fluxpixelthr (needs importing)
        (e.g. 'Fluxperpixel_Module0_VP0-1_ECS_data_ST_').
        - ST (str): Shutter time (e.g. '1s310ms').
        - fromThr (int): Starting threshold of the scan (e.g. 1500).
        - toThr (int): Ending threshold of the scan (e.g. 1600).
        - stepThr (int): Increment of the threshold between scans (e.g. 5).
        - p0 (list of length 4): initial guess of the parameters to fit [f, A, s, E0].  #A: recheck you do need some initial guess for these params, what is it set to and how do we find them??
        - nohits_cut (float between 0 and 1): Default=0.8. Mask the pixel if it has too many zero-valued points 
        (e.g. nohits_cut=0.8 means that the pixel will be masked if the it got no hits for 80% or more of the thresholds).
        - alpha (float): Default=0.05. Level of significance.
        '''

        #Loading FluxperASIC first to find indices of bad thresholds
        #print(FPTfolder+'FluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv')
        fluxASIC=np.loadtxt(FPTfolder+'FluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', dtype=float, delimiter=',')
        
        #If for some thresholds, all acquisitions were bad (which returns nan) set a index to slice all pertinent data to leave the data after the bad thresholds
        #If all thresholds were good, nan_cut is set to 0, which doesn't slice the data
    
        if any(np.isnan(fluxASIC)):  #A: recheck this is just finding where we have a value nan in the flux and getting the position value in the array so that later we can cut that and start ignoring those
            nan_cut=np.where(np.isnan(fluxASIC))[0][-1]+1
           
        else:
            nan_cut=0
        
        #Loading Fluxperpixel and UncertaintyFluxperpixel files. NaN cut is applied by cutting Thrs
        Thrs=list(range(fromThr, toThr+1, stepThr))[nan_cut:]
        FPT=[]
        unc_FPT=[]
        for thr in Thrs:
            f=np.loadtxt(FPTfolder+'Fluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
            u=np.loadtxt(FPTfolder+'UncertaintyFluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
            FPT.append(f)
            unc_FPT.append(u)
        

        #Assigning a 256x256 matrix to the fitted parameters. The 0s are to be replaced by the calculated values. If no or bad fits, the element is unchanged
        E0matrix, fmatrix, smatrix= np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256))
        #Assigning a 256x256 matrix to the uncertainties in the fitted parameters
        uE0matrix, ufmatrix, usmatrix= np.zeros((256, 256)), np.zeros((256, 256)), np.zeros((256, 256))
        #Assigning a 256x256 matrix to fit type. These will express the type (lack of hits (-2), unable to find (-1), bad (0) or good (1) fits) of fit performed on the pixel. One for free fit and one for E0, f and s fixed fit
        fittypematrix, fixedfittypematrix= np.zeros((256, 256)), np.zeros((256, 256)) #A: recheck, first one is for the free fit and the second one is for the fixed fit.
        
        if FPTfolder_baseline != None:
            baseline_matrix = self.baseline_importer(FPTfolder_baseline)
            ubaseline_matrix = self.baseline_importer(FPTfolder_baseline, width = True)
            baseline_ASIC = int(np.mean(baseline_matrix))
        else:
            baseline_matrix = np.zeros((256,256))
            ubaseline_matrix = np.zeros((256,256))
            baseline_ASIC = int(0)


        #Performing fit for flux on ASIC, A: recheck, this being the fit on the whole pixel grid.
        fluxASIC=fluxASIC[nan_cut:] 
        unc_fluxASIC=np.loadtxt(FPTfolder+'UncertaintyFluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', dtype=float, delimiter=',')[nan_cut:]
        
        # Subtract the constant baseline value from every term in the list
        Thrs_ASIC = [thr - baseline_ASIC for thr in Thrs]
        
        poptASIC, pcovASIC=curve_fit(self.fitfunction4, Thrs_ASIC, fluxASIC, sigma=unc_fluxASIC, p0=p0)
        poptASIC[1]=poptASIC[1]/65352 #To rescale to the flux on a single pixel. Used to calculate X2ASIC  #A: recheck X2 refers to chi sqaured test so not to be confused
        E0ASIC=poptASIC[-1] #Renaming for reading ease

        analysistime=[] #Just to make an estimation of total analysis time
        #A: recheck here starting with the individual pixel analysis to perform the fit individually
        #Selecting a single pixel
        good_pixel_counter = 0
        for row in range(256):
            st=time.time()
            print('\nPerfoming the fit for the pixels on row {}'.format(row))
            for column in range(256):
                fluxpixel=[i[row][column] for i in FPT] #obtaining the individual values from the FPT  #A: recheck is it making an array since it is inside []? try it out
                unc_fluxpixel=[i[row][column] for i in unc_FPT] #A: recheck here it is obtaining all the fluxes of the same pixel for all Thrs adn adding it to an array to further analyse it.

                #If fluxpixel contains too many zeros, don't analyze it <-----
                if (len(fluxpixel)-np.count_nonzero(fluxpixel))/len(fluxpixel)>=nohits_cut:
                    #Assign -2 for type of fit if cut for lack of hits
                    fittypematrix[row][column]=-2
                    fixedfittypematrix[row][column]=-2
                    continue
                
                Thrs_pixel = [thr - baseline_matrix[row][column] for thr in Thrs]
                #Cut the hanging zeros on the left and right
                trimcut=np.nonzero(fluxpixel) #A: recheck nonzero?? returns an array/matrix locating the indices of nonzero values
          
                thrs=Thrs_pixel[trimcut[0][0]:trimcut[0][-1]+1] #A: recheck, using the locastion of non zeros values to get the Thrs at which we have an nonzero values
                fluxpixel=np.trim_zeros(fluxpixel) #A: recheck, pretty much doing the same thing but this time just trimming the zzeros not needed to use the location as above
                unc_fluxpixel=unc_fluxpixel[trimcut[0][0]:trimcut[0][-1]+1]
            

                ###Performing fit  #A: recheck i think this is free fit?
                try:  
                    popt, pcov=curve_fit(self.fitfunction4, thrs, fluxpixel, sigma=unc_fluxpixel, p0=p0) #here taking into account sigma = 
                except RuntimeError: 
                    #If fit could not converge, assign -1 to fit type
                    fittypematrix[row][column]=-1
                    continue

                except TypeError: 
                    #If fit could not converge, assign -1 to fit type
                    fittypematrix[row][column]=-1 #this was a -3 to identfiy thos epixels which actually gave this specific error but changed to -1 since now are just fit not found
                    #print('TypeError')
                    continue

                
                #Setting expected values and calculating X2
                E=self.fitfunction4(thrs, *popt)
                X2=self.chisquare(np.array(fluxpixel), E, np.array(unc_fluxpixel)) #A: recheck, this is the chisquared defined by us as a function which is th eone described in the paper, nothing imported from np or scipy
                
                
                #Degrees of freedom for critical X2 (fitting four parameters, so -4)
                dof=len(fluxpixel)-4  #A: recheck why si the len of fluxpixel the number of dof??? RECHECK THS AGAIN I SAW SMTH ON A PAPER TALKING BOUT THIS
                criticalX2=scipy.stats.chi2.ppf(1-alpha, dof) #A: recheck, as said in the paper, this depends explicitly on the confidence parameter alpha

                #Mask pixel if X2>criticalX2, X2==0, or E0, f, s are negative (by leaving the zero in the matrix)
                #If calculating X2 was not possible (division by zero, X2=0), then consider the fit bad
                if X2<criticalX2 and X2!=0 and popt[-1]>0 and popt[0]>0 and popt[1]>0 and popt[2]>0: #A: recheck, criticalX2 is the chi squared value calculated with th enp function to be compared with the calculaated with our analysis
                    #If the fit was good, change the 0 in the 256x256 matrices by the fitted value and its uncertainty
                    #A: recheck, here only cdhanging the values in the corresponding matrices whenever it passes the chi2 test and the values have phsyical significance ie, positive
                    E0matrix[row][column]=popt[-1]
                    uE0matrix[row][column]=np.sqrt(np.diag(pcov))[-1]
                    fmatrix[row][column]=popt[0]
                    ufmatrix[row][column]=np.sqrt(np.diag(pcov))[0]
                    smatrix[row][column]=popt[2]
                    usmatrix[row][column]=np.sqrt(np.diag(pcov))[2]

                    #Assign 1 to fit type if it is a good fit (passed chi square test and parameters are positive)
                    fittypematrix[row][column]=1

                    good_pixel_counter += 1
                    

                
                ###Performing fit fixing E0, f, s  #A: recheck FIXED FIT
                #Fixing E0, f, s parameters in Fitter().fitfunction  
                def E0fsfixedfitfunction(x, A):
                    return self.fitfunction4(x, poptASIC[0], A, poptASIC[2], E0ASIC)

                try:
                    poptfixed, pcovfixed=curve_fit(E0fsfixedfitfunction, thrs, fluxpixel, sigma=unc_fluxpixel, p0=p0[1])
                except RuntimeError:
                    #If fit could not converge, assign -1 to fit type
                    fixedfittypematrix[row][column]=-1
                    continue

                except TypeError: 
                    #If fit could not converge, assign -1 to fit type
                    fittypematrix[row][column]=-1 #this was a -3 to identfiy thos epixels which actually gave this specific error but changed to -1 since now are just fit not found -3
                    print('TypeError')
                    continue
                
                #Setting expected values and calculating X2
                Efixed=E0fsfixedfitfunction(thrs, *poptfixed)
                X2fixed=self.chisquare(np.array(fluxpixel), Efixed, np.array(unc_fluxpixel))

                #Degrees of freedom for critical X2 (fitting one parameter, so -1)
                doffixed=len(fluxpixel)-1
                criticalX2fixed=scipy.stats.chi2.ppf(1-alpha, doffixed)
                
                if X2fixed<criticalX2fixed and X2fixed!=0:
                    #Assign 1 to fit type if it is a good fit (passed chi square test)  #A: recheck other bad fit conditions rewrite the matrix entries above
                    fixedfittypematrix[row][column]=1
                
            
            #Printing the time that analyzing the row took and estimated time (based on average time per row) to finish
            et=time.time()
            analysistime.append(et-st)
            print('Analyzing this row took {:.2f} seconds\nEstimated time to finish analyzing the data set is {:.2f} minutes'.format(et-st, sum(analysistime)/len(analysistime)*(255-row)/60))
        
        print('good_pixel_counter', good_pixel_counter)
        #Saving the eight matrices
        if FPTfolder_baseline != None:
            np.savetxt(FPTfolder+'targetmatrix_ST_'+ST+'_'+str(label)+'.csv', E0matrix, delimiter=',')
            np.savetxt(FPTfolder+'utargetmatrix_ST_'+ST+'_'+str(label)+'.csv', uE0matrix, delimiter=',')
        else:
            np.savetxt(FPTfolder+'E0matrix_ST_'+ST+'_'+str(label)+'.csv', E0matrix, delimiter=',')
            np.savetxt(FPTfolder+'uE0matrix_ST_'+ST+'_'+str(label)+'.csv', uE0matrix, delimiter=',')
        
        np.savetxt(FPTfolder+'fmatrix_ST_'+ST+'_'+str(label)+'.csv', fmatrix, delimiter=',')
        np.savetxt(FPTfolder+'ufmatrix_ST_'+ST+'_'+str(label)+'.csv', ufmatrix, delimiter=',')
        np.savetxt(FPTfolder+'smatrix_ST_'+ST+'_'+str(label)+'.csv', smatrix, delimiter=',')
        np.savetxt(FPTfolder+'usmatrix_ST_'+ST+'_'+str(label)+'.csv', usmatrix, delimiter=',')
        np.savetxt(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv', fittypematrix, delimiter=',')
        np.savetxt(FPTfolder+'Fixedfittypematrix_ST_'+ST+'_'+str(label)+'.csv', fixedfittypematrix, delimiter=',')

        #Printing where the matrices have been saved to make it easier to find them
        if FPTfolder_baseline != None:
            print('\ntargetmatrix has been saved in {}'.format(FPTfolder+'targetmatrix_ST_'+ST+'_'+str(label)+'.csv'))
            print('\nutargetmatrix has been saved in {}'.format(FPTfolder+'utargetmatrix_ST_'+ST+'_'+str(label)+'.csv'))
        else:
            print('\nE0matrix has been saved in {}'.format(FPTfolder+'E0matrix_ST_'+ST+'_'+str(label)+'.csv'))
            print('\nuE0matrix has been saved in {}'.format(FPTfolder+'uE0matrix_ST_'+ST+'_'+str(label)+'.csv'))
    
        print('\nfmatrix has been saved in {}'.format(FPTfolder+'fmatrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nufmatrix has been saved in {}'.format(FPTfolder+'ufmatrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nsmatrix has been saved in {}'.format(FPTfolder+'smatrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nusmatrix has been saved in {}'.format(FPTfolder+'usmatrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nFittypematrix has been saved in {}'.format(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv'))
        print('\nFixedfittypematrix has been saved in {}'.format(FPTfolder+'Fixedfittypematrix_ST_'+ST+'_'+str(label)+'.csv'))

        print('All matrices saved for {} {}'.format(FPTfolder, label))
        
    




class Matrices(organizer.FileOrganizer, Fitter):
    def __init__(self,
            FPTfolder_baseline_l = None,
            s2_opt = None,
            datafolder_l=None,
            FPTfolder_l=None,
            FPTfilename_l=None,
            ST_l=None,
            nacq_l=None,
            fromThr_l=None,
            toThr_l=None,
            stepThr_l=None,
            label_l=None,
            bad_cut_l=None,
            condition_l=None,
            p0=None,
            p0_baseline=None,
            baseline_opt=None,
            separation_opt=None,
            predict_opt_l=None,
            output_opt=None,
            fit_function = None
        ):
        
        """
        
        """
        self.FPTfolder_baseline_l = FPTfolder_baseline_l
        self.s2_opt = s2_opt
        self.datafolder_l = datafolder_l
        self.FPTfolder_l = FPTfolder_l
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

        
        print('Matrices FPTfolder_l', FPTfolder_l.split('/')[-2])
        results_type = FPTfolder_l.split('/')[-2]

        

        

        
    def target_finder(self, FPTfolder_baseline=None, label_l=None, ST_l=None, FPTfolder_l=None):
        
        FPTfolder_baseline = self.FPTfolder_baseline_l
        label_l = self.label_l
        ST = self.ST_l
        FPTfolder_l = self.FPTfolder_l

        if FPTfolder_baseline != None:
            baseline_matrix = Fitter().baseline_importer(FPTfolder_baseline)
            ubaseline_matrix = Fitter().baseline_importer(FPTfolder_baseline, width = True)
            baseline_ASIC = int(np.mean(baseline_matrix))
        
        E0matrix = self.open_csv(prefix='E0', label = label_l, path = FPTfolder_l)[0]
        uE0matrix = self.open_csv(prefix='uE0', label = label_l, path = FPTfolder_l)[0]

        E0matrix=np.ma.masked_where(E0matrix<=0, E0matrix)
        uE0matrix=np.ma.masked_where(E0matrix<=0, uE0matrix)
        

        if np.mean(E0matrix) <= 200 and np.mean(E0matrix)>=90:
            print('Mean of E0 matrix seems to be a target matrix already...')
            targetmatrix = E0matrix
            utargetmatrix = uE0matrix
        else:
            targetmatrix = np.zeros((256,256))
            utargetmatrix = np.zeros((256,256))
            for row in range(256):
                for column in range(256):
                    #Ignoring the masked elements (otherwise will return -- instead of a number)
                    if np.ma.is_masked(E0matrix[row][column]) or baseline_matrix[row][column] == 0:
                        continue

                    targetmatrix[row][column] = E0matrix[row][column] - baseline_matrix[row][column]
                    utargetmatrix[row][column] = np.sqrt((uE0matrix[row][column])**2  + (ubaseline_matrix[row][column])**2) #using error propagation although might be wrong

            targetmatrix=np.ma.masked_where((baseline_matrix==0) | (targetmatrix<=0), targetmatrix)
            utargetmatrix=np.ma.masked_where((baseline_matrix==0) | (targetmatrix<=0), utargetmatrix)

        np.savetxt(FPTfolder_l+'targetmatrix_ST_'+ST+'_'+str(label_l)+'.csv', targetmatrix, delimiter=',')
        np.savetxt(FPTfolder_l+'utargetmatrix_ST_'+ST+'_'+str(label_l)+'.csv', utargetmatrix, delimiter=',')
        print('\ntargetmatrix has been saved in {}'.format(FPTfolder_l+'targetmatrix_ST_'+ST+'_'+str(label_l)+'.csv'))
        print('\nutargetmatrix has been saved in {}'.format(FPTfolder_l+'utargetmatrix_ST_'+ST+'_'+str(label_l)+'.csv'))
    
    def K_finder(self, FPTfolder_baseline = None,baselineASICbeforefit = None, s2_opt = None,FPTfolder = None, ST = None, label = None, E_xrays = 5900, E_ehp = 3.69, E_ehp_err = 0.11):

        FPTfolder = self.FPTfolder_l
        ST = self.ST_l
        label = self.label_l
        s2_opt = self.s2_opt
      
        n_ehp = E_xrays / E_ehp
        n_ehp_err = n_ehp * E_ehp_err / E_ehp
        
        #print('X-rays energy = {} eV, n_ehp = {}'.format(E_xrays, n_ehp))
        target_matrix = self.open_csv(prefix='target', label=label, path=FPTfolder)[0]
        utarget_matrix = self.open_csv(prefix='utarget', label=label, path=FPTfolder)[0]

        # had to add this condition for masking those value where it was 0, not sure why if this filtering already exists within target finder
        target_matrix=np.ma.masked_where(target_matrix<=0, target_matrix)
        utarget_matrix=np.ma.masked_where(target_matrix<=0, utarget_matrix)

        target_matrix=np.ma.masked_where(target_matrix>=1000, target_matrix)
        utarget_matrix=np.ma.masked_where(target_matrix>=1000, utarget_matrix)

        target_mean, target_std = self.find_mean_std(target_matrix, print_opt='target_matrix')
        
        find_ASIC_input = input('Do you want to calculate ASIC and ASICgood values? (y/n): ')
        if find_ASIC_input.lower() == 'y':
            targetASIC, utargetASIC, paramsASIC, poptASIC, poptuncASIC = self.find_ASIC(FPTfolder_baseline = FPTfolder_baseline, baselineASICbeforefit = baselineASICbeforefit, s2_opt=s2_opt)
            targetASICgood, utargetASICgood, paramsASICgood, poptASICgood, poptuncASICgood = self.find_ASICgood(FPTfolder_baseline = FPTfolder_baseline, baselineASICbeforefit = baselineASICbeforefit, s2_opt=s2_opt)
            
            KASIC = n_ehp/targetASIC
            uKASIC = KASIC*np.sqrt((n_ehp_err/n_ehp)**2 + (utargetASIC/targetASIC)**2)

            KASICgood =n_ehp/targetASICgood
            uKASICgood = KASICgood*np.sqrt((n_ehp_err/n_ehp)**2 + (utargetASICgood/targetASICgood)**2)
        else:
            paramsASIC = ['f', 'A', 's', 'E0'] #set as default but depends on the fit function being considered, change when the specific of the fit function is added to the data info file.
            poptASIC = [0,0,0,0]
            poptuncASIC = [0,0,0,0]

            paramsASICgood  = paramsASIC
            poptASICgood = poptASIC
            poptuncASICgood = poptuncASIC

            targetASIC, utargetASIC, paramsASIC, poptASIC, poptuncASIC = 0,0,paramsASIC, poptASIC, poptuncASIC
            targetASICgood, utargetASICgood, paramsASICgood, poptASICgood, poptuncASICgood = 0,0, paramsASICgood, poptASICgood, poptuncASICgood
            KASIC, uKASIC, KASICgood, uKASICgood= 0,0,0,0

    


        K_matrix = n_ehp / target_matrix
       
        K_matrix_err = K_matrix * np.sqrt((n_ehp_err / n_ehp) ** 2 + (utarget_matrix / target_matrix) ** 2)

        K_matrix=np.ma.masked_where(K_matrix>=50, K_matrix)
        K_matrix_err=np.ma.masked_where(K_matrix>=50, K_matrix_err)

        K_mean, K_std = self.find_mean_std(K_matrix, print_opt='K_matrix')
        
  

        print('K_mean, K_std, using meand of the K matrix', K_mean, K_std)
    
        K_matrix = np.ma.filled(K_matrix, 0)
        K_matrix_err = np.ma.filled(K_matrix_err, 0)
        
        flatten_K_matrix = K_matrix.flatten()
        flatten_K_matrix = flatten_K_matrix[flatten_K_matrix > 0]

        pixels_flatten_K_matrix = len(flatten_K_matrix) #getting final number of goodfit pixels after all the previous filtering
        K_mean_err = K_std/np.sqrt(pixels_flatten_K_matrix)

        np.savetxt(FPTfolder + 'K_matrix_ST_' + ST + '_' + str(label) + '.csv', K_matrix, delimiter=',')
        np.savetxt(FPTfolder + 'uK_matrix_ST_' + ST + '_' + str(label) + '.csv', K_matrix_err, delimiter=',')
        
        print('\nK_matrices has been saved in {}'.format(FPTfolder + 'K_matrix_ST_' + ST + '_' + str(label) + '.csv'))
        print('\nFor saving results in results_collection.csv:\n')

        self.results_csv_file(savepath = FPTfolder, baseline = self.baseline_opt, ST = ST, label = label, targetmean = target_mean, utargetmean = target_std,targetASIC=targetASIC, utargetASIC=utargetASIC,targetASICgood=targetASICgood, utargetASICgood=utargetASICgood, Kmean = K_mean, Kstd = K_std, uKmean = K_mean_err, KASIC=KASIC, uKASIC=uKASIC, KASICgood=KASICgood, uKASICgood=uKASICgood)
        self.log_file(FPTfolder, label, pixels_flatten_K_matrix,K_mean, K_std, K_mean_err, paramsASIC, poptASIC, poptuncASIC, paramsASICgood, poptASICgood, poptuncASICgood )


    def find_mean_std(self, matrix, print_opt = False):

        mean = np.ma.mean(matrix)
        std = np.ma.std(matrix)
        
        if print_opt:
            print('{} mean = {:.2f}, std = {:.2f}'.format(print_opt, mean, std))
        
        return mean, std
    
    
    def find_ASIC(self, baselineASICbeforefit = None, s2_opt = None, FPTfolder = None, ST = None, label = None, fromThr = None, toThr = None, stepThr = None, FPTfolder_baseline = None, p0 = None, ASIC_params_output = None, plot_output = None):
        print('Finding ASIC parameters...')

        s2_opt = self.s2_opt
        FPTfolder = self.FPTfolder_l
        ST = self.ST_l
        label = self.label_l
        fromThr = self.fromThr_l
        toThr = self.toThr_l
        stepThr = self.stepThr_l
        
        print('FPTfolder_baseline:', FPTfolder_baseline, 's2_opt:', s2_opt)
        fluxASIC=np.loadtxt(FPTfolder+'FluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', dtype=float, delimiter=',')
        
        #If for some thresholds, all acquisitions were bad (which returns nan) set a index to slice all pertinent data to leave the data after the bad thresholds
        #If all thresholds were good, nan_cut is set to 0, which doesn't slice the data
    
        if any(np.isnan(fluxASIC)):  #A: recheck this is just finding where we have a value nan in the flux and getting the position value in the array so that later we can cut that and start ignoring those
            nan_cut=np.where(np.isnan(fluxASIC))[0][-1]+1
           
        else:
            nan_cut=0

        if FPTfolder == '/data/bfys/apuicerc/N037_new/resultsultimate_newbaseline/' or FPTfolder == '/data/bfys/apuicerc/N037_new/resultscalibration_newbaseline/' :
            nan_cut = 22
        else:
            nan_cut = nan_cut
        
        print('nan_cut', nan_cut)
        
        #Loading Fluxperpixel and UncertaintyFluxperpixel files. NaN cut is applied by cutting Thrs
        Thrs=list(range(fromThr, toThr+1, stepThr))
        
        if FPTfolder_baseline != None:
            baseline_matrix = self.baseline_importer(FPTfolder_baseline)
            ubaseline_matrix = self.baseline_importer(FPTfolder_baseline, width = True)
            baseline_ASIC = self.baseline_importer(FPTfolder_baseline, ASIC_baseline=True)
            ubaselineASIC = self.baseline_importer(FPTfolder_baseline, width = True)
        elif FPTfolder == '/data/bfys/apuicerc/N037_new/resultsultimate_newbaseline/' or FPTfolder == '/data/bfys/apuicerc/N037_new/resultscalibration_newbaseline/' :
            baseline_matrix = np.zeros((256,256))
            ubaseline_matrix = np.zeros((256,256))
            baseline_ASIC = int(0)
            ubaselineASIC = int(0)
            
        else:
            baseline_matrix = np.zeros((256,256))
            ubaseline_matrix = np.zeros((256,256))
            baseline_ASIC = int(0)
            ubaselineASIC = int(0)

        if baselineASICbeforefit != None:
            Thrs_ASIC = [thr - baseline_ASIC for thr in Thrs]
            p0_used = self.p0_baseline
        else:
            Thrs_ASIC = [thr - 0 for thr in Thrs]
            p0_used = self.p0

        #print('p0:', p0_used, 'Thrs_ASIC range',min(Thrs_ASIC), max(Thrs_ASIC), 'baselineASICbeforefit:', baselineASICbeforefit, 's2_opt:', s2_opt, 'FPTfolder', FPTfolder)

        #Performing fit for flux on ASIC, A: recheck, this being the fit on the whole pixel grid.

        fluxASIC=fluxASIC[nan_cut:] 


        unc_fluxASIC=np.loadtxt(FPTfolder+'UncertaintyFluxperASIC_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(fromThr)+'-'+str(toThr)+'-step'+str(stepThr)+'.csv', dtype=float, delimiter=',')[nan_cut:]
        Thrs_ASIC = Thrs_ASIC[nan_cut:]
        # Assuming unc_fluxASIC is a NumPy array or a list
        if all(value == 0 for value in unc_fluxASIC):
            print("unc_fluxASIC contains only zeros, fitting without considering unc_fluxASIC")
            unc_fluxASIC = None
        else:
            print("unc_fluxASIC contains non-zero values")


        signature = inspect.signature(self.fit_function)
        params_l = list(signature.parameters)[1:] 
        n_parameters = len(params_l)

 
        poptASIC, pcovASIC=curve_fit(self.fit_function, Thrs_ASIC, fluxASIC, sigma=unc_fluxASIC, p0=p0_used)
        paramsASIC = params_l

        #poptASIC[1]=poptASIC[1]/65352 #To rescale to the flux on a single pixel. Used to calculate X2ASIC  #A: recheck X2 refers to chi sqaured test so not to be confused
        poptuncASIC=np.sqrt(np.diag(pcovASIC))
        E0ASIC=poptASIC[-1] #Renaming for reading ease
        uE0ASIC = poptuncASIC[-1]


        #print('Thrs_ASIC, fluxASIC, unc_fluxASIC, poptASIC, poptuncASIC')
        #print(Thrs_ASIC, fluxASIC, unc_fluxASIC, poptASIC, poptuncASIC)


        if plot_output != None:
            
            return Thrs_ASIC, fluxASIC, unc_fluxASIC, poptASIC, poptuncASIC
    
        if baselineASICbeforefit != None: #if basleine is considere before fitting then this says its E0 but it is really target
            targetASIC = E0ASIC
            utargetASIC = uE0ASIC

        elif FPTfolder == '/data/bfys/apuicerc/N037_new/resultsultimate_newbaseline/' or FPTfolder == '/data/bfys/apuicerc/N037_new/resultscalibration_newbaseline/' :
            targetASIC = E0ASIC
            utargetASIC = uE0ASIC

        else: #we still want to implement the baselineASIC to find targetASIC, AFTER the computation of E0ASIC
            targetASIC = E0ASIC - baseline_ASIC
            utargetASIC = np.sqrt(uE0ASIC**2 + np.mean(ubaselineASIC[ubaselineASIC>=0])**2)
            
            
            
        print(paramsASIC)
        print(poptASIC)
        print(poptuncASIC)
        params_header = 'poptASIC, Value ,Uncertainty\n'
        params_data = ''
        for i in range(len(poptASIC)):
            params_data += '{},{},{}\n'.format(paramsASIC[i], poptASIC[i], poptuncASIC[i])

        params_data += 'target,{},{}\n'.format(targetASIC, utargetASIC)
        with open(os.path.join(FPTfolder, 'ASIC_params'+str(label)+'.csv'), 'w') as f:
            f.write(params_header)
            f.write(params_data)
        print('ASIC params saved in {}'.format(os.path.join(FPTfolder, 'ASIC_params'+str(label)+'.csv')))
       
        return targetASIC, utargetASIC, paramsASIC, poptASIC, poptuncASIC


    def find_ASICgood(self, baselineASICbeforefit = None, FPTfolder = None, ST = None, label = None, fromThr = None, toThr = None, stepThr = None, FPTfolder_baseline = None, p0 = None, s2_opt = None, plot_output = None):
        s2_opt = self.s2_opt
        FPTfolder = self.FPTfolder_l
        ST = self.ST_l
        label = self.label_l
        fromThr = self.fromThr_l
        toThr = self.toThr_l
        stepThr = self.stepThr_l

        print('Finding ASICgood parameters...')
    
            #here fittyoe is imported to obtain the pixels that found a fit.
        fittype=np.loadtxt(FPTfolder+'Fittypematrix_ST_'+ST+'_'+str(label)+'.csv', dtype=float, delimiter=',')
            #A: recheck these are the matrices (one for free fit and another for fixed fit) containing all the info about the goodness ofthe fit that being
            #wether the values 1 for good fit, 0 -1 and i think -2 too if it could not be found, unphysical solution or computational error maybe recheck this!
        #Loading FluxperASIC first to find indices of bad thresholds
        totalgoodpixels = np.count_nonzero(fittype == 1)
     
        F_ASICgood=np.zeros((int((toThr-fromThr)/stepThr+1), 1))
        unc_F_ASICgood=np.zeros((int((toThr-fromThr)/stepThr+1), 1))

        for thr in range(fromThr, toThr+1, stepThr):
            fluxperpixel=np.loadtxt(FPTfolder+'Fluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
            unc_fluxperpixel=np.loadtxt(FPTfolder+'UncertaintyFluxperpixel_Module0_VP3-1_ECS_data_ST_'+ST+'_'+str(label)+'_THR_'+str(thr)+'.csv', dtype=float, delimiter=',')
            F_pixelgooddummy = []

            avg_dummy = 1
            for row in range(256):
                for column in range(256):
                    if fittype[row][column] == 1:
                        F_pixelgooddummy.append(fluxperpixel[row][column])
                        
                        avg_dummy += 1

            mean = sum(F_pixelgooddummy)/totalgoodpixels
            squared_diff_sum = np.sum((x - mean)**2 for x in F_pixelgooddummy)
            
            acqtime = 2.621
            F_ASICgood[int((thr-fromThr)/stepThr)][0] = mean #(totalgoodpixels)
            unc_F_ASICgood[int((thr-fromThr)/stepThr)][0] = np.sqrt(squared_diff_sum)/totalgoodpixels    #np.sqrt(squared_diff_sum/avg_dummy)/(acqtime*256*256)


        #If for some thresholds, all acquisitions were bad (which returns nan) set a index to slice all pertinent data to leave the data after the bad thresholds
        #If all thresholds were good, nan_cut is set to 0, which doesn't slice the data
        if any(np.isnan(F_ASICgood)):  #A: recheck this is just finding where we have a value nan in the flux and getting the position value in the array so that later we can cut that and start ignoring those
            nan_cut=np.where(np.isnan(F_ASICgood))[0][-1]+1
            #print('nan_cut:', nan_cut)

        else:
            nan_cut=0
        
        if FPTfolder == '/data/bfys/apuicerc/N037_new/resultsultimate_newbaseline/' or FPTfolder == '/data/bfys/apuicerc/N037_new/resultscalibration_newbaseline/' :
            nan_cut = 22
        else:
            nan_cut = nan_cut
        
        #nan_cut = -9 not sure ont he re4ason why -9 is here but when i set to prinbt the nancut in the if condition above it doesnt seem to doit which means that the nan_cut being 0 shoul dbe fine.

        Thrs=list(range(fromThr, toThr+1, stepThr))
        
        if FPTfolder_baseline != None:
            baseline_matrix = self.baseline_importer(FPTfolder_baseline)
            ubaseline_matrix = self.baseline_importer(FPTfolder_baseline, width = True)
            baseline_ASIC = self.baseline_importer(FPTfolder_baseline, ASIC_baseline=True)
            ubaselineASIC = self.baseline_importer(FPTfolder_baseline, width = True)
        elif FPTfolder == '/data/bfys/apuicerc/N037_new/resultsultimate_newbaseline/' or FPTfolder == '/data/bfys/apuicerc/N037_new/resultscalibration_newbaseline/' :
            baseline_matrix = np.zeros((256,256))
            ubaseline_matrix = np.zeros((256,256))
            baseline_ASIC = int(0)
            ubaselineASIC = int(0)
            
        else:
            baseline_matrix = np.zeros((256,256))
            ubaseline_matrix = np.zeros((256,256))
            baseline_ASIC = int(0)
            ubaselineASIC = int(0)

        if baselineASICbeforefit != None:
            Thrs_ASIC = [thr - baseline_ASIC for thr in Thrs]
            p0_used = self.p0_baseline
        else:
            Thrs_ASIC = [thr - 0 for thr in Thrs]
            p0_used = self.p0

        print('p0:', p0_used, 'Thrs_ASIC range',min(Thrs_ASIC), max(Thrs_ASIC), 'baselineASICbeforefit:', baselineASICbeforefit, 's2_opt:', s2_opt, 'FPTfolder', FPTfolder)
        
        #Performing fit for flux on ASIC, A: recheck, this being the fit on the whole pixel grid.

        F_ASICgood=F_ASICgood[nan_cut:] #flatten it so that in curve_fit doesnt raise error
        unc_F_ASICgood=unc_F_ASICgood[nan_cut:]
        Thrs_ASIC = Thrs_ASIC[nan_cut:]

        #Cut the hanging zeros on the left and right
        trimcut=np.nonzero(F_ASICgood)  #A: recheck nonzero?? returns an array/matrix locating the indices of nonzero values
        thrs=Thrs_ASIC[trimcut[0][0]:trimcut[0][-1]+1] #A: recheck, using the locastion of non zeros values to get the Thrs at which we have an nonzero values
        F_ASICgood=np.trim_zeros(F_ASICgood) #A: recheck, pretty much doing the same thing but this time just trimming the zzeros not needed to use the location as above
        unc_F_ASICgood=unc_F_ASICgood[trimcut[0][0]:trimcut[0][-1]+1]

        print('trimcut', trimcut)
        
        F_ASICgood=F_ASICgood.flatten() #flatten it so that in curve_fit doesnt raise error
        unc_F_ASICgood=unc_F_ASICgood.flatten()

        
        #print('unc_F_ASICgood: {}'.format(unc_F_ASICgood))    
        #print('F_ASICgood: {}'.format(F_ASICgood))

        # Assuming unc_fluxASIC is a NumPy array or a list
        if all(value == 0 for value in unc_F_ASICgood):
            print("unc_F_ASICgood contains only zeros, fitting without considering unc_F_ASICgood")
            unc_F_ASICgood = None
        else:
            print("unc_F_ASICgood contains non-zero values")
        
        
        signature = inspect.signature(self.fit_function)
        params_l = list(signature.parameters)[1:] 
        n_parameters = len(params_l)

     
        poptASICgood, pcovASICgood=curve_fit(self.fit_function, thrs, F_ASICgood, sigma=unc_F_ASICgood, p0=p0_used)
        paramsASICgood = params_l
        

        #poptASICgood[1]=poptASICgood[1]/65352 #To rescale to the flux on a single pixel. Used to calculate X2ASIC  #A: recheck X2 refers to chi sqaured test so not to be confused
        poptuncASICgood=np.sqrt(np.diag(pcovASICgood))
        E0ASIC=poptASICgood[-1] #Renaming for reading ease
        uE0ASIC = poptuncASICgood[-1]

        if plot_output != None:
            return thrs, F_ASICgood, unc_F_ASICgood, poptASICgood, poptuncASICgood
        
        
        if baselineASICbeforefit != None: #if basleine is considere before fitting then this says its E0 but it is really target
            targetASIC = E0ASIC
            utargetASIC = uE0ASIC
        elif FPTfolder == '/data/bfys/apuicerc/N037_new/resultsultimate_newbaseline/' or FPTfolder == '/data/bfys/apuicerc/N037_new/resultscalibration_newbaseline/' :
            targetASIC = E0ASIC
            utargetASIC = uE0ASIC
        else: #we still want to implement the baselineASIC to find targetASIC, AFTER the computation of E0ASIC
    
            targetASIC = E0ASIC - baseline_ASIC
            utargetASIC = np.sqrt(uE0ASIC**2 + np.mean(ubaselineASIC[ubaselineASIC>=0])**2)
        

        print(paramsASICgood)
        print(poptASICgood)
        print(poptuncASICgood)
        params_header = 'poptASICgood, Value ,Uncertainty\n'
        params_data = ''
        for i in range(len(poptASICgood)):
            params_data += '{},{},{}\n'.format(paramsASICgood[i], poptASICgood[i], poptuncASICgood[i])
        
        #Adding a row for targetASICgood, utargetASICgood
        params_data += 'target,{},{}\n'.format(targetASIC, utargetASIC)
        with open(os.path.join(FPTfolder, 'ASICgood_params'+str(label)+'.csv'), 'w') as f:
            f.write(params_header)
            f.write(params_data)
        print('ASICgood params saved in {}'.format(os.path.join(FPTfolder, 'ASICgood_params'+str(label)+'.csv')))
       
            
        return targetASIC, utargetASIC,paramsASICgood, poptASICgood, poptuncASICgood




       
    
    def systematic_uncertainty_finder(self, results_selected, listofKmean = None, listofuKmean = None, label = None, reference_dataset = None):
        #with all 4 sets considerd
        #found_sys_unc_warm_cold = [0.24404864851590702, 0.07722621934704439]
        #with 3 for cold since one of the datasets is the same 
        #found_sys_unc_warm_cold = [0.24404864851590702, 0.11578490052660673] #usign the mean of all the systematic uncertainties
        found_sys_unc_warm_cold = [0.4377437173241597, 0.16414074861478997] #using the quadrature propagation
        
        if len(results_selected) == 1:

            return found_sys_unc_warm_cold
        else:
            if label == 'warm':
                return found_sys_unc_warm_cold[0]
            elif label == 'cold':
                return found_sys_unc_warm_cold[1]



        if type(listofKmean) == dict:
            listofKmean_varioustemp = []
            print('dict')
            print(listofKmean)
            for item in listofKmean:
                print(item)
                print(listofKmean[item])
                for mean in listofKmean[item]:
                    listofKmean_varioustemp.append(mean)
            print('listofKmean_varioustemp')
            print(listofKmean_varioustemp)
        
        print('listofKmean', listofKmean)
        
        if reference_dataset == None:
            reference_dataset = 0
            print('Calculating systematic uncertainty with reference dataset: {}'.format(results_selected[reference_dataset]))
            Kmean_ref = listofKmean[0]
            uKmean_ref = listofKmean[0]
        else:
            print('Calculating systematic uncertainty with reference dataset: {}'.format(results_selected[reference_dataset]))

        differences_list = []
        for i in range(len(listofKmean)):
            if reference_dataset == i:
                continue
            diff = abs(Kmean_ref - listofKmean[i] )
            differences_list.append(diff)
        
        print('Systematic uncertainties list:', differences_list)
        mean_diff = np.mean(differences_list)
        print('Mean difference {}'.format(mean_diff))
        
        final_sys_unc = np.sqrt(np.sum(x**2 for x in differences_list))
        print('Final systematic uncertainty: ', final_sys_unc)    
        return final_sys_unc

    def log_file(self, FPTfolder, label, goodpixels,Kmean, Kstd ,uKmean, paramsASIC, poptASIC, poptuncASIC, paramsASICgood, poptASICgood, poptuncASICgood):
        results_name = FPTfolder.split('/')[-2]
        log_filename = os.path.join(FPTfolder, 'log_file_' + str(results_name) + '_' + label + '.csv')

        print('Generating log file for {} ...'.format(results_name))
        total_log_string = ''
        pixels = goodpixels
        percentage_pixels = pixels / (256 ** 2) * 100

        # Prepare the header and initial data
        total_log_string += 'Fittype, Pixels ,Percentage\n'
        total_log_string += 'good,{},{}\n'.format(pixels, percentage_pixels)

        total_log_string += 'Kmean, Kstd ,uKmean\n'
        total_log_string += '{},{},{}\n'.format(Kmean, Kstd ,uKmean)

        # Additional data for paramsASICgood
        total_log_string += 'poptASICgood, Value ,Uncertainty\n'
        for i in range(len(poptASICgood)):
            total_log_string += '{},{},{}\n'.format(paramsASICgood[i], poptASICgood[i], poptuncASICgood[i])

        # Additional data for paramsASIC
        total_log_string += 'poptASIC, Value ,Uncertainty\n'
        for i in range(len(poptASIC)):
            total_log_string +='{},{},{}\n'.format(paramsASIC[i], poptASIC[i], poptuncASIC[i])

        # Read the existing log file if it exists
        if os.path.exists(log_filename):
            with open(log_filename, 'r') as f:
                existing_data = f.read()
        else:
            existing_data = ''


        # Only write to the file if there's a difference
        if total_log_string != existing_data:
            with open(log_filename, 'w') as f:
                f.write(total_log_string)
            print('Log file updated and saved in {}'.format(log_filename))
        else:
            print('No changes detected. Log file not updated.')

