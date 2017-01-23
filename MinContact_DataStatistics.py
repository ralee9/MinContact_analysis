###############################################################################
##    Randy Lee
##    VIA Lab
##    Department of Bioengineering
##    University of Pittsburgh
##
##    Minimum Contact Experiment Analysis
##
##    Function: Crop individual trial data to take statistics
##      of force traces
##
##  Analysis procedure:
##    (1) MinContact_RawDataProcessing to splitt raw data and graph time
##        series and Fourier/Power spectra
##    (2) MinContact_DataStatistics to obtain mean force and std for each
##        experimental condition, and graph dependent variable analysis of
##        spectral power
##
##
##
##     Last updated : 20 October 2016
###############################################################################

import os
import numpy
import math
import time
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


subjectNumber = 'S2008'

#rawDataDirectory = 'C:\\Users\\Randy Lee\\Documents\\VIA Lab\\HHFM Data\\MinimumContact_rawData\\'
rawDataDirectory = 'D:\\Randy Lee\\Documents\\VIA Lab\\HHFM Data\\MinimumContact_rawData'

###############################################################################
##
##    Local script functions
##
###############################################################################

def raw_stats(fullDataArray, trialParameters, trialCounter):
    ''' Calculate statistics of raw force data from GS0-100

    '''
    
    print "Calculating raw data stats for Trial #{}".format(trialCounter+1)
    ## store raw data
    # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
    # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
    # [8] bit2; [9] currentlyInTrial voltage;
    # [10] GS0Force; [11] filteredForce_201; [12] filteredForce_501
    GS0_ForceGrams = fullDataArray[:,10]

    # End of files saved as 2500 ticks AFTER currentlyInTrial pin change
    # visFD turns off 6s before end of trial
    # trim force to noFD section, ignore first second
    endStamp = len(fullDataArray) - 2500
    startStamp = endStamp - 6000
    force_noFD = GS0_ForceGrams[startStamp + 1000 : endStamp]

    ## [0] trial_ID [1] direction [2] mag [3] perceived_target
    direction = trialParameters[trialCounter, 1]    
    magnification = trialParameters[trialCounter, 2]
    perceivedForce = trialParameters[trialCounter, 3]    

    # remove sections of force where HHFM is not in contact with sensor
    force_inContact = exclude_noncontact(force_noFD, trialCounter)

    # save stats
    statistics = numpy.zeros(shape = (7))
    statistics[0] = direction
    statistics[1] = magnification
    statistics[2] = perceivedForce
    statistics[3] = numpy.mean(force_noFD)
    statistics[4] = numpy.std(force_noFD)
    statistics[5] = numpy.mean(force_inContact)
    statistics[6] = numpy.std(force_inContact)

    return statistics

def filter_stats(fullDataArray, trialParameters, trialCounter):
    ''' Calculate statistics of Gaussian filtered force data from GS0-100

    '''
    
    print "Calculating filtered data stats for Trial #{}".format(trialCounter+1)

    ## store raw data
    # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
    # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
    # [8] bit2; [9] currentlyInTrial voltage;
    # [10] GS0Force; [11] filteredForce_201; [12] filteredForce_501
    filteredForce = fullDataArray[:,11]

    # Define start and stop points for trim, extract trial parameters
    # visFD turns off 6s before end of trial
    # trim filtered force to noFD section, ignore first second
    endStamp = len(filteredForce) - 2500
    startStamp = endStamp - 6000
    filteredForce_noFD = filteredForce[startStamp + 1000 : endStamp]

    ## [0] trial_ID [1] direction [2] mag [3] perceived_target
    direction = trialParameters[trialCounter, 1]    
    magnification = trialParameters[trialCounter, 2]
    perceivedForce = trialParameters[trialCounter, 3]     
    
    # save stats
    statistics = numpy.zeros(shape = (5))
    statistics[0] = direction
    statistics[1] = magnification
    statistics[2] = perceivedForce
    statistics[3] = numpy.mean(filteredForce_noFD)
    statistics[4] = numpy.std(filteredForce_noFD)

    return statistics

def linear_fit(fullDataArray, trialParameters, trialCounter):
    ''' Calculate linear regression for entire noFD section of each trial

    '''
    
    print "Calculating linear regression for noFD section of Trial #{}".format(trialCounter+1)
    ## store raw data
    # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
    # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
    # [8] bit2; [9] currentlyInTrial voltage;
    # [10] GS0Force; [11] filteredForce_201; [12] filteredForce_501
    GS0_ForceGrams = fullDataArray[:,10]

    # Define start and stop points for trim, extract trial parameters
    # visFD turns off 6s before end of trial
    # slice trim data to noFD section, ignore first second
    endStamp = len(GS0_ForceGrams) - 2500
    startStamp = endStamp - 6000
    rawForce_noFD = GS0_ForceGrams[startStamp + 1000 : endStamp]

    ## [0] trial_ID [1] direction [2] mag [3] perceived_target
    direction = trialParameters[trialCounter, 1]    
    magnification = trialParameters[trialCounter, 2]
    perceivedForce = trialParameters[trialCounter, 3]  

    # perform linear regression
    x = numpy.arange(len(rawForce_noFD))

    slope_noFD, intercept_noFD, r_noFD, p_val, std_err = stats.linregress(x, 
        rawForce_noFD)

    r2_noFD = r_noFD**2

    # save stats
    statistics = numpy.zeros(shape=(6))
    statistics[0] = direction
    statistics[1] = magnification
    statistics[2] = perceivedForce
    statistics[3] = slope_noFD
    statistics[4] = intercept_noFD
    statistics[5] = r2_noFD

    return statistics

def window_linear_fit(fullDataArray, trialParameters, trialCounter):
    ''' Calculate linear regression for noFD within small overlapping windows
    
    '''
    
    print "Calculating shifting linear regression for Trial #{}".format(trialCounter+1)

    ## store raw data
    # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
    # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
    # [8] bit2; [9] currentlyInTrial voltage;
    # [10] GS0Force; [11] filteredForce_201; [12] filteredForce_501
    GS0_ForceGrams = fullDataArray[:,10]
        
    # Define start and stop points for trim, extract trial parameters
    # visFD turns off 6s before end of trial
    # slice trim data to noFD section, ignore first second
    endStamp = len(GS0_ForceGrams) - 2500
    startStamp = endStamp - 6000    
    rawForce_noFD = GS0_ForceGrams[startStamp + 1000 : endStamp]
    
    # extract trial parameters
    ## [0] trial_ID [1] direction [2] mag [3] perceived_target
    direction = trialParameters[trialCounter, 1]    
    magnification = trialParameters[trialCounter, 2]
    perceivedForce = trialParameters[trialCounter, 3]  
     
    # windowed analysis for noFD section
    # set up window size and shift parameters
    windowSize = 150
    shift = 50
    windowStart = numpy.arange(0, len(rawForce_noFD), shift)

    runningSum_noFD = 0
    reversal_noFD= 0
    last = 0
    results_noFD = numpy.zeros(shape = (len(windowStart), 4))
    for i in range(len(windowStart)):
        if windowStart[i] + windowSize > len(rawForce_noFD):
            # if current section exceeds length of array, slice only until end
            window_force = rawForce_noFD[windowStart[i] : len(rawForce_noFD)]
        else:
            window_force = rawForce_noFD[windowStart[i] : windowStart[i] + windowSize]
        
        x = numpy.arange(len(window_force))
        slope, intercept, r, p_val, std_err = stats.linregress(x, window_force)
        
        # [0] window #; [1] start stamp [2] slope [3] r^2 
        results_noFD[i, :] = [i, windowStart[i], slope, r**2]
        if slope >= 0:
            runningSum_noFD += 1
            if last == -1:
                reversal_noFD += 1
            last = 1
        elif slope < 0:
            runningSum_noFD -= 1
            if last == 1:
                reversal_noFD += 1
            last == -1
    
    # save and store regression results
    savePath = os.getcwd() + '\\win_regression_Trial{0:02d}.txt'.format(trialCounter+1)
    
    numpy.savetxt(savePath, results_noFD, fmt = '%.4f', delimiter = '\t',
                  header = 'window number\t start stamp\t slope\t R^2')
    
    # save stats 
    statistics = numpy.zeros(shape = (5))
    statistics[0] = direction
    statistics[1] = magnification
    statistics[2] = perceivedForce
    statistics[3] = runningSum_noFD
    statistics[4] = reversal_noFD
    
    return statistics

def power_stats(fourierCoeff, fullDataArray, trialParameters, trialCounter):
    '''Calculate and save power band spectra and cumulative power for each trial

    '''
    
    print "Calculating power statistics for Trial #{}".format(trialCounter+1)

    # extract raw data
    fftRawNoFD = fourierCoeff[:,1]
    GS0_ForceGrams = fullDataArray[:,10]

    # extract trial parameters
    # [0] trial_ID [1] direction [2] mag [3] perceived_target
    direction = trialParameters[trialCounter, 1]    
    magnification = trialParameters[trialCounter, 2]
    perceivedForce = trialParameters[trialCounter, 3]  

    power_raw_1_4 = numpy.zeros(shape = (16))
    index = 0
    # calculate power spectra coefficients for freq 1 - 4 Hz as square of
    # fourier coeff
    for i in range(5, 21):
        power_raw_1_4[index] = fftRawNoFD[i] ** 2.
        index += 1

    cumPower_raw_1_4 = numpy.sum(power_raw_1_4)

    # calculate power spectra coefficients for freq 4 - 7 Hz as square of
    # fourier coeff
    power_raw_4_7 = numpy.zeros(shape = (16))
    index = 0
    for i in range(20, 36):
        power_raw_4_7[index] = fftRawNoFD[i] ** 2.
        index += 1

    cumPower_raw_4_7 = numpy.sum(power_raw_4_7)

    # calculate power spectra coefficients for freq 7 - 10 Hz as square of
    # fourier coeff
    power_raw_7_10 = numpy.zeros(shape = (16))
    index = 0
    for i in range(35, 51):
        power_raw_7_10[index] = fftRawNoFD[i] ** 2.
        index += 1

    cumPower_raw_7_10 = numpy.sum(power_raw_7_10)


    ## calculate FFT/power with non-contact sections excluded
    endStamp = len(fullDataArray) - 2500
    startStamp = endStamp - 6000
    force_noFD = GS0_ForceGrams[startStamp + 1000 : endStamp]
    
    # remove sections of force where HHFM is not in contact with sensor
    force_inContact = exclude_noncontact(force_noFD, trialCounter)
    
    if force_inContact.size == 0:
        # in case entire noFD sections 
        fft_inContact = numpy.zeros(shape=60)
    else:
        fft_inContact = numpy.abs(numpy.fft.fft(force_inContact, axis = 0))
    
    #timestep= 0.001
    #freq_inContact = numpy.fft.fftfreq(fft_inContact.size, d = timestep)
    
    # calculate cumulative power by square of FFT coeff
    index = 0
    power_contact_1_4 = numpy.zeros(shape=(16))
    # power over 1-4 Hz band
    for i in range(5,21):
        power_contact_1_4[index] = fft_inContact[i] ** 2.
        index += 1
    
    cumPower_contact_1_4 = numpy.sum(power_contact_1_4)
    
    # power over 4-7Hz band
    index = 0
    power_contact_4_7 = numpy.zeros(shape=(16))
    for i in range(20, 36):
        power_contact_4_7[index] = fft_inContact[i] ** 2.
        index += 1
     
    cumPower_contact_4_7 = numpy.sum(power_contact_4_7)
    
    # power over 7-10Hz band
    index = 0
    power_contact_7_10 = numpy.zeros(shape=(16))
    for i in range(35, 51):
        power_contact_7_10[index] = fft_inContact[i] ** 2.
        index += 1
        
    cumPower_contact_7_10 = numpy.sum(power_contact_7_10)

    # save stats
    statistics = numpy.zeros(shape = (9))
    statistics[0] = direction
    statistics[1] = magnification
    statistics[2] = perceivedForce
    statistics[3] = cumPower_raw_1_4
    statistics[4] = cumPower_raw_4_7
    statistics[5] = cumPower_raw_7_10
    statistics[6] = cumPower_contact_1_4
    statistics[7] = cumPower_contact_4_7
    statistics[8] = cumPower_contact_7_10

    return statistics

def contact_stats(fullDataArray, trialParameters, trialCounter):
    ''' Recover the time spent OUT of contact w/ sensor using Welch's t-test

    Checks first 250 ms of force data against 30ms windows to see if two groups
    are significantly different using Welch's t-test. If p > 0.05, the null 
    hypothesis that the two samples are equal is NOT rejected, so the HHFM is
    out of contact with sensor.
    '''
    
    print "Calculating contact statistics for Trial #{}".format(trialCounter+1)

    ## recover raw force data
    # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
    # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
    # [8] bit2; [9] currentlyInTrial voltage;
    # [10] GS0Force; [11] filteredForce_201; [12] filteredForce_501
    GS0_ForceGrams = fullDataArray[:,10]

    ## [0] trial_ID [1] direction [2] mag [3] perceived_target
    direction = trialParameters[trialCounter, 1]    
    magnification = trialParameters[trialCounter, 2]
    perceivedForce = trialParameters[trialCounter, 3]  

    # trim force trace to noFD section only, ignoring first second
    endStamp = len(GS0_ForceGrams) - 2500
    startStamp = endStamp - 6000
    GS0Force_noFD = GS0_ForceGrams[startStamp + 1000 : endStamp]

    # start of trial, first 250 ms guarenteed no contact by HHFM
    trialStart_raw = GS0_ForceGrams[0:250]
    
    # check raw data in 30ms windows, skip first second of noFD section
    # 30ms b/c it takes 30 samples to produce approx normal distribution
    windowSize = 30
    windowStart_noFD = numpy.arange(0, len(GS0_ForceGrams), windowSize)
    noContact_timestamps = numpy.zeros(shape = windowStart_noFD.shape)
    
    ## calculate and save p values for noFD section only
    sections_noContact = 0
    for i in range(len(windowStart_noFD)):
        # extract 30ms windows from raw force data
        if windowStart_noFD[i] + windowSize >= len(GS0Force_noFD):
            window_force = GS0Force_noFD[windowStart_noFD[i] : len(GS0Force_noFD)]
        else:
            window_force = GS0Force_noFD[windowStart_noFD[i] : windowStart_noFD[i] + windowSize]
        
        # detect contact using Welch's t-test, since n1 != n2
        curr_t, curr_p = stats.ttest_ind(trialStart_raw, window_force,
                                         equal_var = False)
            
        # if p is large, we cannot reject the null hypothesis that the two
        # samples are from the same distribution, so the subject is NOT in 
        # contact with the sensor
        # p = 0.025 chosen empirically as sufficiently sensitive
        if curr_p > 0.025:
            print "Non-contact detected: win_start = {}, curr_p = {}".format(windowStart_noFD[i], curr_p)
            noContact_timestamps[sections_noContact] = windowStart_noFD[i]
            sections_noContact += 1
    
    ## calculate and save p values for entire trial, ignore first 500ms b/c
    # trial data saved as 500ms before trial start
    windowStart_all = numpy.arange(500, len(GS0_ForceGrams), windowSize)
    trial_pVal = numpy.zeros(shape=windowStart_all.shape)
    for i in range(len(windowStart_all)):
        if windowStart_all[i] + windowSize >= len(GS0_ForceGrams):
            window_force = GS0_ForceGrams[windowStart_all[i] : len(GS0_ForceGrams)]
        else:
            window_force = GS0_ForceGrams[windowStart_all[i] : windowStart_all[i] + windowSize]
        
        curr_t, curr_p = stats.ttest_ind(trialStart_raw, window_force,
                                         equal_var = False)
        
        trial_pVal[i] = curr_p

    print "** num 30ms sections_noContact = {} **".format(sections_noContact) 

    if sections_noContact == 0:
        # if always in contact, save timestamp as -1
        noContact_timestamps[0] = -1

    # trim zeros and save timestamps for no contact
    trial_pVal = numpy.column_stack((windowStart_all, trial_pVal))
    noContact_timestamps = numpy.trim_zeros(noContact_timestamps, trim = 'b')
    fileName = os.getcwd() + '\\noContact_timestamps_Trial{0:02d}.txt'.format(trialCounter+1)
    numpy.savetxt(fileName, noContact_timestamps, 
                  header = 'non-contact timestamps')
                  
    # save p values for entire trial for later graphing
    fileName = os.getcwd() + '\\noContact_pVal_Trial{0:02d}.txt'.format(trialCounter+1)
    numpy.savetxt(fileName, trial_pVal,
                  header = 'windowStart\t p_Value', delimiter = '\t')
    
    # convert # of no contact sections to seconds
    noContact_time_sec = (sections_noContact * windowSize)/1000.
    contact_time_sec = 5. - noContact_time_sec
    if contact_time_sec < 0:
        # sometimes 30ms sections give 5.01 sec period out of contact
        contact_time_sec = 0
    print "-->time out of contact = {} seconds".format(noContact_time_sec)
    print "-->time IN contact = {} seconds\n".format(contact_time_sec)

    # save stats
    statistics = numpy.zeros(shape = (6))
    statistics[0] = direction
    statistics[1] = magnification
    statistics[2] = perceivedForce
    statistics[3] = sections_noContact
    statistics[4] = noContact_time_sec
    statistics[5] = contact_time_sec

    return statistics

def graph_window_fits(trialParameters, blockNumberString):
    ''' Graph windowed linear regression fits with visFD/noFD and mag on/off together
    
    '''
    
    # open PDF for plots
    pdfFileName = os.getcwd() + '\\' + blockNumberString + '_windowedFitGraphs.pdf'
    pp = PdfPages(pdfFileName)
    
    # array with ordered list of file names
    blockDataList = numpy.empty(6, dtype = 'str')
    blockDataList = reorder_trials(blockNumberString)
    
    trialCounter = 0
    figCounter = 1
    doneFlag = 0
    print "Graphing windowed linear fits..."
    while(doneFlag == 0):
        
        # extract unmag
        file_name = os.getcwd() + '\\win_regression_{}'.format(blockDataList[trialCounter])
        fullData_nomag = numpy.genfromtxt(file_name, skip_header = 1)
        
        # extract mag=2 data
        file_name = os.getcwd() + '\\win_regression_{}'.format(blockDataList[trialCounter + 1])
        fullData_mag1 = numpy.genfromtxt(file_name, skip_header = 1)
        
        # extract mag=4 data
        file_name = os.getcwd() + '\\win_regression_{}'.format(blockDataList[trialCounter + 2])
        fullData_mag2 = numpy.genfromtxt(file_name, skip_header = 1)
        
        # create figure and begin plotting        
        fig = plt.figure(figCounter)
        print "Figure #{}...\n".format(figCounter)
        figureTitle = 'Block_Half#{}; (Direction:{} Target Force = {})'.format(
                        figCounter, int(trialParameters[trialCounter][1]),
                        5.0)
        fig.suptitle(figureTitle, fontsize = 12, fontweight = 'bold')
                                
        ax1 = plt.subplot(311)
        subplotTitle = 'Magnification = 0x'
        ax1.set_title(subplotTitle, fontsize = 8)
        
        x1 = fullData_nomag[:,1]
        y1 = fullData_nomag[:,2]
        nomag, = plt.plot(x1,y1,label='Magnification = 0x')
        
        plt.setp(ax1.get_xticklabels(), visible = False)
        ax1.set_ylabel('Slope')
        ax1.set_xlim((0,6500))
        ax1.set_ylim((-0.15, 0.15))
        ax1.axhline(y=0, color = 'k', lw = 1.00)
        ax1.grid(True, which = 'both')
        plt.legend(loc = 'lower left', fontsize = 8)

        ax2 = plt.subplot(312)
        subplotTitle = 'Magnification = 5x'
        ax2.set_title(subplotTitle, fontsize = 8)
        
        x2 = fullData_mag1[:,1]
        y2 = fullData_mag1[:,2]
        mag1, = plt.plot(x2,y2,label='Magnification = 5x')
        
        plt.legend(loc = 'lower left', fontsize = 8)
        ax2.set_ylabel('Slope')
        plt.setp(ax2.get_xticklabels(), visible = False)
        #ax2.set_xlabel('Window Start')
        ax2.set_xlim((0, 6500))
        ax2.set_ylim((-0.15, 0.15))
        ax2.axhline(y = 0, color = 'k', lw = 1.00)
        ax2.grid(True, which = 'both')

        ax3 = plt.subplot(313)
        subplotTitle = 'Magnification = 10x'
        ax3.set_title(subplotTitle, fontsize = 8)
        x3 = fullData_mag2[:,1]
        y3 = fullData_mag2[:,2]
        mag2, = plt.plot(x3,y3,label='Magnification = 10x')
        
        plt.legend(loc = 'lower left', fontsize = 8)
        ax3.set_ylabel('Slope')
        ax3.set_xlabel('Window Start')
        ax3.set_xlim((0, 6500))
        ax3.set_ylim((-0.15, 0.15))
        ax3.axhline(y = 0, color = 'k', lw = 1.00)
        ax3.grid(True, which = 'both')
        
        pp.savefig(figCounter)
        plt.close('all')
        
        # advance trial/fig counters
        trialCounter += 3
        figCounter += 1
        
        if trialCounter >= 6:
            plt.close('all')
            pp.close()
            print "\nPDF closed.\n"
            print "{} slope data graphing completed!\n".format(blockNumberString)
            doneFlag = 1
    
    return

def graph_contact(trialParameters, blockNumberString):
    ''' Graph trial data with no-contact sections highlighted.

    '''
    
    # open PDF for plots
    pdfFileName = os.getcwd() + '\\' + blockNumberString + '_noContact_DataPlots.pdf'
    pp = PdfPages(pdfFileName)
    
    # array with ordered list of file names
    blockDataList = numpy.empty(6, dtype = 'str')
    blockDataList = reorder_trials(blockNumberString)

    doneFlag = 0
    trialCounter = 0
    figCounter = 1
    while(doneFlag == 0):
        # Print trial number
        print "Graphing non-contact Trial #{} data...".format(trialCounter+1)

        # Load data
        currentFileName = blockDataList[trialCounter]
        currentFilePath = os.getcwd() + '\\' + currentFileName

        currentData = numpy.genfromtxt(currentFilePath, skip_header=1)
        currentDataSize = currentData.shape

        GS0_ForceGrams = currentData[:,10]

        ## [0] trial_ID [1] direction [2] mag [3] perceived_target
        direction = trialParameters[trialCounter, 1]    
        magnification = trialParameters[trialCounter, 2]
        perceivedForce = trialParameters[trialCounter, 3]

        targetForce = direction * 15

        fig = plt.figure(figCounter)
        print "Figure #{}...\n".format(figCounter)
        figureTitle = 'Trial#{}; Direction = {}; Magnification Gain = {} '\
            'Perceived Target = {})'.format(figCounter, direction, magnification,
                perceivedForce)
        fig.suptitle(figureTitle, fontsize=10, fontweight='bold')

        ax1 = plt.subplot(211)
        subplotTitle = 'GS0-100 Force Sensor'
        ax1.set_title(subplotTitle, fontsize=10)

        ## Plot trial force to top subplot
        x = numpy.arange(currentDataSize[0])

        # set axes and labels
        plt.axis(ymin=-45, ymax=45)
        plt.xticks(numpy.arange(0,len(x),1000),
                   numpy.arange(0,math.ceil(len(x)/1000)+1,1, dtype=numpy.int))
        plt.setp(ax1.get_xticklabels(), visible=False)
        plt.ylabel('Applied Force '+r'$(grams)$', fontsize=10)
        #plt.xlabel('Time (seconds)')
        plt.grid(axis='y', which='major')

        # Plot target force and x axis
        plt.axhline(y=targetForce, color='grey', ls='--', lw=1.25)
        plt.axhline(y=0, color='k', lw=1.0)

        # Plot raw force
        GS0_raw, = plt.plot(x, GS0_ForceGrams, 'k', lw=1)
        ax1_xmin, ax1_xmax = ax1.get_xlim()

        # Plot vertical lines for important time points
        endStamp = currentDataSize[0] - 2500
        startStamp = 500
        plt.axvline(x=startStamp, color='g', ls='--', lw=1.25)
        plt.axvline(x=endStamp, color='r', ls='--', lw=1.25)

        # visFD ends 6 seconds before trial end
        plt.axvline(x=endStamp - 6000, color='r', ls='--', lw=1.25)
        plt.axvline(x=endStamp - 5000, color='grey', ls='--', lw=0.75)
        #plt.text(x=endStamp - 4000, y=35, s='<---(no visual feedback)--->', 
        #         style='italic', fontsize='xx-small')
        
        # Load timestamps, graph in pairs
        fileName = os.getcwd() + '\\noContact_timestamps_Trial{0:02d}.txt'.format(trialCounter+1)
        timestamps = numpy.genfromtxt(fileName, skip_header=1)
        y = numpy.arange(-30,30)
        if timestamps.size == 1:
            if timestamps == -1:
                # check first timestamp in case always in contact
                # if so, break out of loop immediately
                print "timestamp = -1, pass"
                pass
            else:
                left = (endStamp - 5000) + timestamps
                right = left+30
                plt.fill_betweenx(y, x1=left, x2=right, facecolor='y', alpha=0.2)
        else:
            for i in range(0,len(timestamps)):
                if timestamps[i] == -1:
                    # check first timestamp in case always in contact
                    # if so, break out of loop immediately
                    print "break"
                    break
                # time stamps are relative to 1 second past the noFD section
                left = (endStamp - 5000) + timestamps[i]
                right = left+30
                plt.fill_betweenx(y, x1=left, x2=right, facecolor='y', alpha=0.2)

        ax2 = plt.subplot(212)
        subplotTitle = "Welch's t-test p-values"
        ax2.set_title(subplotTitle, fontsize=10)
        
        # set axes and labels
        ax2.set_xlim(left=ax1_xmin, right=ax1_xmax)
        plt.axis(ymin=-0.05, ymax=1)
        plt.xticks(numpy.arange(0,len(x),1000),
                   numpy.arange(0,math.ceil(len(x)/1000)+1,1, dtype=numpy.int))
        plt.xlabel('Time '+r'$(sec)$', fontsize=10)
        plt.ylabel('p-value', fontsize=10)
        plt.grid(axis='y', which='major')
        
        # Load and extract data
        fileName = os.getcwd() + '\\noContact_pVal_Trial{0:02d}.txt'.format(trialCounter+1)
        p_val_data = numpy.genfromtxt(fileName, skip_header=1)

        windowStart = p_val_data[:,0]
        p_values = p_val_data[:,1]
        
        plt.plot(windowStart, p_values, 'k')
        
        # Plot vertical lines for important time points
        endStamp = currentDataSize[0] - 2500
        startStamp = 500
        plt.axvline(x=startStamp, color='g', ls='--', lw=1.25)
        plt.axvline(x=endStamp, color='r', ls='--', lw=1.25)

        # visFD ends 6 seconds before trial end
        plt.axvline(x=endStamp - 6000, color='r', ls='--', lw=1.25)
        plt.axvline(x=endStamp - 5000, color='grey', ls='--', lw=0.75)
        
        plt.axhline(y=0, color='k', lw=1.00)
        plt.axhline(y=0.025, color='grey', ls='--', lw=1.25)        

        # save figure to PDF
        pp.savefig(figCounter)
        plt.close('all')

        trialCounter += 1
        figCounter += 1

        if trialCounter >= 6:
            plt.close('all')
            pp.close()
            print "PDF closed."
            print "{} contact data graphing completed!\n".format(blockNumberString)
            doneFlag = 1

    return

def exclude_noncontact(rawForceArray, trialCounter):
    ''' Removes sections of GS0 force sensor data when HHFM is not in contact
    
    '''

    fileName = os.getcwd() + \
               '\\noContact_timestamps_Trial{0:02d}.txt'.format(trialCounter+1)
    
    timestamps = numpy.genfromtxt(fileName, skip_header=1)
    to_del = numpy.zeros(shape=1)
    
    if timestamps.size == 1:
        left = timestamps
        if left == -1: 
            print "timestamp = -1, pass"
            pass
        else:
            right = left + 30
            to_del = numpy.append(to_del, numpy.arange(left,right+1))
    else:
        for i in range(len(timestamps)):
            left = timestamps[i]
            if left == -1:
                # if always in contact, no data to delete, so escape out of for loop
                break
            right = left + 30
            to_del = numpy.append(to_del, numpy.arange(left,right+1))

    numpy.trim_zeros(to_del, trim='fb')
    if to_del.size == 0:
        # check to see if to_del list is empty. if so, copy force trace
        force_inContact = rawForceArray
    else:
        force_inContact = numpy.delete(rawForceArray, to_del)   
    
    return force_inContact

def reorder_trials(blockNumberString):
    ''' Obtain ordered list of trial data by direction then magnification gain.
    
    From the known random ordering of each block, return ordered trial list.    
    '''
    ## List of prototype numbers, ordered by direction, then magnification gain
    # [1, 2, 3, 4, 5, 6]

    ## blockNumberString passed as S#_BlockX

    if blockNumberString[-6:] == 'Block1':
        orderedTrialList = numpy.array(['Trial03.txt', 'Trial02.txt', 
                                        'Trial04.txt', 'Trial06.txt', 
                                        'Trial01.txt', 'Trial05.txt'])

    elif blockNumberString[-6:] == 'Block2':
        orderedTrialList = numpy.array(['Trial06.txt', 'Trial03.txt',
                                        'Trial04.txt', 'Trial01.txt',
                                        'Trial05.txt', 'Trial02.txt'])

    elif blockNumberString[-6:] == 'Block3':
        orderedTrialList = numpy.array(['Trial06.txt', 'Trial03.txt',
                                        'Trial02.txt', 'Trial05.txt',
                                        'Trial04.txt', 'Trial01.txt'])

    elif blockNumberString[-6:] == 'Block4':
        orderedTrialList = numpy.array(['Trial03.txt', 'Trial01.txt',
                                        'Trial05.txt', 'Trial04.txt',
                                        'Trial06.txt', 'Trial02.txt'])

    elif blockNumberString[-6:] == 'Block5':
        orderedTrialList = numpy.array(['Trial01.txt', 'Trial02.txt',
                                        'Trial05.txt', 'Trial06.txt',
                                        'Trial03.txt', 'Trial04.txt'])

    return orderedTrialList


###############################################################################
##
##    MAIN: Process Data
##
###############################################################################

start = time.time()

print "\n\n*** MinContact_DataStatistics.py - START ***\n\n"

print "Loading raw data folders... \n"

os.chdir(rawDataDirectory)

print "Accessing subject #{} data... \n".format(subjectNumber)
subjectDirectory = rawDataDirectory + '\\' + subjectNumber
os.chdir(subjectNumber)

## Currently in HHFM_Data/subjectNumber/
subjectDataList = os.listdir(os.getcwd())
subjectDataListSize = len(subjectDataList)

## Loop through all five blocks
for block in range(0, 5):
    ## Enter block folder
    blockFolderName = '{}_Block{}'.format(subjectNumber, block + 1)
    blockFolder = os.getcwd() + '\\' + blockFolderName
    os.chdir(blockFolder)

    ## Currently in HHFM_Data/subjectNumber/BlockX/
    ## Load timestamps and trialParameters
    filePath = os.getcwd() + '\\' + blockFolderName + '_timestamps.txt'
    timestamps = numpy.genfromtxt(filePath)

    filePath = blockFolder + '\\' + blockFolderName + '_trialParameters.txt.'
    trialParameters = numpy.genfromtxt(filePath)

    dataStatistics = numpy.zeros(shape = (6,7))
    filterStats = numpy.zeros(shape = (6,5))
    linearFitStats = numpy.zeros(shape = (6,6))
    winLinearFitStats = numpy.zeros(shape = (6,5))
    powerStats = numpy.zeros(shape = (6,9))
    contactStats = numpy.zeros(shape = (6,6))

    blockDataList = numpy.empty(6, dtype = 'str')
    blockDataList = reorder_trials(blockFolderName)

    for trialCounter in range(0,6):
        # Choose trial data from ORDERED list of trial parameters
        fileName = blockDataList[trialCounter]
        print "Accessing {} of {}...".format(fileName, blockFolderName[-6:])
        trialPath = blockFolder + '\\' + fileName
        trialData = numpy.genfromtxt(trialPath, skip_header = 1)

        # Obtain fourier coefficients, already ordered by distal target
        fileName = blockFolderName + '_{}_FourierCoeff.txt'.format(trialCounter+1)
        trialPath = blockFolder + '\\' + fileName
        fourierCoeff = numpy.genfromtxt(trialPath, skip_header = 1)

        # Find time spent out of contact with sensor, record timestamp file
        contactStats[trialCounter] = contact_stats(trialData, trialParameters,
                                                   trialCounter)

        # Return force statistics from cropped trials, taking into account
        # non contact time
        dataStatistics[trialCounter] = raw_stats(trialData, trialParameters,
                                                 trialCounter)

        # Filtered force statistics
        filterStats[trialCounter] = filter_stats(trialData, trialParameters,
                                                 trialCounter)

        # Linear regression analysis
        linearFitStats[trialCounter] = linear_fit(trialData, trialParameters,
                                                  trialCounter)
        
        # Windowed linear regression analysis
        winLinearFitStats[trialCounter] = window_linear_fit(trialData,
                                                trialParameters, trialCounter)

        # Power spectra statistics
        powerStats[trialCounter] = power_stats(fourierCoeff, trialData,
                                               trialParameters, trialCounter)

    # graph windowed linear fit data
    graph_window_fits(trialParameters, blockFolderName)

    # graph raw data with detected non contact regions
    graph_contact(trialParameters, blockFolderName)

    if block == 0:
        # Create and fill combined stats with data statistics from block 1
        combinedStatistics = numpy.zeros(shape = (dataStatistics.shape))
        combinedStatistics = dataStatistics

        combinedFilterStatistics = numpy.zeros(shape = (filterStats.shape))
        combinedFilterStatistics = filterStats

        combinedLinearFits = numpy.zeros(shape = (linearFitStats.shape))
        combinedLinearFits = linearFitStats
        
        combinedWinLinearFitStats = numpy.zeros(shape = (winLinearFitStats.shape))
        combinedWinLinearFitStats = winLinearFitStats

        combinedPowerStatistics = numpy.zeros(shape = (powerStats.shape))
        combinedPowerStatistics = powerStats

        combinedContactStatistics = numpy.zeros(shape = (contactStats.shape))
        combinedContactStatistics = contactStats

    else:
        combinedStatistics = numpy.vstack((combinedStatistics, dataStatistics))

        combinedFilterStatistics = numpy.vstack((combinedFilterStatistics,
                                                 filterStats))

        combinedLinearFits = numpy.vstack((combinedLinearFits, linearFitStats))
        
        combinedWinLinearFitStats = numpy.vstack((combinedWinLinearFitStats,
                                                  winLinearFitStats))

        combinedPowerStatistics = numpy.vstack((combinedPowerStatistics,
                                                powerStats))

        combinedContactStatistics = numpy.vstack((combinedContactStatistics,
                                                  contactStats))

    # Block completed, return to subject directory to continue processing
    os.chdir(subjectDirectory)

## Finalize and save contact time statistics
trialColumn = numpy.array(numpy.arange(1,31))
combinedContactStatistics = numpy.column_stack((trialColumn,
                                                combinedContactStatistics))

contactFileName = subjectDirectory + '\\contactStatistics.txt'
numpy.savetxt(contactFileName, combinedContactStatistics, fmt = '%.4f',
             header = 'Trial Number\t Direction\t Magnification Gain\t'
             'Perceived Target Force\t # 30ms Sections No Contact\t'
             'Time no contact (seconds)\t Time in contact (seconds)',
             delimiter = '\t')

## Finalize and save mean/SD statistics
combinedStatistics = numpy.column_stack((trialColumn, combinedStatistics))

statisticsFileName = subjectDirectory + '\\experimentStatistics.txt'
numpy.savetxt(statisticsFileName, combinedStatistics , fmt= '%.4f',
              header = 'Trial Number\t Direction\t Magnification Gain\t'
              'Perceived Target Force\t Mean Force (raw)\t SD Force (raw)\t'
              'Mean Force (in contact)\t SD Force (in contact)',
              delimiter = '\t')


## Finalize and save filtered force statistics
combinedFilterStatistics = numpy.column_stack((trialColumn, combinedFilterStatistics))

filteredFileName = subjectDirectory + '\\filteredForceStatistics.txt'
numpy.savetxt(filteredFileName, combinedFilterStatistics , fmt= '%.4f',
              header = 'Trial Number\t Direction\t Magnification Gain\t'
              'Perceived Target Force\t Mean Filtered Force (raw)\t'
              'SD Filtered Force (raw)',
              delimiter = '\t')

## Finalize and save linear fit statistics
combinedLinearFits = numpy.column_stack((trialColumn, combinedLinearFits))

linearFitsFileName = subjectDirectory + '\\linearFitStatistics.txt'
numpy.savetxt(linearFitsFileName, combinedLinearFits, fmt = '%.6f',
              header = 'Trial Number\t Direction\t Magnification Gain\t' 
              'Perceived Target\t noFD slope\t noFD intercept\t noFD R^2',
              delimiter = '\t')

## Finalize and save windowed linear fit summary stats
combinedWinLinearFitStats = numpy.column_stack((trialColumn, combinedWinLinearFitStats))

winLinearFitsFileName = subjectDirectory + '\\windowedLinearFitStatistics.txt'
numpy.savetxt(winLinearFitsFileName, combinedWinLinearFitStats, fmt = '%.4f',
              header = 'Trial Number\t Direction\t Magnification Gain\t'
              'Perceived Target Force\t Running Sum (noFD)\t Reversals (noFD)',
              delimiter = '\t')

## Finalize and save power spectra statistics
combinedPowerStatistics = numpy.column_stack((trialColumn, combinedPowerStatistics))

powerFileName = subjectDirectory + '\\powerSpectraStatistics.txt'
numpy.savetxt(powerFileName, combinedPowerStatistics, fmt = '%.4f',
                  header = 'Trial Number\t Direction\t Magnification Gain\t'
                  'Perceived Target Force\t Cumulative Power Raw 1-4Hz\t'
                  'Cumulative Power Raw 4-7Hz\t Cumulative Power Raw 7-10Hz\t'
                  'Cumulative Power Contact 1-4Hz\t Cumulative Power Contact 4-7Hz\t'
                  'Cumulative Power Contact 7-10Hz',
                  delimiter = '\t')


print "{} processing completed!".format(subjectNumber)
print "Time elapsed: {0:03f} seconds\n".format(time.time() - start)
print "\n\n*** MinContact_DataStatistics.py - END ***\n\n"
