###############################################################################
##    Randy Lee
##    VIA Lab
##    Department of Bioengineering
##    University of Pittsburgh
##
##    Low Force Experiment Analysis
##
##    Function: Crop individual trial data to take statistics
##      of force traces
##
##  Analysis procedure:
##    (1) LowForceExp_RawDataProcessing to splitt raw data and graph time
##        series and Fourier/Power spectra
##    (2) LowForceExp_DataStatistics to obtain mean force and std for each
##        experimental condition, and graph dependent variable analysis of
##        spectral power
##
##
##
##     Last updated : 20 September 2016
###############################################################################

import os
import numpy
import math
import time
from scipy import stats
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


subjectNumber = 'S6'

rawDataDirectory = 'C:\\Users\\Randy Lee\\Documents\\VIA Lab\\HHFM Data\\'
#rawDataDirectory = 'D:\\Randy Lee\\Documents\\VIA Lab\\HHFM Data'

###############################################################################
##
##    Local script functions
##
###############################################################################

def raw_stats(fullDataArray, trialParameters, trialCounter):
    "Precondition: full data block set; Postcondition: Block statistics file"

    ## store raw data
    # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
    # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
    # [8] bit2; [9] bit3; [10] currentlyInTrial voltage;
    # [11] GS0Force; [12] filteredForce; [13] drift;
    GS0_ForceGrams = fullDataArray[:,11]

    ## Beginning of trial files saved as 20 ticks BEFORE currentlyInTrial
    ## End of files saved as 2500 ticks AFTER currentlyInTrial
    endStamp = len(fullDataArray) - 2500

    ## startstamp for 6s vis feedback, 6s no vis feedback
    startStamp = endStamp - 12000

    print "Trial #{} --> startStamp = {}, endStamp = {}".format(trialCounter,
        startStamp, endStamp)

    targetForce = trialParameters[trialCounter, 2]
    magnification = trialParameters[trialCounter, 1]
    perceivedForce = trialParameters[trialCounter, 3]    
    #if magnification == 1:
    #    perceivedForce = targetForce * 4.0
    #elif magnification == 0:
    #    perceivedForce = targetForce

    # trim raw data for visFD/ noFD sections, skipping 1st second of each section
    # vis feedback, 6s
    forceVisFeedback = GS0_ForceGrams[startStamp + 1000 : endStamp - 6000]

    # no vis feedback, 6s
    forceNoFeedback = GS0_ForceGrams[endStamp - 5000 : endStamp]

    statistics = numpy.zeros(shape = (7))
    statistics[0] = targetForce
    statistics[1] = perceivedForce
    statistics[2] = magnification
    statistics[3] = numpy.mean(forceVisFeedback)
    statistics[4] = numpy.std(forceVisFeedback)
    statistics[5] = numpy.mean(forceNoFeedback)
    statistics[6] = numpy.std(forceNoFeedback)

    return statistics


def filter_stats(fullDataArray, trialParameters, trialCounter):
    "Precondition: filtered force array. Postcondition: stats to add to dataStatistics"

    GS0_ForceGrams = fullDataArray[:,11]
    filteredForceArray = fullDataArray[:,12]
    #driftArray = fullDataArray[:,13]

    # Define start and stop points for trim, extract trial parameters
    endStamp = len(filteredForceArray) - 2500
    startStamp = endStamp - 12000

    targetForce = trialParameters[trialCounter, 2]
    magnification = trialParameters[trialCounter, 1]
    perceivedForce = trialParameters[trialCounter, 3]    
    #if magnification == 1:
    #    perceivedForce = targetForce * 11.0
    #elif magnification == 0:
    #    perceivedForce = targetForce

    # trim raw data based on visFD/ noFD
    filteredForce_visFD = filteredForceArray[startStamp : endStamp - 6000]
    filteredForce_noFD = filteredForceArray[endStamp - 6000 : endStamp]

    #driftForceVisFeedback = driftArray[startStamp : endStamp - 6000]
    #driftForceNoFeedback = driftArray[endStamp - 6000 : endStamp]

    # count time out of contact with sensor
    trialStart_raw_mean = numpy.mean(GS0_ForceGrams[0:250])
    trialStart_raw_std = numpy.std(GS0_ForceGrams[0:250])
    trialStart_filter_mean = numpy.mean(filteredForceArray[0:250])
    trialStart_filter_std= numpy.std(filteredForceArray[0:250])
    
    # in 100ms windows, check to see if filtered force within raw mean & std
    windowSize = 50
    windowStart_visFD = numpy.arange(1000, len(filteredForce_visFD), windowSize)
    sections_noContact_visFD = 0
    for i in range(len(windowStart_visFD)):
        if windowStart_visFD[i] + windowSize > len(filteredForce_visFD):
            window_force = filteredForce_visFD[windowStart_visFD[i] : len(filteredForce_noFD)]
        else:
            window_force = filteredForce_visFD[windowStart_visFD[i] : windowStart_visFD[i] + windowSize]
        
        time_noContact_visFD = 0        
        for j in range(len(window_force)):
            if abs(window_force[j] - trialStart_filter_mean) <= trialStart_filter_std:
                # within mean & std of beginning force 
                time_noContact_visFD += 1
        
        if numpy.sum(time_noContact_visFD) == len(window_force):
            # entire window was out of contact with sensor
            sections_noContact_visFD += 1
    
    windowStart_noFD = numpy.arange(1000, len(filteredForce_noFD), windowSize)
    sections_noContact_noFD = 0
    for i in range(len(windowStart_noFD)):
        if windowStart_noFD[i] + windowSize > len(filteredForce_noFD):
            window_force = filteredForce_noFD[windowStart_noFD[i] : len(filteredForce_noFD)]
        else:
            window_force = filteredForce_noFD[windowStart_noFD[i] : windowStart_noFD[i] + windowSize]
        
        time_noContact_noFD = 0        
        for j in range(len(window_force)):
            if abs(window_force[j] - trialStart_filter_mean) <= trialStart_filter_std:
                # within mean & std of beginning force 
                time_noContact_noFD += 1
                
        if numpy.sum(time_noContact_noFD) == len(window_force):
            # entire window was out of contact with sensor
            sections_noContact_noFD += 1

    print "sections_noContact_visFD = {}\nsections_noContact_noFD = {}".format(sections_noContact_visFD, sections_noContact_noFD)    
    
    statistics = numpy.zeros(shape = (11))
    statistics[0] = targetForce
    statistics[1] = perceivedForce
    statistics[2] = magnification
    statistics[3] = numpy.mean(filteredForce_visFD[1000:])
    statistics[4] = numpy.std(filteredForce_visFD[1000:])
    statistics[5] = numpy.mean(filteredForce_noFD[1000:])
    statistics[6] = numpy.std(filteredForce_noFD[1000:])
    statistics[7] = sections_noContact_visFD
    statistics[8] = sections_noContact_visFD/len(windowStart_visFD)
    statistics[9] = sections_noContact_noFD
    statistics[10] = sections_noContact_noFD/len(windowStart_noFD)
    #statistics[7] = numpy.mean(driftForceVisFeedback[1000:])
    #statistics[8] = numpy.std(driftForceVisFeedback[1000:])
    #statistics[9] = numpy.mean(driftForceNoFeedback[1000:])
    #statistics[10] = numpy.std(driftForceNoFeedback[1000:])

    return statistics

def linear_fit(fullDataArray, trialParameters, trialCounter):
    ''' Calculate linear regression for visFD/ noFD halves of each trial

    '''

    # store raw data
    GS0_ForceGrams = fullDataArray[:,11]

    # extract trial parameters
    endStamp = len(GS0_ForceGrams) - 2500
    startStamp = endStamp - 12000

    targetForce = trialParameters[trialCounter, 2]
    magnification = trialParameters[trialCounter, 1]
    perceivedForce = trialParameters[trialCounter, 3]
    #if magnification == 1:
    #    perceivedForce = targetForce * 4.0
    #elif magnification == 0:
    #    perceivedForce = targetForce

    # skip first second of each half
    rawForce_visFD = GS0_ForceGrams[startStamp + 1000 : endStamp - 6000]
    rawForce_noFD = GS0_ForceGrams[endStamp - 5000 : endStamp]

    x1 = numpy.arange(len(rawForce_visFD))
    x2 = numpy.arange(len(rawForce_noFD))

    slope_visFD, intercept_visFD, r_visFD, p_val, std_err = stats.linregress(x1, rawForce_visFD)
    slope_noFD, intercept_noFD, r_noFD, p_val, std_err = stats.linregress(x2, rawForce_noFD)

    r2_visFD = r_visFD**2
    r2_noFD = r_noFD**2

    statistics = numpy.zeros(shape=(9))
    statistics[0] = targetForce
    statistics[1] = perceivedForce
    statistics[2] = magnification
    statistics[3] = slope_visFD
    statistics[4] = intercept_visFD
    statistics[5] = r2_visFD
    statistics[6] = slope_noFD
    statistics[7] = intercept_noFD
    statistics[8] = r2_noFD

    return statistics

def window_linear_fit(fullDataArray, trialParameters, trialCounter):
    ''' Calculate linear regression for visFD/noFD within small windows
    
    '''
    
    GS0_ForceGrams = fullDataArray[:,11]
    
    # set up window size and shift parameters
    windowSize = 150
    shift = 50
    
    # extract trial parameters
    endStamp = len(GS0_ForceGrams) - 2500
    startStamp = endStamp - 12000
    
    targetForce = trialParameters[trialCounter, 2]
    magnification = trialParameters[trialCounter, 1]
    perceivedForce = trialParameters[trialCounter, 3]
    #if magnification == 1:
    #    perceivedForce = targetForce * 4.0
    #elif magnification == 0:
    #    perceivedForce = targetForce
    
    rawForce_visFD = GS0_ForceGrams[startStamp : endStamp - 6000]
    rawForce_noFD = GS0_ForceGrams[endStamp - 6000 : endStamp]
    
    # analysis for visFD
    windowStart = numpy.arange(0, len(rawForce_visFD), shift)
    
    # [0] window #; [1] start stamp [2] slope [3] r^2  
    runningSum_visFD = 0
    reversal_visFD = 0
    last = 0
    results_visFD = numpy.zeros(shape = (len(windowStart), 5))
    for i in range(len(windowStart)):
        if windowStart[i] + windowSize > len(rawForce_visFD):
            # if current section exceeds length of array, slice only until end
            window_force = rawForce_visFD[windowStart[i] : len(rawForce_visFD)]
        else:
            window_force = rawForce_visFD[windowStart[i] : windowStart[i] + windowSize]
        
        x = numpy.arange(len(window_force))
        slope, intercept, r, p_val, std_err = stats.linregress(x, window_force)
        
        # [0] visFD status, [1] window start stamp [2] lin regression slope
        # [3] r^2
        results_visFD[i, :] = [1, i, windowStart[i], slope, r**2]
        
        if slope >= 0: 
            runningSum_visFD += 1
            # if last window had negative slope, but now positive, count reversal
            if last == -1:
                reversal_visFD += 1
            last = 1
        elif slope < 0:
            runningSum_visFD -= 1
            # if last window had positive slope, but now negative, count reversal
            if last == 1:
                reversal_visFD += 1
            last = -1
            
    # analysis for noFD
    windowStart = numpy.arange(0, len(rawForce_noFD), shift)

    runningSum_noFD = 0
    reversal_noFD= 0
    results_noFD = numpy.zeros(shape = (len(windowStart), 5))
    for i in range(len(windowStart)):
        if windowStart[i] + windowSize > len(rawForce_noFD):
            # if current section exceeds length of array, slice only until end
            window_force = rawForce_noFD[windowStart[i] : len(rawForce_noFD)]
        else:
            window_force = rawForce_noFD[windowStart[i] : windowStart[i] + windowSize]
        
        x = numpy.arange(len(window_force))
        slope, intercept, r, p_val, std_err = stats.linregress(x, window_force)
        
        results_noFD[i, :] = [0, i, windowStart[i], slope, r**2]
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
    results = numpy.vstack((results_visFD, results_noFD))
    savePath = os.getcwd() + '\\win_regression_Trial{0:02d}.txt'.format(trialCounter+1)
    
    numpy.savetxt(savePath, results, fmt = '%.4f', delimiter = '\t',
                  header = 'visFD\t window number\t start stamp\t slope\t R^2')
    
    # return statistics as running sums and reversal totals 
    statistics = numpy.zeros(shape = (7))
    statistics[0] = targetForce
    statistics[1] = perceivedForce
    statistics[2] = magnification
    statistics[3] = runningSum_visFD
    statistics[4] = reversal_visFD
    statistics[5] = runningSum_noFD
    statistics[6] = reversal_noFD
    
    return statistics


def power_stats(fourierCoeff, trialParameters, trialCounter):
    '''Calculate and save power band spectra and cumulative power for each trial

    '''

    # store raw data
    fftRawVisFD = fourierCoeff[:,1]
    fftRawNoFD = fourierCoeff[:,3]
    #freq = fourierCoeff[:,0]

    # extract trial parameters
    targetForce = trialParameters[trialCounter, 2]
    magnification = trialParameters[trialCounter, 1]
    perceivedForce = trialParameters[trialCounter, 3]
    #if magnification == 1:
    #    perceivedForce = targetForce * 4.0
    #elif magnification == 0:
    #    perceivedForce = targetForce

    powerVisFD1_4 = numpy.zeros(shape = (19))
    powerNoFD1_4 = numpy.zeros(shape = (19))
    index = 0
    # calculate power spectra coefficients for freq 1 - 4 Hz as square of
    # fourier coeff
    for i in range(6, 25):
        powerVisFD1_4[index] = fftRawVisFD[i] ** 2
        powerNoFD1_4[index] = fftRawNoFD[i] ** 2
        index += 1

    for index in range(len(powerVisFD1_4)):
        if index == 0:
            cumPowerVisFD1_4 = powerVisFD1_4[index]
            cumPowerNoFD1_4 = powerNoFD1_4[index]
        else:
            cumPowerVisFD1_4 += powerVisFD1_4[index]
            cumPowerNoFD1_4 += powerNoFD1_4[index]

    # calculate power spectra coefficients for freq 4 - 7 Hz as square of
    # fourier coeff
    powerVisFD4_7 = numpy.zeros(shape = (19))
    powerNoFD4_7 = numpy.zeros(shape = (19))
    index = 0
    for i in range(25, 43):
        powerVisFD4_7[index] = fftRawVisFD[i] ** 2
        powerNoFD4_7[index] = fftRawNoFD[i] ** 2
        index += 1

    for index in range (len(powerVisFD1_4)):
        if index == 0:
            cumPowerVisFD4_7 = powerVisFD4_7[index]
            cumPowerNoFD4_7 = powerNoFD4_7[index]
        else:
            cumPowerVisFD4_7 += powerVisFD4_7[index]
            cumPowerNoFD4_7 += powerNoFD4_7[index]

    # calculate power spectra coefficients for freq 7 - 10 Hz as square of
    # fourier coeff
    powerVisFD7_10 = numpy.zeros(shape = (19))
    powerNoFD7_10 = numpy.zeros(shape = (19))
    index = 0
    for i in range(43, 61):
        powerVisFD7_10[index] = fftRawVisFD[i] ** 2
        powerNoFD7_10[index] = fftRawNoFD[i] ** 2
        index += 1

    for index in range (len(powerVisFD1_4)):
        if index == 0:
            cumPowerVisFD7_10 = powerVisFD7_10[index]
            cumPowerNoFD7_10 = powerNoFD7_10[index]
        else:
            cumPowerVisFD7_10 += powerVisFD7_10[index]
            cumPowerNoFD7_10 += powerNoFD7_10[index]

    statistics = numpy.zeros(shape = (9))
    statistics[0] = targetForce
    statistics[1] = perceivedForce
    statistics[2] = magnification
    statistics[3] = cumPowerVisFD1_4
    statistics[4] = cumPowerNoFD1_4
    statistics[5] = cumPowerVisFD4_7
    statistics[6] = cumPowerNoFD4_7
    statistics[7] = cumPowerVisFD7_10
    statistics[8] = cumPowerNoFD7_10

    return statistics


def window_time_stats(fullDataArray, trialParameters, trialCounter):
    ''' Calculate the amount of time spent within the target zone for visFD/noFD

    '''
    GS0_ForceGrams = fullDataArray[:,11]

    # extract trial parameters and set up upper/lower bounds
    targetForce = trialParameters[trialCounter, 2]
    magnification = trialParameters[trialCounter, 1]
    perceivedForce = trialParameters[trialCounter, 3]
    #if magnification == 1:
    #    perceivedForce = targetForce * 4.0
    #elif magnification == 0:
    #    perceivedForce = targetForce

    upperBound = targetForce + 0.5
    lowerBound = targetForce - 0.5

    # set up slicing of raw data into visFD/noFD
    endStamp = len(GS0_ForceGrams) - 2500
    startStamp = endStamp - 12000

    GS0_visFD = GS0_ForceGrams[startStamp : endStamp - 6000]
    GS0_noFD = GS0_ForceGrams[endStamp - 6000 : endStamp]

    # count time in window, visFD
    timeInWindow = 0.0
    for timestamp in range(len(GS0_visFD)):
        if (lowerBound <= GS0_visFD[timestamp] <= upperBound) == True:
            timeInWindow += 1.0
    
    print 'time_inWindow_visFD: {}'.format(timeInWindow)
    inWinProportion_visFD = float(timeInWindow/len(GS0_visFD))

    # count time in window, noFD
    timeInWindow = 0.0
    for timestamp in range(len(GS0_noFD)):
        if (lowerBound <= GS0_noFD[timestamp] <= upperBound) == True:
            timeInWindow += 1.0
    
    print 'time_inWindow_noFD: {}\n'.format(timeInWindow)
    inWinProportion_noFD = float(timeInWindow/len(GS0_noFD))

    statistics = numpy.zeros(shape = (5))

    statistics[0] = targetForce
    statistics[1] = perceivedForce
    statistics[2] = magnification
    statistics[3] = inWinProportion_visFD
    statistics[4] = inWinProportion_noFD

    return statistics
    
def RMS_error(trialData, trialParameters, trialCounter):
    
    GS0_ForceGrams = trialData[:,11]

    targetForce = trialParameters[trialCounter, 2]
    magnification = trialParameters[trialCounter, 1]
    perceivedForce = trialParameters[trialCounter, 3]
    #if magnification == 1:
    #    perceivedForce = targetForce * 4.0
    #elif magnification == 0:
    #    perceivedForce = targetForce
    
    endStamp = len(GS0_ForceGrams) - 2500
    startStamp = endStamp - 12000
    
    # slice raw force data with first second of each half skipped
    GS0Force_visFD = GS0_ForceGrams[startStamp + 1000 : endStamp - 6000]    
    GS0Force_noFD = GS0_ForceGrams[endStamp - 5000 : endStamp]
    visFD_mean = numpy.mean(GS0Force_visFD)
    noFD_mean = numpy.mean(GS0Force_noFD)
    
    error_visFD = numpy.zeros(shape = len(GS0Force_visFD))
    errorSq_visFD = numpy.zeros(shape = len(GS0Force_visFD))
    error_noFD = numpy.zeros(shape = len(GS0Force_noFD))
    errorSq_noFD = numpy.zeros(shape = len(GS0Force_noFD))
    
    errorSqSum_visFD = 0.0
    for i in range(len(errorSq_visFD)):
        error_visFD[i] = GS0Force_visFD[i] - visFD_mean
        errorSq_visFD[i] = error_visFD[i]**2.0
        errorSqSum_visFD += errorSq_visFD[i]
    
    errorSqSum_noFD = 0.0
    for i in range(len(errorSq_noFD)):
        error_noFD[i] = GS0Force_noFD[i] - noFD_mean
        errorSq_noFD[i] = error_noFD[i]**2.0
        errorSqSum_noFD += errorSq_noFD[i]
    
    errorRMS_visFD = math.sqrt(errorSqSum_visFD/5000.0)
    errorRMS_noFD = math.sqrt(errorSqSum_noFD/5000.0)
    
    results = numpy.zeros(shape = (9))
    results[0] = targetForce
    results[1] = perceivedForce
    results[2] = magnification
    results[3] = numpy.mean(error_visFD)
    results[4] = numpy.std(error_visFD)
    results[5] = errorRMS_visFD
    results[6] = numpy.mean(error_noFD)
    results[7] = numpy.std(error_noFD)
    results[8] = errorRMS_noFD
    
    return results


def graph_window_fits(trialParameters, blockNumberString):
    ''' Graph windowed linear regression fits with visFD/noFD and mag on/off together
    
    '''
    
    # open PDF for plots
    pdfFileName = os.getcwd() + '\\' + blockNumberString + '_windowedFitGraphs.pdf'
    pp = PdfPages(pdfFileName)
    
    # array with ordered list of file names
    blockDataList = numpy.empty(16, dtype = 'str')
    blockDataList = reorder_trials(blockNumberString)
    
    trialCounter = 0
    figCounter = 1
    doneFlag = 0
    print "Graphing windowed linear fits..."
    while(doneFlag == 0):
        
        # extract magnified data
        file_name = os.getcwd() + '\\win_regression_{}'.format(blockDataList[trialCounter])
        fullData_mag = numpy.genfromtxt(file_name, skip_header = 1)
        
        # extract unmagnified data
        file_name = os.getcwd() + '\\win_regression_{}'.format(blockDataList[trialCounter + 1])
        fullData_nomag = numpy.genfromtxt(file_name, skip_header = 1)
        
        # seperate visFD from noFD for mag and nomag
        conditions = [fullData_mag[:,0] == 1, fullData_mag[:,0] == 1]
        choices = [fullData_mag[:,2], fullData_mag[:,3]]
        indx_mag_visFD = numpy.trim_zeros(numpy.select(conditions, choices), 'b')
        #print indx_mag_visFD
        
        conditions = [fullData_mag[:,0] == 0, fullData_mag[:,0] == 0]
        indx_mag_noFD = numpy.trim_zeros(numpy.select(conditions, choices), 'b')
        #print indx_mag_noFD
        
        conditions = [fullData_nomag[:,0] == 1]
        choices = [fullData_nomag[:,2]] 
        indx_nomag_visFD = numpy.trim_zeros(numpy.select(conditions, choices), 'b')
        
        conditions = [fullData_nomag[:,0] == 0]
        indx_nomag_noFD = numpy.trim_zeros(numpy.select(conditions, choices), 'b')
        
        # create figure and begin plotting        
        fig = plt.figure(figCounter)
        print "Figure #{}...\n".format(figCounter)
        figureTitle = 'Trial#{}; (Prototype #{:d}&{:d}: Target Force = {})'.format(
                        figCounter, int(trialParameters[trialCounter][0]),
                        int(trialParameters[trialCounter+1][0]),
                        trialParameters[trialCounter][2])
        fig.suptitle(figureTitle, fontsize = 12, fontweight = 'bold')
                                
        ax1 = plt.subplot(211)
        
        subplotTitle = 'Magnification On'
        ax1.set_title(subplotTitle, fontsize = 8)
        
        x1 = fullData_mag[0:len(indx_mag_visFD),2]
        y1 = fullData_mag[0:len(indx_mag_visFD),3]
        x2 = fullData_mag[len(indx_mag_visFD) : len(indx_mag_visFD) + len(indx_mag_noFD),2]
        y2 = fullData_mag[len(indx_mag_visFD) : len(indx_mag_visFD) + len(indx_mag_noFD),3]
        
        mag_visFD, = plt.plot(x1, y1, 'b-', lw = 1.00) 
        mag_noFD, = plt.plot(x2, y2, 'r-', lw = 1.00)
        
        plt.setp(ax1.get_xticklabels(), visible = False)
        ax1.set_ylabel('Slope')
        ax1.set_xlim((0,6500))
        ax1.set_ylim((-0.15, 0.15))
        ax1.axhline(y=0, color = 'k', lw = 1.00)
        ax1.grid(True, which = 'both')
        plt.legend((mag_visFD, mag_noFD), ('visFD', 'noFD'), 
                   loc = 'lower left', fontsize = 8)

        ax2 = plt.subplot(212)
        subplotTitle = 'Magnification Off'
        ax2.set_title(subplotTitle, fontsize = 8)
        
        x1 = fullData_nomag[0:len(indx_nomag_visFD),2]
        y1 = fullData_nomag[0:len(indx_nomag_visFD),3]
        x2 = fullData_nomag[len(indx_nomag_visFD) : len(indx_nomag_visFD) + len(indx_nomag_noFD),2]
        y2 = fullData_nomag[len(indx_nomag_visFD) : len(indx_nomag_visFD) + len(indx_nomag_noFD),3]
        
        nomag_visFD, = plt.plot(x1, y1, 'b-', lw = 1.00)
        nomag_noFD, = plt.plot(x2, y2, 'r-', lw = 1.00)
        
        plt.legend((nomag_visFD, nomag_noFD), ('visFD', 'noFD'),
                   loc = 'lower left', fontsize = 8)
        ax2.set_ylabel('Slope')
        ax2.set_xlabel('Window Start')
        ax2.set_xlim((0, 6500))
        ax2.set_ylim((-0.15, 0.15))
        ax2.axhline(y = 0, color = 'k', lw = 1.00)
        ax2.grid(True, which = 'both')
        
        pp.savefig(figCounter)
        plt.close('all')
        
        # advance trial/fig counters
        trialCounter += 2
        figCounter += 1
        
        if trialCounter >= 16:
            plt.close('all')
            pp.close()
            print "\nPDF closed.\n"
            print "{} slope data graphing completed!\n".format(blockNumberString)
            doneFlag = 1
    
    
    return


def reorder_trials(blockNumberString):

    ## List of prototype numbers, ordered by force, then unmag/mag
    # [1, 9 , 2, 10, 3, 11, 4, 12, 5, 13, 6, 14, 7, 15, 8, 16]

    ## blockNumberString passed as S#_BlockX

    if blockNumberString[-6:] == 'Block1':
        orderedTrialList = numpy.array(['Trial01.txt', 'Trial06.txt', 'Trial16.txt',
                            'Trial03.txt', 'Trial15.txt', 'Trial07.txt',
                            'Trial05.txt', 'Trial04.txt', 'Trial11.txt',
                            'Trial08.txt', 'Trial10.txt', 'Trial14.txt',
                            'Trial02.txt', 'Trial13.txt', 'Trial12.txt',
                            'Trial09.txt'])


    elif blockNumberString[-6:] == 'Block2':
        orderedTrialList = numpy.array(['Trial12.txt', 'Trial06.txt', 'Trial03.txt',
                            'Trial14.txt', 'Trial09.txt', 'Trial13.txt',
                            'Trial04.txt', 'Trial02.txt', 'Trial08.txt',
                            'Trial01.txt', 'Trial11.txt', 'Trial16.txt',
                            'Trial05.txt', 'Trial07.txt', 'Trial10.txt',
                            'Trial15.txt'])

    elif blockNumberString[-6:] == 'Block3':
        orderedTrialList = numpy.array(['Trial14.txt', 'Trial12.txt', 'Trial10.txt',
                            'Trial13.txt', 'Trial02.txt', 'Trial07.txt',
                            'Trial01.txt', 'Trial15.txt', 'Trial09.txt',
                            'Trial04.txt', 'Trial05.txt', 'Trial06.txt',
                            'Trial16.txt', 'Trial08.txt', 'Trial11.txt',
                            'Trial03.txt'])

    elif blockNumberString[-6:] == 'Block4':
        orderedTrialList = numpy.array(['Trial08.txt', 'Trial13.txt', 'Trial02.txt',
                            'Trial01.txt', 'Trial15.txt', 'Trial11.txt',
                            'Trial16.txt', 'Trial07.txt', 'Trial09.txt',
                            'Trial14.txt', 'Trial05.txt', 'Trial03.txt',
                            'Trial06.txt', 'Trial12.txt', 'Trial04.txt',
                            'Trial10.txt'])

    elif blockNumberString[-6:] == 'Block5':
        orderedTrialList = numpy.array(['Trial01.txt', 'Trial03.txt', 'Trial16.txt',
                            'Trial10.txt', 'Trial08.txt', 'Trial12.txt',
                            'Trial07.txt', 'Trial14.txt', 'Trial02.txt',
                            'Trial05.txt', 'Trial11.txt', 'Trial06.txt',
                            'Trial13.txt', 'Trial04.txt', 'Trial15.txt',
                            'Trial09.txt'])

    return orderedTrialList


###############################################################################
##
##    MAIN: Process Data
##
###############################################################################

start = time.time()

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
    #blockDataFiles = os.listdir( os.getcwd() )

    trialCounter = 0
    dataStatistics = numpy.zeros(shape = (16,7))
    filterStats = numpy.zeros(shape = (16,11))
    linearFitStats = numpy.zeros(shape = (16,9))
    winLinearFitStats = numpy.zeros(shape = (16,7))
    powerStats = numpy.zeros(shape = (16,9))
    windowTime = numpy.zeros(shape = (16,5))
    RMSerror = numpy.zeros(shape = (16,9))

    blockDataList = numpy.empty(16, dtype = 'str')
    blockDataList = reorder_trials(blockFolderName)

    for trialCounter in range(0,16):
        # Choose trial data from ORDERED list of trial parameters
        fileName = blockDataList[trialCounter]
        print "Accessing {} of {}...".format(fileName, blockFolderName[-6:])
        trialPath = blockFolder + '\\' + fileName
        trialData = numpy.genfromtxt(trialPath, skip_header = 1)

        # Obtain fourier coefficients, already ordered by distal target
        fileName = blockFolderName + '_{}_FourierCoeff.txt'.format(trialCounter + 1)
        trialPath = blockFolder + '\\' + fileName
        fourierCoeff = numpy.genfromtxt(trialPath, skip_header = 1)

        # Return force statistics from cropped trials
        dataStatistics[trialCounter] = raw_stats(trialData,
                                                 trialParameters, trialCounter)

        # Filtered force statistics
        filterStats[trialCounter] = filter_stats(trialData, 
                                                 trialParameters, trialCounter)

        # Linear regression analysis
        linearFitStats[trialCounter] = linear_fit(trialData, trialParameters,
                                                  trialCounter)
        
        # Windowed linear regression analysis
        winLinearFitStats[trialCounter] = window_linear_fit(trialData,
                                                trialParameters, trialCounter)

        # Power spectra statistics
        powerStats[trialCounter] = power_stats(fourierCoeff, trialParameters,
                                               trialCounter)

        # Time spent in target window
        windowTime[trialCounter] = window_time_stats(trialData, trialParameters,
                                                trialCounter)
        
        # Calculate RMS error
        RMSerror[trialCounter] = RMS_error(trialData, trialParameters, trialCounter)

    # graph windowed linear fit data
    graph_window_fits(trialParameters, blockFolderName)

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

        combinedWindowTimeStatistics = numpy.zeros(shape = (windowTime.shape))
        combinedWindowTimeStatistics = windowTime
        
        combinedRMSError = numpy.zeros(shape = (RMSerror.shape))
        combinedRMSError = RMSerror

    else:
        combinedStatistics = numpy.vstack((combinedStatistics, dataStatistics))

        combinedFilterStatistics = numpy.vstack((combinedFilterStatistics,
            filterStats))

        combinedLinearFits = numpy.vstack((combinedLinearFits, linearFitStats))
        
        combinedWinLinearFitStats = numpy.vstack((combinedWinLinearFitStats,
                                                  winLinearFitStats))

        combinedPowerStatistics = numpy.vstack((combinedPowerStatistics,
                                                powerStats))

        combinedWindowTimeStatistics = numpy.vstack((combinedWindowTimeStatistics,
                                                     windowTime))
        
        combinedRMSError = numpy.vstack((combinedRMSError, RMSerror))

    # Graph stats data for entire block
    #graphBlockData(dataStatistics)
    #graphFilteredData()
    os.chdir(subjectDirectory)

## Finalize and save window time statistics
trialColumn = numpy.array(numpy.arange(1,81))
combinedWindowTimeStatistics = numpy.column_stack((trialColumn,
                                                   combinedWindowTimeStatistics))

windowTimeFileName = subjectDirectory + '\\windowTimeStatistics.txt'
with open(windowTimeFileName, 'w+') as windowTimeFile:
    numpy.savetxt(windowTimeFileName, combinedWindowTimeStatistics, fmt = '%.4f',
                  header = 'Trial Number\t Distal Target Force\t'
                  'Magnified Target Force\t Magnification On(1)/Off(0)\t'
                  '% Time in Window (visFD)\t % Time in Window (noFD)',
                  delimiter = '\t')

## Finalize and save experiment statistics
combinedStatistics = numpy.column_stack((trialColumn, combinedWindowTimeStatistics[:,4:],
                                         combinedStatistics))

statisticsFileName = subjectDirectory + '\\experimentStatistics.txt'
#statisticsFile = open(statisticsFileName , 'w+' )
with open(statisticsFileName , 'w+' ) as statisticsFile:
    numpy.savetxt(statisticsFileName, combinedStatistics , fmt= '%.4f',
                  header = 'Trial Number\t % Time in Window (visFD)\t'
                  '% Time in Window (noFD)\t Distal Target Force\t Magnified Target Force\t'
                  'Magnification On(1)/Off(0)\t Mean Force (vis feedback)\t'
                  'SD Force (vis feedback)\t Mean Force (no feedback)\t'
                  'SD Force (no feedback)',
                  delimiter = '\t')
#statisticsFile.close()


## Finalize and save filtered force statistics
combinedFilterStatistics = numpy.column_stack((trialColumn, combinedWindowTimeStatistics[:,4:],
                                               combinedFilterStatistics))

filteredFileName = subjectDirectory + '\\filteredForceStatistics.txt'
#filteredFile = open(filteredFileName, 'w+')
with open(filteredFileName, 'w+') as filteredFile:
    numpy.savetxt(filteredFileName, combinedFilterStatistics , fmt= '%.4f',
                  header = 'Trial Number\t % Time in Window (visFD)\t'
                  '% Time in Window (noFD)\t Distal Target Force\t Magnified Target Force\t'
                  'Magnification On(1)/Off(0)\t Mean Force (vis feedback)\t'
                  'SD Force (vis feedback)\t Mean Force (no feedback)\t'
                  ' SD Force (no feedback)\t sections_noContact visFD\t'
                  '% sections_noContact visFD\t sections_noContact noFD\t'
                  '% sections_noContact noFD',
                  delimiter = '\t')
#filteredFile.close()

## Finalize and save linear fit statistics
combinedLinearFits = numpy.column_stack((trialColumn, combinedWindowTimeStatistics[:,4:],
                                         combinedLinearFits))

linearFitsFileName = subjectDirectory + '\\linearFitStatistics.txt'
#linearFitsFile = open(linearFitsFileName, 'w+')
with open(linearFitsFileName, 'w+') as linearFitsFile:
    numpy.savetxt(linearFitsFileName, combinedLinearFits, fmt = '%.6f',
                  header = 'Trial Number\t %Time in Window (visFD)\t'
                  '% Time in Window (noFD)\t Distal Target Force\t'
                  'Magnified Target Force\t Magnification On(1)/Off(0)\t'
                  'visFD slope\t visFD intercept\t visFD R^2\t'
                  'noFD slope\t noFD intercept\t noFD R^2',
                  delimiter = '\t')
#linearFitsFile.close()

## Finalize and save windowed linear fit summary stats
combinedWinLinearFitStats = numpy.column_stack((trialColumn, combinedWindowTimeStatistics[:,4:],
                                                combinedWinLinearFitStats))

winLinearFitsFileName = subjectDirectory + '\\windowedLinearFitStatistics.txt'
with open(winLinearFitsFileName, 'w+') as winLinearFitsFile:
    numpy.savetxt(winLinearFitsFileName, combinedWinLinearFitStats, fmt = '%.4f',
                  header = 'Trial Number\t % Time in Window (visFD)\t'
                  '% Time in Window (noFD)\t Target Force\t Perceived Force\t'
                  'Magnification On(1)/ Off(0)\t Running Sum (visFD)\t'
                  'Reversals (visFD)\t Running Sum (noFD)\t Reversals (noFD)',
                  delimiter = '\t')

## Finalize and save power spectra statistics
combinedPowerStatistics = numpy.column_stack((trialColumn, combinedWindowTimeStatistics[:,4:],
                                              combinedPowerStatistics))
powerFileName = subjectDirectory + '\\powerSpectraStatistics.txt'
#powerFile = open(powerFileName, 'w+')
with open(powerFileName, 'w+') as powerFile:
    numpy.savetxt(powerFileName, combinedPowerStatistics, fmt = '%.4f',
                  header = 'Trial Number\t % Time in Window (visFD)\t '
                  '% Time in Window (noFD)\t Distal Target Force\t '
                  'Magnified Target Force\t Magnification On(1)/Off(0)\t'
                  'Cumulative Power- visFD 1-4Hz\t Cumulative Power- noFD 1-4Hz\t'
                  'Cumulative Power- visFD 4-7Hz\t Cumulative Power- noFD 4-7Hz\t'
                  'Cumulative Power- visFD 7-10Hz\t Cumulative Power- noFD 7-10Hz',
                  delimiter = '\t')
#powerFile.close()

## Finalize and save RMS error statistics
combinedRMSError = numpy.column_stack((trialColumn, combinedWindowTimeStatistics[:,4:],
                                       combinedRMSError))
RMSErrorFileName = subjectDirectory + '\\RMS_errorStats.txt'
with open(RMSErrorFileName, 'w+') as windowTimeFile:
    numpy.savetxt(RMSErrorFileName, combinedRMSError, fmt = '%.4f',
                  header = 'Trial Number\t % Time in Window (visFD)\t '
                  '% Time in Window (noFD)\t Distal Target Force\t '
                  'Magnified Target Force\t Magnification On(1)/Off(0)\t'
                  'RMS error visFD_mean\t RMS error visFD_std\t RMS error visFD sum\t'
                  'RMS error noFD_mean\t RMS error noFD_std\t RMS error noFD sum',
                  delimiter = '\t')

print "{} processing completed!\n".format(subjectNumber)
print "Time elapsed: {0:03f} seconds\n".format(time.time() - start)
