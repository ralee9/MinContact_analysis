##############################################################################
##  Randy Lee
##  VIA Lab
##  Department of Bioengineering
##  University of Pittsburgh
##
##  GS0 Low Force Experiment Raw Data Processing
##
##  Function: Separate and plot trial data for blocks of low force exp trials
##
##  Analysis procedure:
##    (1) LowForceExp_RawDataProcessing to splitt raw data and graph time
##        series and Fourier/Power spectra
##    (2) LowForceExp_DataStatistics to obtain mean force and std for each
##        experimental condition, and graph dependent variable analysis of
##        spectral power
##
##   NI-6009 Analog input:
##  - AI0.0 : FS0-100 force sensor; from ext. hardware (0 - 8V)
##  - AI0.1 : HHFM solenoid voltage; from HHFM hardware (-10V - 10V)
##  - AI0.2 : HHFM force sensor voltage; from HHFM hardware (0 - 10V)
##  - AI0.3 : Trial number bit 0; from exp. ADUC7026
##  - AI0.4 : Trial number bit 1; from exp. ADUC7026
##  - AI0.5 : Trial number bit 2; from exp. ADUC7026
##  - AI0.6 : Trial number bit 3; from exp. ADUC7026
##  - AI0.7 : Currently in trial; from exp. ADUC7026
##            (i.e. 0 = before/after trial; 1 = currently in trial)
##                  or
##            Threshold felt; from exp. ADUC7026
##
##  Last updated: 16 September 2016
##############################################################################

import os
import math
import time
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

subjectNumber = 'S5' # Only need to change subjectNumber to process all data

rawDataDirectory = 'C:\\Users\\Randy Lee\\Documents\\VIA Lab\\HHFM Data\\LFE_testing\\LFE_max_4g\\'
#rawDataDirectory = 'D:\\Randy Lee\\Documents\\VIA Lab\\HHFM Data\\'

##############################################################################
##
## Local script functions
##
##############################################################################


def split_trials(fullDataArray, blockNumber):
    '''Create text files for individual trials using raw block data.

    Uses 2D array obtained from raw data files and splits using voltage on
    currentlyInTrial pin (col 9). Copies section from raw data and includes
    conversion to force from GS0 sensor, filtered force, and drift data.
    Trials are unordered -- #1-16 are in the order as they were presented in
    the block.
    '''

    ## Generate empty data arrays and flags
    fullDataSize = fullDataArray.shape
    timestampArray = numpy.zeros(shape = (33, 1))  # at most 32 (16 trials)
    GS0Force = numpy.zeros(shape = (fullDataSize[0], 1))
    timestampCounter = 0
    #lastInTrialBinary = 0

    ## Scan through full raw data set to look for start/end of trials
    for i in range(0, fullDataSize[0]):
        # Convert GS0 voltage to grams force
        # Abs value because voltage can dip to negative small numbers
        # Calibration = (20g/200mV) to quantify force (cal factor for push)
        # Negative to match direction of target force, i.e. neg = push, pos = pull
        currentInTrialVoltage = abs(numpy.float(fullDataArray[i, 9]))
        GS0Force[i] = (numpy.float(fullDataArray[i][2])) * -(20/0.200)

        ## Identify timestamps where trials start/stop
        if(i == 0):
            # To prevent false alarm on first data point
            lastInTrialVoltage = currentInTrialVoltage

        inTrialDiff = lastInTrialVoltage - currentInTrialVoltage

        if(inTrialDiff >= 2):
            # ex. 3.3V --> 0; trial just started
            # Save as 20 timestamps before actual change
            timestampArray[timestampCounter] = (i - 20)
            print "startTime = {}\t".format(timestampArray[timestampCounter])
            timestampCounter = timestampCounter + 1

        if(inTrialDiff <= -2):
            # ex. 0V --> 3.3; trial just ended
            # Save as 2500 timestamps after actual end
            timestampArray[timestampCounter] = (i + 2500)
            print "endTime = {}\n".format(timestampArray[timestampCounter])
            timestampCounter = timestampCounter + 1

        lastInTrialVoltage = currentInTrialVoltage
        
        if timestampCounter > 32:
            # stop splitting once we read all 16 trials
            break

    # Add GS0Force (units = grams) as column [10]
    fullDataArray = numpy.column_stack((fullDataArray, GS0Force))

    # Filter forces with 11 element gaussian (Takes full raw data array)
    # Drift is raw - filtered
    print "Performing gaussian filtering..."
    filteredForces, drift = filter_forces(fullDataArray)
    print "Gaussian filtering done.\n"

    # Filtered force (grams) as column 11, and drift (grams) as column 12
    # Get new data size, with these new
    tempStack = numpy.column_stack((filteredForces, drift))
    fullDataArray = numpy.column_stack((fullDataArray, tempStack))
    #fullDataArray = numpy.column_stack((fullDataArray, filteredForces))
    fullDataSize = fullDataArray.shape

    print "Writing individual trial files...\n"

    # Reset timestamp location counter for each new block
    timestampCounter = 0
    ## Create individual trial files for this block
    for fileCounter in range(1, 17):
        ## 16 trials
        ## start going through timestamp holder
        startTime = numpy.int(timestampArray[timestampCounter])
        endTime = numpy.int(timestampArray[timestampCounter + 1])

        # rows = newDataSize
        # columns = fullDataSize[1] + 1 because replace 2 NaN with block# and
        # trial#, timestamp, then rest of voltage/force data
        newDataSize = endTime - startTime
        newData = numpy.zeros(shape = (newDataSize, fullDataSize[1] + 1))

        ## Begin writing the data files
        trialDataFileName = newBlockDir + '\\Trial{0:02d}.txt'.format(fileCounter)
        trialDataFile = open(trialDataFileName, 'w+')
        for newDataLocation in range(0, newDataSize):
            # Raw data location
            rawDataLocation = startTime + newDataLocation

            ## Copy data into new trial file
            # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
            # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
            # [8] bit2; [9] bit3; [10] currentlyInTrial voltage;
            # [11] GS0Force; [12] filteredForce; [13] drift;
            newData[newDataLocation][0] = blockNumber
            newData[newDataLocation][1] = fileCounter
            #newData[ newDataLocation ][ 0 ] = fullDataArray[ startTime+newDataLocation ][ 0 ]    #date-->Nan
            #newData[newDataLocation][1] = fullDataArray[startTime + newDataLocation][1]    #time-->Nan
            newData[newDataLocation][2] = numpy.int(newDataLocation)        # timestamp
            newData[newDataLocation][3] = numpy.float(fullDataArray[rawDataLocation][2])   # GS0-100 Force Voltage
            newData[newDataLocation][4] = numpy.float(fullDataArray[rawDataLocation][3])   # HHFM Solenoid Voltage
            newData[newDataLocation][5] = numpy.float(fullDataArray[rawDataLocation][4])   # HHFM Force Sensor
            newData[newDataLocation][6] = numpy.float(fullDataArray[rawDataLocation][5])   # trialNumber bit0
            newData[newDataLocation][7] = numpy.float(fullDataArray[rawDataLocation][6])   # trialNumber bit1
            newData[newDataLocation][8] = numpy.float(fullDataArray[rawDataLocation][7])   # trialNumber bit2
            newData[newDataLocation][9] = numpy.float(fullDataArray[rawDataLocation][8])   # trialNumber bit3
            newData[newDataLocation][10] = numpy.float(fullDataArray[rawDataLocation][9])   # currentlyInTrial
            newData[newDataLocation][11] = numpy.float(fullDataArray[rawDataLocation][10]) # calibration = (100g/8V) for GS0-100
            newData[newDataLocation][12] = numpy.float(fullDataArray[rawDataLocation][11]) # Force (grams) filtered with gaussian
            newData[newDataLocation][13] = numpy.float(fullDataArray[rawDataLocation][12]) # Drift from local average (grams)

        numpy.savetxt(trialDataFile, newData, fmt='%3.6f',
            header = 'Block Number\t Trial Number\t Timestamp\t'
                'GS0-100 Voltage\t HHFM Actuator Voltage\t'
                'HHFM Sensor Voltage\t trialNum bit0\t trialNum bit1\t'
                'trialNum bit2\t trialNum bit3\t currentlyInTrial\t'
                'GS0-100 Force (grams)\t Filtered force (grams)\t'
                'High Pass Noise (grams)\t',
            delimiter = '\t')
        trialDataFile.write('\n')
        trialDataFile.close()

        # Advance timestampCounter to next start/end pair
        timestampCounter += 2

    print "Block{} individual trials complete.\n".format(blockNumber)
    return fullDataArray, timestampArray

def get_trial_parameters(trialData):
    '''Return parameter array from binary representation of trial number

    Takes trial-segmented data to extract the prototype ID and return an array
    describing the trial parameters: prototype number, magnification state,
    distal force, and perceived force.
    '''

    # [0] trialPrototypeNum, [1] magnification, [2] targetForce, [3] perceivedForce
    parameters = numpy.zeros(shape = 4)
    magnification = 2.0

    ## Extract binary encoding at 25 timestamps after inTrial pulses (i.e. @ +25)
    ## 3.3V = 1; 0V = 0
    bit0 = getBinaryFromVoltage(trialData[75][6])    # AI0.3
    bit1 = getBinaryFromVoltage(trialData[75][7])    # AI0.4
    bit2 = getBinaryFromVoltage(trialData[75][8])    # AI0.5
    bit3 = getBinaryFromVoltage(trialData[75][9])    # AI0.6

    ## Concatenate bits to obtain trial number
    trialNumber = (bit0 * (1.)) + (bit1 * (2.)) + (bit2 * (4.)) + (bit3 * (8.))

    # Add 1 to trialNumber since binary rep is 1 less
    parameters[0] = trialNumber + 1
    print "bit0 = {}, bit1 = {}, bit2 = {}, bit3 = {}".format(bit0, bit1, bit2, bit3)
    print "--> trialNumber = {}\n".format(parameters[0])

    ## Trials grouped by distal force, neg = push, pos = pull
    ## 1g @ sensor, push, mag ON --> perceived = 11g
    if(parameters[0] == 1):
        parameters[2] = (-1.0)
        parameters[1] = (1)
        parameters[3] = (magnification + 1) * parameters[2]
    ## 1g @ sensor, pull, mag ON --> perceived = 11g
    if(parameters[0] == 2):
        parameters[2] = (1.0)
        parameters[1] = (1)
        parameters[3] = (magnification + 1) * parameters[2]
    ## 2g @ sensor, push, mag ON --> perceived = 22g
    if(parameters[0] == 3):
        parameters[2] = (-2.0)
        parameters[1] = (1)
        parameters[3] = (magnification + 1) * parameters[2]
    ## 2g @ sensor, pull, mag ON --> perceived = 22g
    if(parameters[0] == 4):
        parameters[2] = (2.0)
        parameters[1] = (1)
        parameters[3] = (magnification + 1) * parameters[2]
    ## 3g @ sensor, push, mag ON --> perceived = 33g
    if(parameters[0] == 5):
        parameters[2] = (-3.0)
        parameters[1] = (1)
        parameters[3] = (magnification + 1) * parameters[2]
    ## 3g @ sensor, pull, mag ON --> perceived = 33g
    if(parameters[0] == 6):
        parameters[2] = (3.0)
        parameters[1] = (1)
        parameters[3] = (magnification + 1) * parameters[2]
    ## 4g @ sensor, push, mag ON --> perceived = 44g
    if(parameters[0] == 7):
        parameters[2] = (-4.0)
        parameters[1] = (1)
        parameters[3] = (magnification + 1) * parameters[2]
    ## 4g @ sensor, pull, mag ON --> perceived = 44g
    if(parameters[0] == 8):
        parameters[2] = (4.0)
        parameters[1] = (1)
        parameters[3] = (magnification + 1) * parameters[2]
    ## 1g @ sensor, push, mag OFF --> perceived = 1g
    if(parameters[0] == 9):
        parameters[2] = (-1.0)
        parameters[1] = (0)
        parameters[3] = parameters[2]
    ## 1g @ sensor, pull, mag OFF --> perceived = 1g
    if(parameters[0] == 10):
        parameters[2] = (1.0)
        parameters[1] = (0)
        parameters[3] = parameters[2]
    ## 2g @ sensor, push, mag OFF --> perceived = 2g
    if(parameters[0] == 11):
        parameters[2] = (-2.0)
        parameters[1] = (0)
        parameters[3] = parameters[2]
    ## 2g @ sensor, pull, mag OFF --> perceived = 2g
    if(parameters[0] == 12):
        parameters[2] = (2.0)
        parameters[1] = (0)
        parameters[3] = parameters[2]
    ## 3g @ sensor, push, mag OFF --> perceived = 3g
    if(parameters[0] == 13):
        parameters[2] = (-3.0)
        parameters[1] = (0)
        parameters[3] = parameters[2]
    ## 3g @ sensor, pull, mag OFF --> perceived = 3g
    if(parameters[0] == 14):
        parameters[2] = (3.0)
        parameters[1] = (0)
        parameters[3] = parameters[2]
    ## 4g @ sensor, push, mag OFF --> perceived = 4g
    if(parameters[0] == 15):
        parameters[2] = (-4.0)
        parameters[1] = (0)
        parameters[3] = parameters[2]
    ## 4g @ sensor, pull, mag OFF --> perceived = 4g
    if(parameters[0] == 16):
        parameters[2] = (4.0)
        parameters[1] = (0)
        parameters[3] = parameters[2]

    return parameters


def process_split_trials(pdfSaveFolderPath, timestampArray, blockNumberString):
    '''Graph previously split trials using array of start/stop timestamps'''

    ## currently in HHFM_Data/subjectNumber/Blockx/
    ## make list of block directory
    # 1-D array with list of all file names (split trials)
    #blockDataList = os.listdir(os.getcwd())

    # 1-D array with ORDERED list of split file names
    print "blockNumberString = {}\n".format(blockNumberString[2:])
    blockDataList = numpy.empty(16, dtype = 'str')
    blockDataList = reorder_trials(blockNumberString)
    print "blockDataList = {}\n".format(blockDataList)

    # Obtain number of files in subjectData directory
    blockDataListSize = len(blockDataList)

    ## open PDF for plots and figures
    pdfFileName = os.getcwd() +'\\'+ blockNumberString + '_DataPlots.pdf'
    pp = PdfPages(pdfFileName)            #Pointer to pdf file

    ## create empty arrays and flags
    print "Reading trial data...\n"

    # 16 trials, 3 columns each
    # [0] trialNumber [1] mag [2] targetForce [3] perceivedForce
    trialParameters = numpy.zeros(shape = (17,4))

    # Flag to iterate through all files in block
    doneFlag = 0
    trialCounter = 0

    ## Graph and obtain trial parameters for all trials in block
    # fileCounter = 1 because blockX_timestamps comes first
    # fileCounter = 0 with ordered list because only trials in list
    fileCounter = 0
    figCounter = 1
    while(doneFlag == 0):
        # Print trial number to screen
        print "Trial #{} graphing...".format(fileCounter + 1)

        ## Begin processing each trial in the block
        # Obtain name & path of current data file to process
        currentFileName = blockDataList[fileCounter]
        currentFilePath = os.getcwd() + '\\' + currentFileName

        ## Load TRIAL data to memory, save txt file as currentData array
        currentData = numpy.genfromtxt(currentFilePath, skip_header = 1)
        currentDataSize = currentData.shape

        ## Obtain trial parameters from currentData
        # Obtain trial parameters from current trial data
        # Input = currentData array
        # Output = completed row of trial parameters for trialNumber
        trialParameters[trialCounter, :] = get_trial_parameters(currentData)

        print "(Trial template #{})\n".format(trialParameters[trialCounter][0])

        ## calculate the upper and lower bounds per trial
        ## [0] trialNumber; [1] magnification; [2] targetForce;

        # III: If visual feedback is tied to absolute range around distal force
        lowerBoundGrams = trialParameters[trialCounter][2] - 0.5
        upperBoundGrams = trialParameters[trialCounter][2] + 0.5

        ## Load data from current trial
        # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
        # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
        # [8] bit2; [9] bit3; [10] currentlyInTrial voltage;
        # [11] GS0Force; [12] filteredForce; [13] drift;
        #time = currentData[:,2]
        #GS0_ForceVoltage = currentData[:,3]   # GS0 Force sensor voltage, AI0.0
        GS0_ForceGrams = currentData[:,11]    # GS0 Force sensor (grams)
        HHFM_SolenoidForce = currentData[:,4] # HHFM Solenoid voltage, AI0.1
        HHFM_SensorForce = currentData[:,5]   # HHFM Sensor voltage, AI0.2
        #inTrialVoltage = currentData[:,10]     # inTrial analog voltage, AI0.7
        filteredForce = currentData[:,12]     # filtered by gaussian
        driftForce = currentData[:,13]        # raw - filtered

        ## Begin graphing data
        figFlag = 0
        while figFlag == 0:
            fig = plt.figure(figCounter)
            print "Figure #{}...\n".format(figCounter)
            figureTitle = 'Trial#{}; (Prototype #{}: Mag = {}; '\
                   'Target Force = {})'.format(figCounter,
                   trialParameters[trialCounter][0],
                   trialParameters[trialCounter][1],
                   trialParameters[trialCounter][2])
            fig.suptitle(figureTitle, fontsize = 12, fontweight = 'bold')


            ax1 = plt.subplot(211)
            subplotTitle = 'GS0 Force Sensor'
            ax1.set_title(subplotTitle, fontsize = 10)
            
            ## Plot trial force & filtered force to top subplot
            x = numpy.arange(currentDataSize[0])
            
            #Plot lower and upper bounds, x-axis
            #plt.axhline(y = lowerBoundGrams, color = 'g', ls = '--', lw = 0.50)
            #plt.axhline(y = upperBoundGrams, color = 'r', ls = '--', lw = 0.50)
            targetZone = plt.fill_between(x, y1 = lowerBoundGrams, y2 = upperBoundGrams,
                             alpha = 0.30, color = 'y')
            plt.axhline(y = 0, lw = 1.0, color = 'k')
            
            # plot data
            GS0Line, = plt.plot(x, GS0_ForceGrams, 'k', lw = 1)
            #FilterLine, = plt.plot(x, filteredForce, 'r', lw = 0.75, 
            #                       label = 'Filtered Force')
            #DriftLine, = plt.plot(x, driftForce, 'k', label = 'Drift Force',
            #    lw = 0.75)
                
            # set axis ticks and labels
            plt.axis(ymin = -30, ymax = 30)
            plt.setp(ax1.get_xticklabels(), visible = False)
            plt.xticks(numpy.arange(0,len(x),1000))
            plt.ylabel('Applied Force '+r'$(grams)$', fontsize = 10)            
            plt.grid(axis = 'y', which = 'major')
            
            plt.legend((GS0Line, targetZone), ('GS0-100 Force (raw)', 'Target Window'),
                       loc = 'lower right', fontsize = 8, handlelength = 3)

            ## Graph vertical lines to mark time points in trial
            # Actual begin of trial stimulus 20 pts after start of file
            #plt.axvline(x = 20, color = 'k', ls = '--', lw = 1.25)

            # End of trial marker = 2500 points before end of file
            endStamp = currentDataSize[0] - 2500
            plt.axvline(x = endStamp,
                color = 'r', ls = '--', lw = 1.25)

            # Start of trial happened 12s before end of trial
            # 6s visual feedback, 6s no visual feedback
            plt.axvline(x = endStamp - 12000, color = 'g', ls = '--', lw = 1.25)

            # End of visual feedback, 6s before trial end
            plt.axvline(x = (endStamp - 6000), color = 'r', ls = '--', lw= 1.25)
            plt.text(x = (endStamp - 4800), y = 22, s = '<---(no visual feedback)--->',
                style = 'italic', fontsize = 8)
            
            ## Plot trial HHFM voltages to bottom subplot
            ax2L = plt.subplot(212, sharex = ax1)
            #plt.axis(ymin = 3, ymax = 7)
            subplotTitle = 'HHFM'
            plt.title(subplotTitle, fontsize = 10)

            x = numpy.arange(currentDataSize[0])
            HHFM_out, = plt.plot(x, HHFM_SolenoidForce, 'b', label = 'HHFM Actuator')
            ax2L.set_xlabel('Time '+r'$(sec)$', fontsize = 10)
            ax2L.set_ylabel('Actuator Voltage', fontsize = 10)

            # setup right hand axes for sensor voltage
            ax2R = ax2L.twinx()
            ax2R.set_ylabel('Sensor Voltage', fontsize = 10)
            HHFM_in, = ax2R.plot(x, HHFM_SensorForce, 'm', label = 'HHFM Sensor')
            #plt.axis(ymin = 3, ymax = 7)
            plt.xticks(numpy.arange(0,len(x),1000), 
                       numpy.arange(0,math.ceil(len(x)/1000)+1,1, dtype = numpy.int))
            plt.tick_params(axis = 'x', labelsize = 8)
            plt.grid(axis = 'y', which = 'major')

            plt.legend((HHFM_out, HHFM_in), ('HHFM Actuator', 'HHFM Sensor'),
                loc = 'lower right', fontsize = 8)

            ## Save figure to PDF
            pp.savefig(figCounter)
            plt.close('all')
            figCounter = figCounter + 1
            figFlag = 1

        ## Done graphing, advance to next file/ trial
        fileCounter = fileCounter + 1
        trialCounter = trialCounter + 1

        if fileCounter >= blockDataListSize:
            plt.close('all')
            pp.close()
            print "\nPDF closed.\n"
            print "{} split data graphing completed!\n".format(blockNumberString)
            # newDataFile.close()
            # print "Data file closed.\n"
            doneFlag = 1

    ## Return (16, 4) array
    return trialParameters

def graph_full_data(fullData, timestamps, trialParameters, blockNumberString):
    '''Graph the full raw data with annotations.

    Draw two figures: (1) the full raw data and (2) the reordered data, both
    with annotations.
    '''
    # Extract the reordered list
    orderedBlockDataList = numpy.empty(16, dtype = 'str')
    orderedBlockDataList = reorder_trials(blockNumberString)

    # Overall size of raw data
    fullDataSize = fullData.shape

    ## Open PDF for plots and figures
    pdfFileName = os.getcwd() +'\\'+ blockNumberString + '_FullDataPlots.pdf'
    pp = PdfPages(pdfFileName)            # Pointer to pdf file

    ## Load data into function scope from current trial
    # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
    # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
    # [8] bit2; [9] bit3; [10] currentlyInTrial voltage;
    # [11] GS0Force; [12] filteredForce; [13] drift;
    GS0_ForceGrams = fullData[:,10]

    ## Begin graphing raw data, trials in order as they were presented
    fig = plt.figure(1)
    figureTitle = '{} Raw Data'.format(blockNumberString)
    fig.suptitle(figureTitle, fontsize = 12, fontweight = 'bold')

    ax1 = plt.subplot(111)
    ax1.set_ylabel('Grams')
    ax1.set_xlabel('Time')

    subplotTitle = 'Applied Force (GS0 Force Sensor)'
    plt.title(subplotTitle, fontsize = 10)
    plt.axis(ymin = -30, ymax = 30)

    x = numpy.arange(fullDataSize[0])
    plt.plot(x, GS0_ForceGrams, 'b', label = 'GS0-100 Force')
    plt.legend(loc = 'lower right', fontsize = 'xx-small')

    ## Add color target force/magnification color code
    ## [0] trialNumber; [1] magnification; [2] targetForce;
    timestampCounter = 0
    for trialNumber in range(0, 16):

        startStamp = timestamps[timestampCounter]
        endStamp = timestamps[timestampCounter + 1]

        trialPath = os.getcwd() + '\\' + 'Trial{0:02d}.txt'.format(trialNumber + 1)
        currTrialData = numpy.genfromtxt(trialPath, skip_header = 1)

        ## obtain trialParameters from fullData timestamps
        currParameters = numpy.zeros(shape=(4))
        currParameters = get_trial_parameters(currTrialData)

        ## Assign color background based on distal force
        if abs(currParameters[2]) == 1.0:
            currentColor = 'yellow'
        elif abs(currParameters[2]) == 2.0:
            currentColor = 'green'
        elif abs(currParameters[2]) == 3.0:
            currentColor = 'violet'
        elif abs(currParameters[2]) == 4.0:
            currentColor = 'red'
        else:
            currentColor = 'cyan'

        ## Assign opacity based on magnified/unmagnified
        if currParameters[1] == 1:
            currentAlpha = 1.00
        elif currParameters[1] == 0:
            currentAlpha = 0.30

        plt.axvspan(xmin = startStamp, xmax = endStamp,
                    color = currentColor, alpha = currentAlpha)
        plt.hlines(y = (currParameters[2]),
                   xmin = startStamp, xmax = endStamp, colors = 'k',
                   linestyles = 'dashed', lw = 1.00)

        ## label target forces and magnification status
        if (trialNumber % 2) == 0:
            plt.text(x = timestamps[timestampCounter] + 50, y = 47,
                s = 'TF = {}\nMag = {}'.format(
                    currParameters[2], int(currParameters[1])),
                fontsize = 'xx-small')
        else:
            plt.text(x = timestamps[timestampCounter] + 50, y = 54,
                s = 'TF = {}\nMag = {}'.format(
                    currParameters[2],
                    int(currParameters[1])),
                fontsize = 'xx-small')

        # Advance past last used timestamp pair
        timestampCounter = timestampCounter + 2

    ## Save figure to PDF
    pp.savefig()
    plt.close('all')            # close figure
    print "\nPDF closed.\nFull data graphing 1 completed!\n"

    print "Graphing reordered trials..."
    fig = plt.figure(2)
    figureTitle = '{} Raw Data (Reordered)'.format(blockNumberString)
    fig.suptitle(figureTitle, fontsize = 12, fontweight = 'bold')

    ax = plt.subplot(111)
    ax.set_xlabel('Grams')
    ax.set_ylabel('Time')

    subplotTitle = 'Applied Force (GS0 Force Sensor)'
    plt.title(subplotTitle, fontsize = 10)
    plt.axis(ymin = -30, ymax = 30)

    # Begin graphing raw data, ordered so mag/unmag trials are adjacent
    lastEndStamp = 0
    for trialNumber in range(0,16):

        trialName = orderedBlockDataList[trialNumber]
        trialPath = os.getcwd() + '\\' + trialName

        currTrialData = numpy.genfromtxt(trialPath, skip_header = 1)
        GS0_ForceGramsTrial = currTrialData[:,11]

        startStamp = lastEndStamp + 1
        endStamp = startStamp + len(GS0_ForceGramsTrial)

        x = numpy.arange(startStamp, endStamp)
        plt.plot(x, GS0_ForceGramsTrial, 'b', label = 'GS0-100 Force')

        currParameters = trialParameters[trialNumber,:]

        ## Assign color background based on distal force
        if abs(currParameters[2]) == 1.0:
            currentColor = 'yellow'
        elif abs(currParameters[2]) == 2.0:
            currentColor = 'green'
        elif abs(currParameters[2]) == 3.0:
            currentColor = 'violet'
        elif abs(currParameters[2]) == 4.0:
            currentColor = 'red'
        else:
            currentColor = 'cyan'

        ## Assign opacity based on magnified/unmagnified
        if currParameters[1] == 1:
            currentAlpha = 1.00
        elif currParameters[1] == 0:
            currentAlpha = 0.30

        plt.axvspan(xmin = startStamp, xmax = endStamp,
                    color = currentColor, alpha = currentAlpha)
        plt.hlines(y = (currParameters[2]),
                   xmin = startStamp, xmax = endStamp, colors = 'k',
                   linestyles = 'dashed', lw = 1.00)

        ## label target forces and magnification status
        if (trialNumber % 2) == 0:
            plt.text(x = startStamp + 50, y = 47,
                s = 'TF = {}\nMag = {}'.format(
                    currParameters[2], int(currParameters[1])),
                fontsize = 'xx-small')
        else:
            plt.text(x = startStamp + 50, y = 54,
                s = 'TF = {}\nMag = {}'.format(
                    currParameters[2], int(currParameters[1])),
                fontsize = 'xx-small')

        lastEndStamp = endStamp

    ## Save figure to PDF
    pp.savefig()
    plt.close('all')
    pp.close()
    print "\nPDF closed. \nFull data graphing 2 completed!\n"

    return


def graph_fourier(trialParameters, blockNumberString):
    '''Calculate and graph FFT/ power spectra of raw GS0-100 force data

    Calculates, saves, and plots FFT of raw GS0 Force sensor data. Calculates
    and plots the power spectrum, the integral of the squared coefficients.
    Plots show FFT for f = [1 - 10] Hz, and power for f = [1 - 4] Hz.
    '''

    # 1-D array with ORDERED list of file names
    blockDataList = numpy.empty(16, dtype = 'str')
    blockDataList = reorder_trials(blockNumberString)
    print 'blockDataList = {}\n'.format(blockDataList)

    # Obtain number of files in subjectData directory
    blockDataListSize = len(blockDataList)

    ## open PDF for plots and figures
    pdfFileName = os.getcwd() +'\\'+ blockNumberString + '_FourierPlots.pdf'
    pp = PdfPages(pdfFileName)            #Pointer to pdf file

    ## create empty arrays and flags
    print "Reading trial data to graph FFT...\n"

    # Flag to iterate through all files in block
    doneFlag = 0

    ## Graph and obtain trial parameters for all trials in block
    trialCounter = 0
    fileCounter = 0
    figCounter = 1
    while(doneFlag == 0):
        # Print trial number to screen
        print "Trial #{} FFT graphing...".format(fileCounter + 1)

        ## Begin processing each trial in the block
        # Path to data file currently processing
        currentFileName = blockDataList[fileCounter]
        currentFilePath = os.getcwd() + '\\' + currentFileName

        ## Load TRIAL data to memory, save txt file as currentData array
        currentData = numpy.genfromtxt(currentFilePath, skip_header = 1)
        #currentDataSize = currentData.shape

        # 16 trials, 3 columns each
        # [0] trialNumber [1] mag [2] targetForce [3] perceivedForce
        print "(Trial template #{})\n".format(trialParameters[trialCounter,0])

        ## Load data from current trial
        # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
        # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
        # [8] bit2; [9] bit3; [10] currentlyInTrial voltage;
        # [11] GS0Force; [12] filteredForce; [13] drift;
        #time = currentData[:,2]
        #GS0_ForceVoltage = currentData[:,3]   # GS0 Force sensor voltage, AI0.0
        GS0_ForceGrams = currentData[:,11]    # GS0 Force sensor (grams)
        #HHFM_SolenoidForce = currentData[:,4] # HHFM Solenoid voltage, AI0.1
        #HHFM_SensorForce = currentData[:,5]   # HHFM Sensor voltage, AI0.2
        #inTrialVoltage = currentData[:,10]     # inTrial analog voltage, AI0.7
        #filteredForce = currentData[:,12]     # filtered by gaussian
        #driftForce = currentData[:,13]        # raw - filtered

        ## Begin graphing data
        # Set up output figures
        fig = plt.figure(figCounter)
        print "Figure #{}...\n".format(figCounter)
        figureTitle = 'Trial#{}; (Prototype #{}: Mag = {};'\
                'Target Force = {})'.format(figCounter,
                    trialParameters[trialCounter, 0],
                    trialParameters[trialCounter, 1],
                    trialParameters[trialCounter, 2])
        fig.suptitle(figureTitle, fontsize = 12, fontweight = 'bold')

        ax1L = plt.subplot(211)
        subplotTitle = 'Trial #{}; FFT Applied Force '\
            '(GS0 Force Sensor)'.format(figCounter)
        ax1L.set_title(subplotTitle, fontsize = 8)

        ## Define start and stop points for trim, calculate FFT & freq
        # split raw data trims raw data 2500 points after the end of trial sig
        endStamp = len(GS0_ForceGrams) - 2500
        startStamp = endStamp - 12000

        rawForceVisFD = GS0_ForceGrams[startStamp : endStamp - 6000]
        rawForceNoFD = GS0_ForceGrams[endStamp - 6000 : endStamp]

        fftRawVisFD = numpy.abs(numpy.fft.fft(rawForceVisFD, axis = 0))
        fftRawNoFD = numpy.abs(numpy.fft.fft(rawForceNoFD, axis = 0))

        # timestep in units of seconds to obtain cycles per second
        timestep = 0.001
        freqVisFD = numpy.fft.fftfreq(rawForceVisFD.size, d = timestep)
        freqNoFD = numpy.fft.fftfreq(rawForceNoFD.size, d = timestep)

        fftVisFD, = plt.plot(freqVisFD, fftRawVisFD, 'b',
            label = 'FFT Raw Force (visual feedback)', lw = 1.50)

        # Set up axes now that figure has been plotted
        ax1L.set_xlim(left = 1, right = 10)
        ax1L.set_ylim(bottom = 0, top = 10000)
        plt.tick_params(axis = 'both', which = 'major', labelsize = 8)

        y_label = r'$\mathcal{F}(Force)$'
        ax1L.set_ylabel(y_label, fontsize = 8)
        #ax1L.set_ylabel('FFT(Force_Grams), vis FD', fontsize = 8)
        ax1L.set_xlabel('Frequency (Hz)', fontsize = 8)

        # Second set of axes for no FD FFT
        ax1R = ax1L.twinx()

        fftNoFD, = ax1R.plot(freqNoFD, fftRawNoFD, 'r',
            label = 'FFT Raw Force (no feedback)', lw = 0.75)

        #ax1R.set_ylabel('FFT(Force_Grams), no FD', fontsize = 8)
        #y_label = r'$\mathcal{F}(Force), no FD$'
        #ax1R.set_ylabel(y_label, fontsize = 8)
        ax1R.set_xlim(left = 1, right = 10)
        ax1R.set_ylim(bottom = 0, top = 10000)
        #plt.tick_params(axis = 'both', which = 'major', labelsize = 8)
        plt.setp(ax1R.get_yticklabels(), visible = False)

        # Set up legend for Plot1
        plt.legend((fftVisFD, fftNoFD), ('FFT Raw Force (visual feedback)',
                'FFT Raw Force (no feedback)'),
            loc = 'best', fontsize = 6)

        # Save FFT as text file
        temp1 = numpy.column_stack((freqVisFD, fftRawVisFD))
        temp2 = numpy.column_stack((freqNoFD, fftRawNoFD))
        temp = numpy.column_stack((temp1, temp2))

        fileName = os.getcwd() +'\\'+ blockNumberString + '_{}_FourierCoeff.txt'.format(fileCounter + 1)

        numpy.savetxt(fileName, temp, fmt = '%.4e',
             header = 'Freq VisFD\t fftRawVisFD\t Freq NoFD\t fftRawNoFD',
             delimiter = '\t')

        # take out all negative frequencies
        conditions = [freqVisFD >= 0]
        choices = [freqVisFD]
        posFreq = numpy.trim_zeros(numpy.select(conditions,choices), 'b')

        # upper bound for indices
        indices = numpy.array(numpy.nonzero(posFreq <= 4))
        indices = numpy.transpose(indices)

        #print "indices - ", indices
        #print "freq - {}".format(posFreq[indices])
        #print "coefficients - \nvisFD: {}\n noFD: {}\n".format(
        #    fftRawVisFD[indices], fftRawNoFD[indices])

        # set up integral
        index = 0
        pwrVisFD = numpy.zeros(indices.shape)
        pwrNoFD = numpy.zeros(indices.shape)
        pwrDiff = numpy.zeros(indices.shape)
        sumVisFD = numpy.zeros(indices.shape)
        sumNoFD = numpy.zeros(indices.shape)
        sumDiff = numpy.zeros(indices.shape)
        diff = fftRawVisFD - fftRawNoFD

        # multiply by conjugate because int{ |F(w)|^2 } = int{ |f(t)|^2 }
        # by perseval's formula. Returns measure of "energy" in signal
        # from index 6 -> 25 because those match with posFreq = 1 -> 4 Hz
        for i in range(6, 25):
            pwrVisFD[index] = (fftRawVisFD[i] * numpy.conjugate(fftRawVisFD[i]))
            pwrNoFD[index] = (fftRawNoFD[i] * numpy.conjugate(fftRawNoFD[i]))
            pwrDiff[index] = diff[i] * numpy.conjugate(diff[i])
            index += 1

        # calculate integral for pwr from 1 - 4 Hz
        for index in range(len(pwrVisFD)):
            if index == 0:
                sumVisFD[index] = pwrVisFD[index]
                sumNoFD[index] = pwrNoFD[index]
                sumDiff[index] = pwrDiff[index]
            else:
                sumVisFD[index] += pwrVisFD[index] + sumVisFD[index - 1]
                sumNoFD[index] += pwrNoFD[index] + sumNoFD[index - 1]
                sumDiff[index] += pwrDiff[index] + sumDiff[index - 1]

        #print sumVisFD
        #print sumNoFD
        #print sumDiff

        #cumulative = max(numpy.amax(sumVisFD), numpy.amax(sumNoFD))

        ## Second plot, for integral of FFT
        ax2L = plt.subplot(212)

        plt.xlabel('Frequency (Hz)', fontsize = 8)
        y_label = r'$\int\/\/|\mathcal{F}(Force)|^2 $'
        ax2L.set_ylabel(y_label, fontsize = 8)

        intFFTVisFD, = ax2L.plot(posFreq[indices], sumVisFD, 'b',
            label = 'Cumulative Power Spectrum (vis FD)', lw = 1.50)
        ax2L.set_xlim(left = 1, right = 4)
        ax2L.set_ylim(bottom = 0, top = 1e8)
        plt.tick_params(axis = 'both', which = 'major', labelsize = 8)

        ax2R = ax2L.twinx()

        intFFTNoFD, = ax2R.plot(posFreq[indices], sumNoFD, 'r',
            label = 'Cumulative Power Spectrum (no FD)', lw = 0.75)

        intFFTDiff, = ax2R.plot(posFreq[indices], sumDiff, 'k',
            label = 'Cumulative Power Spectrum (diff)', lw = 1.00)

        ax2R.set_xlim(left = 1, right = 4)
        ax2R.set_ylim(bottom = 0, top = 1e8)
        plt.setp(ax2R.get_yticklabels(), visible = False)
        plt.tick_params(axis = 'both', which = 'major', labelsize = 8)

        plt.legend((intFFTVisFD, intFFTNoFD, intFFTDiff),
            ('Cumulative Power Spectrum (visual feedback)',
             'Cumulative Power Spectrum (no feedback)',
             'Cumulative Power Spectrum (diff)'),
            loc = 'best', fontsize = 6)

        #plt.tight_layout()

        ## Save figure to PDF
        pp.savefig(figCounter)
        plt.close('all')
        figCounter = figCounter + 1

        ## Done graphing, advance to next file/ trial
        fileCounter = fileCounter + 1
        trialCounter = trialCounter + 1

        if fileCounter >= blockDataListSize:
            plt.close('all')
            pp.close()
            print "\nPDF closed.\n"
            print "{} split data graphing completed!\n".format(blockNumberString)
            # newDataFile.close()
            # print "Data file closed.\n"
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


def getBinaryFromVoltage(voltage):
    '''Infer binary state from float voltage data.

    Returns 1 if voltage > 1.00'''
    if voltage > 1.00:
        # one bit is max @ ~1.50V
        result = 1.
    else:
        result = 0.

    return result

def filter_forces(fullDataArray):
    "Precondition: raw force traces from 1 trial. Postcondition: filtered forces."

    # 21 element gaussian, sigma = 4
    #gaussian = numpy.array([0.000272337, 0.00089296, 0.002583865, 0.00659813,
    #                        0.014869116, 0.029570767, 0.051898313, 0.080381679,
    #                        0.109868729, 0.132526984, 0.14107424, 0.132526984,
    #                        0.109868729, 0.080381679, 0.051898313, 0.029570767,
    #                        0.014869116, 0.00659813, 0.002583865, 0.00089296,
    #                        0.000272337])

    # 201 element gaussian file, sigma = 4:
    #fileName = rawDataDirectory + '\\' + '201element_normGaussian.txt'
    #gaussian = numpy.genfromtxt(fileName, skip_header = 1)

    # 501 element gaussian file, sigma = 200:
    fileName = rawDataDirectory + '\\' + '501element_normGaussian.txt'
    gaussian = numpy.genfromtxt(fileName, skip_header = 1)

    #filteredBlock = numpy.
    GS0_ForceGrams = fullDataArray[:,10]
    size = GS0_ForceGrams.shape
    filteredForces = numpy.zeros(shape = size)

    # Loop for 21 element gaussian
    #for timestamp in range(11, size[0] - 10):
    #    rawChunk = GS0_ForceGrams[timestamp - 11 : timestamp + 10]
    #
    #    pointAverage = numpy.convolve(rawChunk, gaussian, mode = 'valid')
    #
    #    filteredForces[timestamp] = pointAverage

    # Loop for 201 element gaussian
    #for timestamp in range(101, size[0] - 100):
    #    rawChunk = GS0_ForceGrams[timestamp - 101 : timestamp + 100]
    #
    #    pointAverage = numpy.convolve(rawChunk, gaussian[:,1], mode = 'valid')
    #
    #    filteredForces[timestamp] = pointAverage

    # Loop for 501 element gaussian
    for timestamp in range(251, size[0] - 250):
        rawChunk = GS0_ForceGrams[timestamp - 251 : timestamp + 250]

        pointAverage = numpy.convolve(rawChunk, gaussian[:,1], mode = 'valid')

        filteredForces[timestamp] = pointAverage

    drift = GS0_ForceGrams - filteredForces

    return filteredForces, drift

###############################################################################
##
## MAIN : Process Data
##
###############################################################################

# Call time to set as start time
start = time.time()

print "Loading raw data folders... \n"

os.chdir(rawDataDirectory)

print "Accessing subject #{} data...\n".format(subjectNumber)

# Change to directory with subject data files
subjectDirectory = rawDataDirectory + '\\' + subjectNumber
os.chdir(subjectDirectory)

# 1D array with list of all raw data file names
# Block data 1-5 first, then thresholds
subjectDataList = os.listdir(os.getcwd())
subjectDataListSize = len(subjectDataList)

print "Extracting trial parameters and splitting block data...\n"

## Loop through all five blocks
for i in range(0, 5):
    ## Create block folder in HHFM_Data/subjectNumber/
    newBlockFolderName = '{}_Block{}'.format(subjectNumber, i + 1)
    newBlockDir = os.getcwd() + '\\' + newBlockFolderName
    os.mkdir(newBlockDir)
    os.chdir(newBlockDir)

    ## Currently in HHFM_Data/subjectNumber/Blockx/
    # Timestamps = numpy.zeros(shape = ( 33, 1 ))
    # Init/ reset timestamps array
    subjectDataFullPath = subjectDirectory + '\\' + subjectDataList[i]
    initFullData = numpy.genfromtxt(subjectDataFullPath, skip_header = 1,
        skip_footer = 1)
    initFullDataSize = initFullData.shape

    ## Create fullData matrix, with extra column for GS0 Force calculation
    # Creates text files with split data
    # Input = initial full data, block number (i + 1)
    # Output = fullData w/ GS0Force, timestampArray
    #fullData = numpy.zeros(shape = (initFullDataSize[0], initFullDataSize[1] + 1))
    fullData, timestamps = split_trials(initFullData, i + 1)

    timestampFilePath = os.getcwd() + '\\' + newBlockFolderName + \
        '_timestamps.txt'
    numpy.savetxt(timestampFilePath, timestamps, fmt = '%d')
    print "Block{} timestamps saved.\n".format(i + 1)

    ## Process split trials and get list of trial parameters for the current block
    # Creates data graphs into PDF format
    # Input = PDF save path, timestampArray, blockNameString
    print "Block{} trial processing...\n".format(i + 1)
    savePath = os.getcwd()
    trialParameters = process_split_trials(savePath, timestamps, newBlockFolderName)

    # Save obtained trialParameters array to text file
    trialParametersFilePath = os.getcwd() + '\\' + newBlockFolderName + \
        '_trialParameters.txt'
    numpy.savetxt(trialParametersFilePath, trialParameters, fmt = '%.2f')

    ## Summary graph of full raw data with annotations
    # Input = arrays: fullData, timestamps, trialParameters, blockNameString
    # Output = PDF with fullData graphs
    graph_full_data(fullData, timestamps, trialParameters, newBlockFolderName)

    ## Graph fourier transforms
    graph_fourier(trialParameters, newBlockFolderName)

    print "Block{} done!\n".format(i + 1)

    ## Go back to HHFM_Data/subjectNumber
    os.chdir(subjectDirectory)


## Loop through remaining threshold files
#thresholdCounter = 1
#for i in range( 5 , subjectDataListSize ):

#    thresholdDataPath = subjectDirectory + '\\' + subjectDataList[ i ]
#    thresholdData = numpy.genfromtxt( thresholdDataPath , skip_header = 1 , skip_footer = 1 )

#    graphThreshold( thresholdData , thresholdCounter )
#    print "Threshold{} graphing completed!\n".format(thresholdCounter)
#    thresholdCounter = thresholdCounter + 1


print "Processing Completed!\nTime elapsed: {0:03f} seconds".format(time.time() - start)
