################################################################################
##  Randy Lee
##  VIA Lab
##  Department of Bioengineering
##  University of Pittsburgh
##
##  Minimum Contact Raw Data Processing
##
##  Function: Separate and plot trial data for blocks of low force exp trials
##
##  Analysis procedure:
##    (1) MinContact_RawDataProcessing to split raw data and graph time
##        series and Fourier/Power spectra
##    (2) MinContact_DataStatistics to obtain mean force and std for each
##        experimental condition, and graph dependent variable analysis of
##        spectral power
##
##   NI-6009 Analog input:
##  - AI0.0 : GS0-100 force sensor; from ext. hardware (0 - 8V)
##  - AI0.1 : HHFM solenoid voltage; from HHFM hardware (-10V - 10V)
##  - AI0.2 : HHFM force sensor voltage; from HHFM hardware (0 - 10V)
##  - AI0.3 : Trial number bit 0; from exp. ADUC7026 P0.0
##  - AI0.4 : Trial number bit 1; from exp. ADUC7026 P0.1
##  - AI0.5 : Trial number bit 2; from exp. ADUC7026 P0.2
##  - AI0.6 : N/C
##  - AI0.7 : Currently in trial; from exp. ADUC7026
##            (i.e. 0 = before/after trial; 1 = currently in trial)
##                  or
##            Threshold felt; from exp. ADUC7026
##
##  Last updated: 20 October 2016
################################################################################

import os
import math
import time
import numpy
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator

subjectNumber = 'S2007' # Only need to change subjectNumber to process all data

rawDataDirectory = 'C:\\Users\\Randy Lee\\Documents\\VIA Lab\\HHFM Data\\MinimumContact_rawData'
#rawDataDirectory = 'D:\\Randy Lee\\Documents\\VIA Lab\\HHFM Data\\MinimumContact_rawData'

# font sizes for axis labels, legends, and other text
label_fontSize = 16
legend_fontSize = 14
text_fontSize = 12
tick_fontSize = 16

# set tick labels to Times New Roman font to match LaTeX axis labels
ticks_font = matplotlib.font_manager.FontProperties(family = 'Times New Roman',
                                                    style = 'normal',
                                                    size = tick_fontSize)

################################################################################
##
## Local script functions
##
################################################################################

def split_trials(rawData, blockNumber):
    '''Create txt files for individual trials & timestamps using raw block data.

    Uses 2D array obtained from raw data files and splits using voltage on
    currentlyInTrial pin (col 9). Copies section from raw data and includes
    conversion to force from GS0 sensor, filtered force, and drift data.
    Trials are unordered -- #1-16 are in the order as they were presented in
    the block.
    
    Return raw data w/ extra data columns and trial timestamps
    '''

    ## Generate empty data arrays and flags
    rawDataSize = rawData.shape
    
    # at most 12 timestamps (6 trials w/ 2 stamps each)    
    timestampArray = numpy.zeros(shape = (13, 1))  
    GS0Force = numpy.zeros(shape = (rawDataSize[0],1))

    ## Scan through full raw data set to look for start/end of trials
    timestampCounter = 0
    for i in range(0, rawDataSize[0]):
        ## Convert GS0 voltage to grams force
        # Abs value because voltage can dip to negative small numbers
        # Calibration = (20g/200mV) to quantify force (cal factor for push)
        # Negative to match "direction" of target, i.e. neg = push, pos = pull
        currentInTrialVoltage = abs(numpy.float(rawData[i, 9]))
        GS0Force[i] = (numpy.float(rawData[i][2])) * -(20/0.200)

        ## Identify timestamps where trials start/stop
        if(i == 0):
            # To prevent false alarm on first data point
            lastInTrialVoltage = currentInTrialVoltage

        inTrialDiff = lastInTrialVoltage - currentInTrialVoltage

        if(inTrialDiff >= 2):
            ## ex. 3.3V --> 0; trial just started
            # Save as 500 ms before actual change, to capture all GPIO ID pins
            timestampArray[timestampCounter] = (i - 500)
            print "startTime = {}\t".format(timestampArray[timestampCounter])
            timestampCounter += 1
        elif(inTrialDiff <= -2):
            ## ex. 0V --> 3.3; trial just ended
            # Save as 2500 timestamps after actual end
            timestampArray[timestampCounter] = (i + 2500)
            print "endTime = {}\n".format(timestampArray[timestampCounter])
            timestampCounter += 1

        lastInTrialVoltage = currentInTrialVoltage
        
        if timestampCounter > 12:
            # stop splitting once we read all 12 trials
            break

    # Add GS0Force (units = grams) as column [10]
    rawData = numpy.column_stack((rawData, GS0Force))

    ## Filter forces with 201 & 501 element gaussian (Takes full raw data array)
    print "Performing gaussian filtering..."
    filteredForces201, filteredForces501 = filter_forces(rawData)
    print "Gaussian filtering done.\n"

    # Filtered force (grams) as column [11] & [12]
    # Get new data size, with these 2 new data columns
    rawData = numpy.column_stack((rawData, filteredForces201, filteredForces501))
    rawDataSize = rawData.shape

    print "Writing individual trial files...\n"

    # Reset timestamp location counter for each new block
    timestampCounter = 0
    ## Create & save individual trial files for this block
    for fileCounter in range(1, 7):
        ## 6 trials
        ## start going through timestamp holder in pairs
        startTime = numpy.int(timestampArray[timestampCounter])
        endTime = numpy.int(timestampArray[timestampCounter + 1])

        # rows = newDataSize
        # columns = 13
        newDataSize = endTime - startTime
        newData = numpy.zeros(shape = (newDataSize, 13))

        ## Begin writing the data files
        trialDataFileName = newBlockDir + \
                            '\\Trial{0:02d}.txt'.format(fileCounter)
        for newDataLocation in range(0, newDataSize):
            # Raw data location
            rawDataLocation = startTime + newDataLocation

            ## Copy data into new trial file
            # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
            # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
            # [8] bit2; [9] currentlyInTrial voltage;
            # [10] GS0Force; [11] filteredForce_201; [12] filteredForce_501
            newData[newDataLocation][0] = blockNumber
            newData[newDataLocation][1] = fileCounter
            newData[newDataLocation][2] = numpy.int(newDataLocation)        # timestamp
            newData[newDataLocation][3] = numpy.float(rawData[rawDataLocation][2])   # GS0-100 Force Voltage
            newData[newDataLocation][4] = numpy.float(rawData[rawDataLocation][3])   # HHFM Solenoid Voltage
            newData[newDataLocation][5] = numpy.float(rawData[rawDataLocation][4])   # HHFM Force Sensor
            newData[newDataLocation][6] = numpy.float(rawData[rawDataLocation][5])   # trialNumber bit0
            newData[newDataLocation][7] = numpy.float(rawData[rawDataLocation][6])   # trialNumber bit1
            newData[newDataLocation][8] = numpy.float(rawData[rawDataLocation][7])   # trialNumber bit2
            newData[newDataLocation][9] = numpy.float(rawData[rawDataLocation][9])   # currentlyInTrial
            newData[newDataLocation][10] = numpy.float(rawData[rawDataLocation][10]) # calibration = (100g/8V) for GS0-100
            newData[newDataLocation][11] = numpy.float(rawData[rawDataLocation][11]) # Force (grams) filtered with 201 element gaussian
            newData[newDataLocation][12] = numpy.float(rawData[rawDataLocation][12]) # Force (grams) filtered with 501 element gaussian

        numpy.savetxt(trialDataFileName, newData, fmt='%3.6f',
            header = 'Block Number\t Trial Number\t Timestamp\t'
                     'GS0-100 Voltage\t HHFM Actuator Voltage\t'
                     'HHFM Sensor Voltage\t trialNum bit0\t trialNum bit1\t'
                     'trialNum bit2\t currentlyInTrial\t'
                     'GS0-100 Force (grams)\t Filtered force (grams) 201\t'
                     'Filtered Force (grams) 501',
            delimiter = '\t')

        # Advance timestampCounter to next start/end pair
        timestampCounter += 2

    print "Block{} individual trials complete.\n".format(blockNumber)
    return rawData, timestampArray

def get_trial_parameters(trialData):
    '''Return parameter array from binary representation of trial number

    Takes trial-segmented data to extract the prototype ID and return an array
    describing the trial parameters: prototype number, magnification state,
    distal force, and perceived force.
    '''

    ## [0] trial_ID [1] direction [2] mag [3] perceived_target
    parameters = numpy.zeros(shape = 4, dtype='float')

    ## Extract binary encoding at 200 timestamps after GPIO ID pins already set,
    ## (before inTrial pulses)
    ## 3.3V = 1; 0V = 0
    bit0 = getBinaryFromVoltage(trialData[200][6])    # AI0.3
    bit1 = getBinaryFromVoltage(trialData[200][7])    # AI0.4
    bit2 = getBinaryFromVoltage(trialData[200][8])    # AI0.5

    ## Concatenate bits to obtain trial number
    trialNumber = (bit0 * (1.)) + (bit1 * (2.)) + (bit2 * (4.))

    # Add 1 to trialNumber since binary rep is 1 less
    parameters[0] = trialNumber + 1
    print "bit0 = {}, bit1 = {}, bit2 = {}".format(bit0, bit1, bit2)
    print "--> trialNumber = {}\n".format(parameters[0])

    targetForce = 15.
    ## Trials grouped by distal force, neg = push, pos = pull
    ## [0] trial_ID [1] direction [2] mag [3] perceived_target
    if(parameters[0] == 1):
        parameters[1] = -1.
        parameters[2] = 0.
    elif(parameters[0] == 2):
        parameters[1] = -1.
        parameters[2] = 2.
    elif(parameters[0] == 3):
        parameters[1] = -1.
        parameters[2] = 4.
    elif(parameters[0] == 4):
        parameters[1] = 1.
        parameters[2] = 0.
    elif(parameters[0] == 5):
        parameters[1] = 1.
        parameters[2] = 2.
    elif(parameters[0] == 6):
        parameters[1] = 1.
        parameters[2] = 4.
    
    # perceived target force = TF * gain
    magnification = parameters[2]
    parameters[3] = (magnification + 1) * parameters[1] * targetForce

    return parameters

def process_split_trials(pdfSaveFolderPath, timestampArray, blockNumberString):
    '''Graph previously split trials using array of start/stop timestamps
    
    '''

    print "Reading trial data...\n"
    
    ## currently in HHFM_Data/subjectNumber/Blockx/
    ## make list of block directory
    # 1-D array with list of all file names (split trials)
    #blockDataList = os.listdir(os.getcwd())

    # 1-D array with ORDERED list of split file names
    print "blockNumberString = {}\n".format(blockNumberString[2:])
    blockDataList = numpy.empty(6, dtype='str')
    blockDataList = reorder_trials(blockNumberString)
    print "blockDataList = {}\n".format(blockDataList)

    # Obtain number of files in subjectData directory
    blockDataListSize = len(blockDataList)

    ## open PDF for plots and figures
    #pdfFileName = os.getcwd() +'\\'+ blockNumberString + '_DataPlots.pdf'
    pdfFileName = pdfSaveFolderPath + '\\' + blockNumberString + '_DataPlots.pdf'
    pp = PdfPages(pdfFileName)            #Pointer to pdf file

    # 6 trials, 3 columns each
    # [0] trialNumber [1] mag [2] targetForce [3] perceivedTarget
    trialParameters = numpy.zeros(shape = (6,4))

    ## Graph and obtain trial parameters for all trials in block
    # fileCounter = 0 with ordered list because only trials in list
    allDone = 0
    trialCounter = 0
    fileCounter = 0
    figCounter = 1
    while(allDone == 0):
        # Print trial number to screen
        print "Trial #{} graphing...".format(fileCounter + 1)

        ## Begin processing each trial in the block
        # Obtain name & path of current data file to process
        currentFileName = blockDataList[fileCounter]
        currentFilePath = os.getcwd() + '\\' + currentFileName

        ## Load TRIAL data to memory, save txt file as currentData array
        currentData = numpy.genfromtxt(currentFilePath, skip_header=1)
        currentDataSize = currentData.shape

        ## Obtain trial parameters from currentData
        # Obtain trial parameters from current trial data
        # Input = currentData array
        # Output = completed row of trial parameters for trialNumber
        trialParameters[trialCounter, :] = get_trial_parameters(currentData)

        print "(Trial template #{})\n".format(trialParameters[trialCounter][0])

        ## Trials grouped by direction, neg = push, pos = pull
        ## [0] trial_ID [1] direction [2] mag [3] perceived_target
        # get signed target force; neg = push, pos = pull
        targetForce = trialParameters[trialCounter][1] * 15.

        ## Load data from current trial
        # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
        # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
        # [8] bit2; [9] currentlyInTrial voltage;
        # [10] GS0Force; [11] filteredForce_201; [12] filteredForce_501
        GS0_ForceGrams = currentData[:,10]      # GS0 Force sensor (grams)
        HHFM_SolenoidForce = currentData[:,4]   # HHFM Solenoid voltage, AI0.1
        HHFM_SensorForce = currentData[:,5]     # HHFM Sensor voltage, AI0.2
        filteredForce_201 = currentData[:,11]   # filtered by 201 element gaussian
        filteredForce_501 = currentData[:,12]   # filtered by 501 element gaussian
        
        saturation_timestamps = saturation_check(HHFM_SolenoidForce)
        fileName = os.getcwd() +'\\'+ blockNumberString + \
                   '_{}_saturationTimestamps.txt'.format(fileCounter + 1)
        numpy.savetxt(fileName, saturation_timestamps, fmt = '%d', 
                      header = '#saturation timestamp', delimiter = '\t')

        ## Begin graphing data
        figFlag = 0
        while figFlag == 0:
            fig = plt.figure(figCounter)
            print "Figure #{}...\n".format(figCounter)
            
            if trialParameters[trialCounter,1] == -1:
                direction = 'Push'
            elif trialParameters[trialCounter,1] == 1:
                direction = 'Pull'
            figureTitle = 'Trial#{} - {} (Direction: {}({}); Mag = {}x)'.format(
                   figCounter, currentFileName, direction, 
                   int(trialParameters[trialCounter, 1]), 
                   int(trialParameters[trialCounter, 2]))
            fig.suptitle(figureTitle, fontsize=8, fontweight='bold')

            ax1 = plt.subplot(211)
            subplotTitle = r'$\mathrm{GS0}$'+'-'+r'$\mathrm{100\ Force\ Sensor}$'
            ax1.set_title(subplotTitle, fontsize=label_fontSize)
            
            ## Plot trial force & filtered force to top subplot
            x = numpy.arange(currentDataSize[0])

            # set axis ticks and labels
            plt.axis(ymin=-45, ymax=45)
            plt.setp(ax1.get_xticklabels(), visible=False)
            plt.xticks(numpy.arange(0,len(x),1000))
            plt.ylabel(r'$\mathrm{Applied\ Force}\ \mathit{(grams)}$',
                       fontsize=label_fontSize+2)            
            plt.grid(axis='y', which='major') 

            # update tick labels to correct font
            for label in ax1.get_xticklabels():
                label.set_fontproperties(ticks_font)
    
            for label in ax1.get_yticklabels():
                label.set_fontproperties(ticks_font)
            
            # Plot target force and x axis
            targetLine = plt.axhline(y=targetForce, color='grey', ls='--', 
                                     lw=0.75, label=r'$\mathrm{Target\ Force}$')
            plt.axhline(y=0, lw=1.0, color='k')
            
            # Plot data
            GS0Line, = plt.plot(x, GS0_ForceGrams, 'k', lw=1)
            #FilterLine201, = plt.plot(x, filteredForce_201, 'r', lw=0.75)
            #FilterLine501, = plt.plot(x, filteredForce_501, 'c--', lw=0.75)
            
            # Plot saturation sections
            saturation = numpy.arange(-30,30)
            if saturation_timestamps[0] == -1:
                print "no saturation to graph!"
                pass
            else:
                for i in range(len(saturation_timestamps)):
                    left = saturation_timestamps[i]
                    right = left + 20
                    plt.fill_betweenx(saturation, x1=left, x2=right, 
                                      facecolor='y', alpha=0.1)
            
#            plt.legend((GS0Line, FilterLine201, FilterLine501), 
#                       ('GS0-100 Force (raw)', 
#                        'Filtered Force (201 element gaussian, sigma=100)',
#                        'Filtered Force (501 element gaussian, sigma=200)'),
#                       loc='lower right', fontsize=6)


            ## Graph vertical lines to mark time points in trial
            # Start of trial marker = 500 points after start of file
            # End of trial marker = 2500 points before end of file
            endStamp = currentDataSize[0] - 2500
            startStamp = 500
            startLine = plt.axvline(x=startStamp, color='g', ls='--', lw=1.25)
            endLine = plt.axvline(x=endStamp, color='r', ls='--', lw=1.25)

            # End of visual feedback, 6s before trial end
            visLine = plt.axvline(x=endStamp - 6000, color='r', ls='--', lw=1.25)
            plt.text(x=endStamp - 5300, y=32, 
                     s=r'$\longleftarrow \mathit{(no\ visual\ feedback)} \longrightarrow$',
                     fontsize=text_fontSize)
            
            plt.legend((GS0Line, targetLine),
                       (r'$\mathrm{GS0}$'+'-'+r'$\mathrm{100\ Force\ (raw)}$',
                        r'$\mathrm{Target\ Force}$'), handlelength=3,
                       loc='lower right', fontsize=legend_fontSize)

            ## Plot trial HHFM voltages to bottom subplot
            ax2L = plt.subplot(212, sharex = ax1)
            plt.axis(ymin = 2, ymax = 7)
            subplotTitle = r'$\mathrm{HHFM}$'
            plt.title(subplotTitle, fontsize=label_fontSize)

            x = numpy.arange(currentDataSize[0])
            HHFM_out, = plt.plot(x, HHFM_SolenoidForce, 'b', 
                                 label=r'$\mathrm{HHFM\ Actuator}$')
            
            # Plot saturation sections
            if saturation_timestamps[0] == -1:
                print "no saturation to graph!"
                pass
            else:
                for i in range(len(saturation_timestamps)):
                    left = saturation_timestamps[i]
                    right = left + 20
                    plt.fill_betweenx(saturation, x1=left, x2=right, facecolor='y', 
                                     alpha=0.1)

            ax2L.set_xlabel(r'$\mathrm{Time}\ \mathit{(sec)}$', fontsize=label_fontSize+3)
            ax2L.set_ylabel(r'$\mathrm{Actuator\ Voltage}$', fontsize=label_fontSize+2)
            
            # update tick labels to correct font
            for label in ax2L.get_xticklabels():
                label.set_fontproperties(ticks_font)
    
            for label in ax2L.get_yticklabels():
                label.set_fontproperties(ticks_font)            

            # setup right hand axes for sensor voltage
            ax2R = ax2L.twinx()
            ax2R.set_ylabel(r'$\mathrm{Sensor\ Voltage}$', fontsize=label_fontSize+2)
            HHFM_in, = ax2R.plot(x, HHFM_SensorForce, 'm', 
                                 label=r'$\mathrm{HHFM\ Sensor}$')
                                     
            plt.axis(ymin = 5, ymax = 6)
            plt.xticks(numpy.arange(0,len(x),1000), 
                       numpy.arange(0,math.ceil(len(x)/1000)+1,1, 
                                    dtype = numpy.int))
            plt.tick_params(axis='x', labelsize=tick_fontSize)
            plt.grid(axis='y', which='major')
            
            # update tick labels to correct font
            for label in ax2R.get_xticklabels():
                label.set_fontproperties(ticks_font)
    
            for label in ax2R.get_yticklabels():
                label.set_fontproperties(ticks_font)            

            plt.legend((HHFM_out, HHFM_in), 
                       (r'$\mathrm{HHFM\ Actuator}$',
                        r'$\mathrm{HHFM\ Sensor}$'),
                       loc='lower right', fontsize=legend_fontSize)

            ## Save figure to PDF
            pp.savefig(figCounter)
            plt.close('all')
            figCounter += 1
            figFlag = 1

        ## Done graphing, advance to next file/ trial
        fileCounter += 1
        trialCounter += 1

        if fileCounter >= blockDataListSize:
            plt.close('all')
            pp.close()
            print "\nPDF closed.\n"
            print "{} split data graphing completed!\n".format(blockNumberString)
            # newDataFile.close()
            # print "Data file closed.\n"
            allDone = 1

    ## Return (6, 4) array
    return trialParameters

def graph_full_data(fullData, timestamps, trialParameters, blockNumberString):
    '''Graph the full raw data with annotations.

    Draw two figures: (1) the full raw data and (2) the reordered data, both
    with annotations.
    '''
    # Extract the reordered list
    orderedBlockDataList = numpy.empty(6, dtype = 'str')
    orderedBlockDataList = reorder_trials(blockNumberString)

    # Overall size of raw data
    fullDataSize = fullData.shape

    ## Open PDF for plots and figures
    pdfFileName = os.getcwd() +'\\'+ blockNumberString + '_FullDataPlots.pdf'
    pp = PdfPages(pdfFileName)            # Pointer to pdf file

    ## Load data into function scope from current trial
    # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
    # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
    # [8] bit2; [9] currentlyInTrial voltage;
    # [10] GS0Force; [11] filteredForce;
    GS0_ForceGrams = fullData[:,10]

    ## Begin graphing raw data, trials in order as they were presented
    fig = plt.figure(1)
    figureTitle = '{} Raw Data'.format(blockNumberString)
    fig.suptitle(figureTitle, fontsize=12, fontweight='bold')

    ax1 = plt.subplot(111)
    ax1.set_ylabel('Grams')
    ax1.set_xlabel('Time')

    subplotTitle = 'Applied Force (GS0 Force Sensor)'
    plt.title(subplotTitle, fontsize=10)
    plt.axis(ymin=-45, ymax=45)

    x = numpy.arange(fullDataSize[0])
    plt.plot(x, GS0_ForceGrams, 'b', label='GS0-100 Force')
    plt.legend(loc='lower right', fontsize='xx-small')

    ## Add color target force/magnification color code
    ## [0] trial_ID [1] direction [2] mag [3] perceived_target
    timestampCounter = 0
    for trialNumber in range(0, 6):
        startStamp = timestamps[timestampCounter]
        endStamp = timestamps[timestampCounter + 1]

        trialPath = os.getcwd() + '\\' + \
                    'Trial{0:02d}.txt'.format(trialNumber + 1)
        currTrialData = numpy.genfromtxt(trialPath, skip_header=1)

        ## obtain trialParameters from fullData timestamps
        currParameters = numpy.zeros(shape=(4))
        currParameters = get_trial_parameters(currTrialData)

        ## Assign color background based on gain
        if abs(currParameters[2]) == 0:
            currentColor = 'yellow'
        elif abs(currParameters[2]) == 2.0:
            currentColor = 'green'
        elif abs(currParameters[2]) == 4.0:
            currentColor = 'violet'
        else:
            #should never happen
            currentColor = 'cyan'

        ## Assign opacity based on direction
        if currParameters[1] == 1:
            currentAlpha = 1.00
        elif currParameters[1] == -1:
            currentAlpha = 0.30

        plt.axvspan(xmin = startStamp, xmax = endStamp, color = currentColor, 
                    alpha = currentAlpha)

        # Advance to next timestamp pair
        timestampCounter += 2

    ## Save figure to PDF
    pp.savefig()
    plt.close('all')            # close figure
    print "\nPDF closed.\nFull data graphing 1 completed!\n"

    print "Graphing reordered trials..."
    fig = plt.figure(2)
    figureTitle = '{} Raw Data (Reordered)'.format(blockNumberString)
    fig.suptitle(figureTitle, fontsize = 12, fontweight = 'bold')

    ax = plt.subplot(111)
    ax.set_xlabel('Time')
    ax.set_ylabel('Grams')

    subplotTitle = 'Applied Force (GS0 Force Sensor)'
    plt.title(subplotTitle, fontsize = 10)
    plt.axis(ymin = -30, ymax = 30)

    # Begin graphing raw data, ordered so mag/unmag trials are adjacent
    lastEndStamp = 0
    for trialNumber in range(0, 6):
        trialName = orderedBlockDataList[trialNumber]
        trialPath = os.getcwd() + '\\' + trialName

        currTrialData = numpy.genfromtxt(trialPath, skip_header=1)
        GS0_ForceGramsTrial = currTrialData[:,11]

        startStamp = lastEndStamp + 1
        endStamp = startStamp + len(GS0_ForceGramsTrial)

        x = numpy.arange(startStamp, endStamp)
        plt.plot(x, GS0_ForceGramsTrial, 'b', label='GS0-100 Force')

        currParameters = trialParameters[trialNumber,:]

        ## Assign color background based on distal force
        if abs(currParameters[2]) == 0:
            currentColor = 'yellow'
        elif abs(currParameters[2]) == 2.0:
            currentColor = 'green'
        elif abs(currParameters[2]) == 4.0:
            currentColor = 'violet'
        else:
            currentColor = 'cyan'

        ## Assign opacity based on magnified/unmagnified
        if currParameters[1] == 1:
            currentAlpha = 1.00
        elif currParameters[1] == -1:
            currentAlpha = 0.30

        plt.axvspan(xmin = startStamp, xmax = endStamp, color = currentColor, 
                    alpha = currentAlpha)

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
    blockDataList = numpy.empty(6, dtype = 'str')
    blockDataList = reorder_trials(blockNumberString)
    print 'blockDataList = {}\n'.format(blockDataList)

    # Obtain number of files in subjectData directory
    blockDataListSize = len(blockDataList)

    ## open PDF for plots and figures
    pdfFileName = os.getcwd() +'\\'+ blockNumberString + '_FourierPlots.pdf'
    pp = PdfPages(pdfFileName)            #Pointer to pdf file

    ## create empty arrays and flags
    print "Reading trial data to graph FFT...\n"

    ## Graph and obtain trial parameters for all trials in block
    doneFlag = 0    
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
        currentData = numpy.genfromtxt(currentFilePath, skip_header=1)

        # 16 trials, 3 columns each
        # [0] trial_ID [1] direction [2] mag [3] perceived_target
        print "(Trial template #{})\n".format(trialParameters[trialCounter,0])

        ## Load data from current trial, calculate FFT and freq
        # [0] block#; [1] trial#; [2] trial timestamp; [3] GS0Voltage;
        # [4] HHFM_out Voltage; [5] HHFM_in Voltage; [6] bit0; [7] bit1;
        # [8] bit2; [9] currentlyInTrial voltage;
        # [10] GS0Force; [11] filteredForce;
        GS0_ForceGrams = currentData[:,10]

        endStamp = len(GS0_ForceGrams) - 2500
        startStamp = endStamp - 6000
        rawForceNoFD = GS0_ForceGrams[startStamp + 1000 : endStamp]

        fftRawNoFD = numpy.abs(numpy.fft.fft(rawForceNoFD, axis = 0))

        # timestep in units of seconds to obtain freq in cycles per second
        timestep = 0.001
        freqNoFD = numpy.fft.fftfreq(rawForceNoFD.size, d = timestep)

        # Set up figures/axes/labels
        fig = plt.figure(figCounter)
        print "Figure #{}...\n".format(figCounter)
        
        if trialParameters[trialCounter,1] == -1:
            direction = 'Push'
        elif trialParameters[trialCounter,1] == 1:
            direction = 'Pull'
        figureTitle = 'Trial#{} - {} (Direction: {}({}); Mag = {}x)'.format(
                    figCounter, currentFileName, direction, 
                    int(trialParameters[trialCounter, 1]), 
                    int(trialParameters[trialCounter, 2]))
        fig.suptitle(figureTitle, fontsize = 12, fontweight = 'bold')

        ax1L = plt.subplot(211)
        subplotTitle = 'Trial #{}; FFT Applied Force '\
                        '(GS0 Force Sensor)'.format(figCounter)
        ax1L.set_title(subplotTitle, fontsize = 8)
        ax1L.set_xlim(left = 1, right = 10)
        ax1L.set_ylim(bottom = 0, top = 5000)
        plt.tick_params(axis = 'both', which = 'major', labelsize = 8)

        y_label = r'$\mathcal{F}(Force)$'
        ax1L.set_ylabel(y_label, fontsize = 8)
        ax1L.set_xlabel('Frequency (Hz)', fontsize = 8)

        ## Plot FFT
        fftNoFD, = ax1L.plot(freqNoFD, fftRawNoFD, 'b', 
                             label = 'FFT Raw Force (no feedback)', lw = 1.00)
        
        plt.legend(loc = 'best', fontsize = 8)

        # Save FFT as text file
        temp = numpy.column_stack((freqNoFD, fftRawNoFD))

        fileName = os.getcwd() +'\\'+ blockNumberString + \
                   '_{}_FourierCoeff.txt'.format(fileCounter + 1)

        numpy.savetxt(fileName, temp, fmt = '%.4e', 
                      header = 'Freq NoFD\t fftRawNoFD', delimiter = '\t')

        # take out all negative frequencies
        # @TODO: maybe not needed
        conditions = [freqNoFD >= 0]
        choices = [freqNoFD]
        posFreq = numpy.trim_zeros(numpy.select(conditions,choices), 'b')
        
        # multiply by conjugate because int{ |F(w)|^2 } = int{ |f(t)|^2 }
        # by perseval's formula. Returns measure of "energy" in signal
        # from index 5 - 21 because those match with posFreq = 1 -> 4 Hz
        index = 0
        pwrNoFD = numpy.zeros(shape=(16))
        sumNoFD = numpy.zeros(shape=(16))
        for i in range(5, 21):
            pwrNoFD[index] = (fftRawNoFD[i] * numpy.conjugate(fftRawNoFD[i]))
            index += 1

        # calculate integral for pwr from 1 - 4 Hz
        for index in range(len(pwrNoFD)):
            if index == 0:
                sumNoFD[index] = pwrNoFD[index]
            else:
                sumNoFD[index] += pwrNoFD[index] + sumNoFD[index - 1]

        ## Second plot, for integral of FFT
        ax2L = plt.subplot(212)

        y_label = r'$\int\/\/|\mathcal{F}(Force)|^2 $'
        ax2L.set_ylabel(y_label, fontsize = 8)
        plt.xlabel('Frequency (Hz)', fontsize = 8)

        ax2L.set_xlim(left = 1, right = 4)
        ax2L.set_ylim(bottom = 0, top = 2e8)
        plt.tick_params(axis = 'both', which = 'major', labelsize = 8)

        intFFTNoFD, = ax2L.plot(posFreq[5:21], sumNoFD, 'b', lw = 1.00,
                                label = 'Cumulative Power Spectrum (no FD)')

        plt.legend(loc = 'best', fontsize = 8)

        ## Save figure to PDF
        pp.savefig(figCounter)
        plt.close('all')
        figCounter += 1

        ## Done graphing, advance to next file/ trial
        fileCounter += 1
        trialCounter += 1

        if fileCounter >= blockDataListSize:
            plt.close('all')
            pp.close()
            print "\nPDF closed.\n"
            print "{} split data graphing completed!\n".format(blockNumberString)
            doneFlag = 1

    return

def reorder_trials(blockNumberString):

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

def getBinaryFromVoltage(voltage):
    '''Infer binary state from float voltage data.

    Returns 1 if voltage > 1.00'''
    if voltage > 1.00:
        # one pin is max @ ~1.50V
        result = 1.
    else:
        result = 0.

    return result

def saturation_check(actuatorVoltage):
    ''' Check and save timestamps where actuator has saturated
    
    '''
    
    windowSize = 20
    window_start = numpy.arange(0,len(actuatorVoltage), windowSize)
        
    timestamps = numpy.empty(0)
    for i in range(len(window_start)):
        left = window_start[i]
        right = left + 20        
        window = actuatorVoltage[left : right]
        if window.mean() < 3.0 or window.mean() > 6.0:
            timestamps = numpy.append(timestamps, left)
    
    # if there is no saturation, array will still be empty
    if timestamps.size == 0:
        timestamps = numpy.append(timestamps, -1)
        
    return timestamps

def filter_forces(fullDataArray):
    ''' Convolve raw force trace with 201 and 501 element normalized Gaussians.
    
    '''

    # 21 element gaussian, sigma = 4
    #gaussian = numpy.array([0.000272337, 0.00089296, 0.002583865, 0.00659813,
    #                        0.014869116, 0.029570767, 0.051898313, 0.080381679,
    #                        0.109868729, 0.132526984, 0.14107424, 0.132526984,
    #                        0.109868729, 0.080381679, 0.051898313, 0.029570767,
    #                        0.014869116, 0.00659813, 0.002583865, 0.00089296,
    #                        0.000272337])

    # load 201 element gaussian file, sigma = 100:
    fileName = rawDataDirectory + '\\' + '201element_normGaussian.txt'
    gaussian201 = numpy.genfromtxt(fileName, skip_header = 1)

    # load 501 element gaussian file, sigma = 200:
    fileName = rawDataDirectory + '\\' + '501element_normGaussian.txt'
    gaussian501 = numpy.genfromtxt(fileName, skip_header = 1)

    GS0_ForceGrams = fullDataArray[:,10]
    size = GS0_ForceGrams.shape
    filteredForces201 = numpy.zeros(shape = size)
    filteredForces501 = numpy.zeros(shape = size)
    
    # Loop for 21 element gaussian
    #for timestamp in range(11, size[0] - 10):
    #    rawChunk = GS0_ForceGrams[timestamp - 11 : timestamp + 10]
    #
    #    pointAverage = numpy.convolve(rawChunk, gaussian, mode = 'valid')
    #
    #    filteredForces[timestamp] = pointAverage

    # Loop for 201 element gaussian convolution
    for timestamp in range(101, size[0] - 100):
       rawChunk = GS0_ForceGrams[timestamp - 101 : timestamp + 100]
    
       pointAverage = numpy.convolve(rawChunk, gaussian201[:,1], mode = 'valid')
    
       filteredForces201[timestamp] = pointAverage

    # Loop for 501 element gaussian convolution
    for timestamp in range(251, size[0] - 250):
        rawChunk = GS0_ForceGrams[timestamp - 251 : timestamp + 250]

        pointAverage = numpy.convolve(rawChunk, gaussian501[:,1], mode = 'valid')

        filteredForces501[timestamp] = pointAverage

    return filteredForces201, filteredForces501

###############################################################################
##
## MAIN : Process Data
##
###############################################################################

# Call time to set as start time
start = time.time()

print "\n\n*** MinContact_RawDataProcessing.py - START ***\n\n"

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
    newBlockFolderName = '{}_Block{}'.format(subjectNumber, i+1)
    newBlockDir = os.getcwd() + '\\' + newBlockFolderName
    os.mkdir(newBlockDir)
    os.chdir(newBlockDir)

    ## Currently in HHFM_Data/subjectNumber/Blockx/
    # extract raw data
    subjectDataFullPath = subjectDirectory + '\\' + subjectDataList[i]
    initFullData = numpy.genfromtxt(subjectDataFullPath, skip_header = 1,
                                    skip_footer = 1)
    initFullDataSize = initFullData.shape

    ## Create fullData matrix, with extra column for GS0 Force calculation
    # Creates individual text files for each trial
    # Extract start/stop timestamps
    fullData, timestamps = split_trials(initFullData, i+1)

    timestampFilePath = os.getcwd() + '\\' + newBlockFolderName + \
                        '_timestamps.txt'
    numpy.savetxt(timestampFilePath, timestamps, fmt = '%d')
    print "Block{} timestamps saved.\n".format(i+1)

    ## Process split trials and get list of trial parameters for the current block
    # Creates data graphs for each trial in PDF, save parameters for block
    print "Block{} trial processing...\n".format(i+1)
    savePath = os.getcwd()
    trialParameters = process_split_trials(savePath, timestamps, newBlockFolderName)

    # Save trialParameters array to text file
    trialParametersFilePath = os.getcwd() + '\\' + newBlockFolderName + \
                              '_trialParameters.txt'
    numpy.savetxt(trialParametersFilePath, trialParameters, fmt = '%.2f')

    ## Summary graph of full raw data with annotations
    graph_full_data(fullData, timestamps, trialParameters, newBlockFolderName)

    ## Graph and save fourier transforms, power spectra
    graph_fourier(trialParameters, newBlockFolderName)

    print "Block{} done!\n".format(i + 1)

    ## Go back to HHFM_Data/subjectNumber
    os.chdir(subjectDirectory)


print "{} processing Completed!".format(subjectNumber)
print "Time elapsed: {0:03f} seconds".format(time.time() - start)
print "\n\n*** MinContact_RawDataProcessing.py - END ***\n\n"
