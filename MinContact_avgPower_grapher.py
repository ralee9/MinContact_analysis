# -*- coding: utf-8 -*-
"""
Created on Fri Aug 05 23:38:40 2016

@author: Randy Lee
"""
##############################################################################
##  Randy Lee
##  VIA Lab 
##  Department of Bioengineering
##  University of Pittsburgh
##
##  Plot contact time CDF and histograms by trial z scores
##
##  Last updated: 05 August 2016
###############################################################################

import os
import numpy
import time
import math
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# edit these constants to access figure
#rootFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts'
#dataFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts\data'
#outputFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts\output'
rootFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts'
dataFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts\data'
outputFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts\output'

data = '\\MinContact_avgPower.txt'

pdfSaveName = '\\MinContact_powerPercentage.pdf'

# font sizes for axis labels, legends, and other text
label_fontSize = 16
legend_fontSize = 12
text_fontSize = 14
tick_fontSize = 12

# set tick labels to Times New Roman font to match LaTeX axis labels
ticks_font = matplotlib.font_manager.FontProperties(family = 'Times New Roman',
                                                    style = 'normal',
                                                    size = 14)


def plot_data(data, title, y_label, savePath):
    
    N = 6
    
    # extract data
    mean = data[:,0]/5.
    SD = data[:,1]/5.
    SEM = SD/math.sqrt(N)
    
    gain = numpy.array([1, 2, 3], dtype = 'int64')
    ticks = [r'$\mathrm{1-4\ Hz}$', r'$\mathrm{4-7\ Hz}$', 
             r'$\mathrm{7-10\ Hz}$']
    #width = 0.05
    
    pp = PdfPages(savePath)

    fig = plt.figure()
    fig.suptitle(title, fontsize = 8, fontweight = 'bold', y = 0.99)   
    
    ax = plt.subplot(111)
    plt.errorbar(gain, mean, yerr = SEM, color = 'k', marker = 'o',
                 linestyle = '-', mec = 'k', ms = 4.0,
                 mew = 1.1, capsize=0)
     
    # set up rest of figure & labels  
    plt.grid(which = 'major', axis = 'y')
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    plt.xticks(gain, ticks)
    ax.tick_params(axis = 'both', labelsize = tick_fontSize)
    plt.xlim(xmin = 0, xmax = 4)
    #plt.axvline(x = 0, linewidth = 1.0, color = 'k')
    plt.xlabel(r'$\mathrm{Frequency\ Band}$', fontsize = label_fontSize)
    plt.ylabel(y_label, fontsize = label_fontSize)
    plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
    
    # update tick labels to correct font
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    

    # remove error bars in legend
    #handles, labels = ax.get_legend_handles_labels()
    #handles = [h[0] for h in handles]
    #ax.legend(handles, labels, loc = 'upper right', fontsize = 8, numpoints = 1)    
    plt.tight_layout()
    
    pp.savefig(fig)
    plt.close('all')
    pp.close()    
    
    return




###############################################################################
##  MAIN
###############################################################################

start = time.time()

print "\n\n*** MinContact_avgPower_grapher.py - START *** \n\n"

os.chdir(dataFolder)

rawData = numpy.genfromtxt(dataFolder+data, skip_header = 1)

savePath = outputFolder + pdfSaveName


#plot data
plot_data(rawData,' ', r'$\mathrm{Power}\ (grams^2 \; Hz)$', savePath)


print "\nGraphing done!!"
print "Time elapsed: {0:03f} seconds".format(time.time() - start)
print "\n\n*** MinContact_avgPower_grapher.py - END ***\n\n"
