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
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# edit these constants to access figure
rootFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts'
dataFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts\data'
#rootFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts'
#dataFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts\data'
mag0_contact = '\\MinContact_mag0_contact_zscore.txt'
mag2_contact = '\\MinContact_mag2_contact_zscore.txt'
mag4_contact = '\\MinContact_mag4_contact_zscore.txt'

pdfSaveName_combined = '\\MinContact_contact_combined_zscore.pdf'
pdfSaveName_CDF = '\\MinContact_contact_CDF_zscore.pdf'
pdfSaveName_hist = '\\MinContact_contact_hist_zscore.pdf'


# edit these constants to customize figure. Data should be organized in columns
title = 'Time in Contact'
x_label = 'z-score'
CDF_ylabel = 'Cumulative Probabilty'
Hist_ylabel = 'Frequency'


def plot_data(mag0, mag2, mag4, xlabel, ylabel1, ylabel2, savePath):
    
    pp = PdfPages(savePath)

    fig = plt.figure()
    fig.suptitle(title, fontsize = 8, fontweight = 'bold', y = 0.99)
    
    # plot histograms
    ax1 = plt.subplot(121)
    
#    #mag 0
#    n, mag0_bins, patches = plt.hist(mag0, bins = 50, color = 'b', 
#                                     histtype = 'bar', label = 'Gain = 0')
#    #mag 2
#    plt.hist(mag2, bins = mag0_bins, color = 'r', label = 'Gain = 2', 
#             alpha = 0.50, histtype = 'bar')
#    #mag 4
#    plt.hist(mag4, bins = mag0_bins, color = 'orange', label = 'Gain = 4', 
#             alpha = 0.75, histtype = 'bar')
    plt.hist([mag0,mag2,mag4], bins = 30, histtype = 'bar', range = (-5.,5.),
             label = ['Gain = 0', 'Gain = 2', 'Gain = 4'], fill = True)
    
    # set up rest of figure & labels
    plt.grid(True)
    plt.xlabel(xlabel, fontsize = 8)
    plt.ylabel(ylabel1, fontsize = 8)
    plt.xlim(xmin = -5, xmax = 5)
    plt.ylim(ymin = 0, ymax = 40)
    plt.legend(loc = 'upper right', fontsize = 6)
    
    ## plot CDF
    ax2 = plt.subplot(122)
    
    N = mag0.size
    sort_mag0 = numpy.sort(mag0)
    sort_mag2 = numpy.sort(mag2)
    sort_mag4 = numpy.sort(mag4)
    y2 = numpy.array(range(N))/float(N)
    
    plt.plot(sort_mag0, y2, color = 'b', label = 'Gain = 0')
             
    plt.plot(sort_mag2, y2, color = 'r', label = 'Gain = 2')
    
    plt.plot(sort_mag4, y2, color = 'g', label = 'Gain = 4')

    # set up rest of figure & labels  
    plt.grid(True)  
    plt.xlim(xmin = -5, xmax = 5)
    plt.ylim(ymin = 0, ymax = 1)
    plt.legend(loc = 'lower right', fontsize = 6)
    plt.xlabel(xlabel, fontsize = 8)
    plt.ylabel(ylabel2, fontsize = 8)

    plt.tight_layout()
    
    pp.savefig(fig)
    plt.close('all')
    pp.close()    
    
    return

def plot_CDF(mag0, mag2, mag4, x_label, y_label, savePath):
    pp = PdfPages(savePath)   
    
    #plot CDF together, push
    fig = plt.figure()
    title = 'Time in Contact, Cumulative Distribution Function by Z-score (push)'
    fig.suptitle(title, fontsize = 8, fontweight = 'bold')
    
    sort_mag0_push = numpy.sort(mag0[0:65])
    sort_mag2_push = numpy.sort(mag2[0:65])
    sort_mag4_push = numpy.sort(mag4[0:65])
    N = sort_mag0_push.size
    y = numpy.array(range(N))/float(N)
    
    plt.plot(sort_mag0_push, y, color = 'b', label = 'Gain = 0')
             
    plt.plot(sort_mag2_push, y, color = 'r', label = 'Gain = 2')
             
    plt.plot(sort_mag4_push, y, color = 'g', label = 'Gain = 4')

    # set up rest of figure & labels  
    plt.grid(True)  
    plt.xlim(xmin = -5, xmax = 5)
    plt.ylim(ymin = 0, ymax = 1)
    plt.legend(loc = 'lower right', fontsize = 8)
    plt.xlabel(x_label, fontsize = 8)
    plt.ylabel(y_label, fontsize = 8)
    plt.tight_layout()
    pp.savefig(fig)
    
    ## pull CDF
    fig2 = plt.figure()
    title = 'Time in Contact, Cumulative Distribution Function by Z-score (pull)'
    fig2.suptitle(title, fontsize = 8, fontweight = 'bold')
    
    sort_mag0_pull = numpy.sort(mag0[65:])
    sort_mag2_pull = numpy.sort(mag2[65:])
    sort_mag4_pull = numpy.sort(mag4[65:])
    N = sort_mag0_pull.size
    y = numpy.array(range(N))/float(N)
    
    plt.plot(sort_mag0_pull, y, color = 'b', label = 'Gain = 0')
             
    plt.plot(sort_mag2_pull, y, color = 'r', label = 'Gain = 2')
             
    plt.plot(sort_mag4_pull, y, color = 'g', label = 'Gain = 4')

    # set up rest of figure & labels  
    plt.grid(True)  
    plt.xlim(xmin = -5, xmax = 5)
    plt.ylim(ymin = 0, ymax = 1)
    plt.legend(loc = 'lower right', fontsize = 8)
    plt.xlabel(x_label, fontsize = 8)
    plt.ylabel(y_label, fontsize = 8)
    plt.tight_layout()
    
    pp.savefig(fig2)
    
    # combined CDF, by mag
    fig3 = plt.figure()
    title = 'Time in Contact, Cumulative Distribution Function by Z-score (combined)'
    fig3.suptitle(title, fontsize = 8, fontweight = 'bold')
    
    sort_mag0_all = numpy.sort(mag0)
    sort_mag2_all = numpy.sort(mag2)
    sort_mag4_all = numpy.sort(mag4)
    N = sort_mag0_all.size
    y = numpy.array(range(N))/float(N)
    
    plt.plot(sort_mag0_all, y, color = 'b', label = 'Gain = 0')
                 
    plt.plot(sort_mag2_all, y, color = 'r', label = 'Gain = 2')
             
    plt.plot(sort_mag4_all, y, color = 'g', label = 'Gain = 4')
    
    plt.legend(loc = 'lower right', fontsize = 8)
    plt.grid(True)
    plt.xlim(xmin = -5, xmax = 5)
    plt.ylim(ymin = 0, ymax = 1)
    plt.xlabel(x_label, fontsize = 8)
    plt.ylabel(y_label, fontsize = 8)
    plt.tight_layout()
    
    pp.savefig(fig3)
    plt.close('all')
    pp.close() 
    
    # calculate and print K-S test stats
    print "**PUSH ONLY**"
    print "K-S mag0 vs mag2:"
    D_0vs2, p_0vs2 = stats.ks_2samp(sort_mag0_push, sort_mag2_push)
    print "D_0vs2 = {};\tp_0vs2 = {}".format(D_0vs2, p_0vs2)
    
    print "K-S mag0 vs mag4:"
    D_0vs4, p_0vs4 = stats.ks_2samp(sort_mag0_push, sort_mag4_push)
    print "D_0vs4 = {};\tp_0vs4 = {}".format(D_0vs4, p_0vs4)
    
    print "K-S mag2 vs mag4:"
    D_2vs4, p_2vs4 = stats.ks_2samp(sort_mag2_push, sort_mag4_push)
    print "D_2vs4 = {};\tp_2vs4 = {}".format(D_2vs4, p_2vs4)    
    
    print "\n\n**PULL ONLY**"
    print "K-S mag0 vs mag2:"
    D_0vs2, p_0vs2 = stats.ks_2samp(sort_mag0_pull, sort_mag2_pull)
    print "D_0vs2 = {};\tp_0vs2 = {}".format(D_0vs2, p_0vs2)
    
    print "K-S mag0 vs mag4:"
    D_0vs4, p_0vs4 = stats.ks_2samp(sort_mag0_pull, sort_mag4_pull)
    print "D_0vs4 = {};\tp_0vs4 = {}".format(D_0vs4, p_0vs4)
    
    print "K-S mag2 vs mag4:"
    D_2vs4, p_2vs4 = stats.ks_2samp(sort_mag2_pull, sort_mag4_pull)
    print "D_2vs4 = {};\tp_2vs4 = {}".format(D_2vs4, p_2vs4)
    
    print "\n\n**PUSH vs PULL**"
    print "K-S push vs pull @ mag0:"
    D_pushpull_0, p_pushpull_0 = stats.ks_2samp(sort_mag0_push, sort_mag0_pull)
    print "D_pushpull_0 = {};\tp_pushpull_0 = {}".format(D_pushpull_0, p_pushpull_0)
    
    print "K-S push vs pull @ mag2:"
    D_pushpull_2, p_pushpull_2 = stats.ks_2samp(sort_mag2_push, sort_mag2_pull)
    print "D_pushpull_2 = {};\tp_pushpull_2 = {}".format(D_pushpull_2, p_pushpull_2)
    
    print "K-S push vs pull @ mag4:"
    D_pushpull_4, p_pushpull_4 = stats.ks_2samp(sort_mag4_push, sort_mag4_pull)
    print "D_pushpull_4 = {};\tp_pushpull_4 = {}".format(D_pushpull_4, p_pushpull_4)
    
    print "\n\n**COMBINED, by MAG**"
    print "K-S mag0 vs mag4:"
    D_0vs2_all, p_0vs2_all = stats.ks_2samp(sort_mag0_all, sort_mag2_all)
    print "D_0vs2 = {};\tp_0vs2 = {}".format(D_0vs2_all, p_0vs2_all)
    
    print "K-S mag0 vs mag4:"
    D_0vs4_all, p_0vs4_all = stats.ks_2samp(sort_mag0_all, sort_mag4_all)
    print "D_0vs4 = {};\tp_0vs4 = {}".format(D_0vs4_all, p_0vs4_all)
    
    print "K-S mag2 vs mag4:"
    D_2vs4_all, p_2vs4_all = stats.ks_2samp(sort_mag2_all, sort_mag4_all)
    print "D_2vs4 = {};\tp_2vs4 = {}".format(D_2vs4_all, p_2vs4_all)
        
    return

def plot_histograms(mag0, mag2, mag4, x_label, y_label, savePath):
    
    pp = PdfPages(savePath)

    # histogram subplots
    fig1 = plt.figure()
    fig1.suptitle(title, fontsize = 8, fontweight = 'bold', y = 0.99)
    
    plt.subplot(311)
    n, mag0_bins, patches = plt.hist(mag0, bins = 30, color = 'b', 
                                     range = (-5.,5.), label = 'Gain = 0',
                                     histtype = 'stepfilled')
    plt.xlim(xmin = -5, xmax = 5)
    plt.ylim(ymin = 0, ymax = 40)
    plt.ylabel(y_label, fontsize = 8)
    plt.grid(True)
    plt.legend(loc = 'upper right', fontsize = 8)
    
    plt.subplot(312)
    plt.hist(mag2, bins = mag0_bins, color = 'r', label = 'Gain = 2', 
             histtype = 'stepfilled')
    plt.xlim(xmin = -5, xmax = 5)
    plt.ylim(ymin = 0, ymax = 40)
    plt.ylabel(y_label, fontsize = 8)
    plt.grid(True)
    plt.legend(loc = 'upper right', fontsize = 8)
    
    plt.subplot(313)
    plt.hist(mag4, bins = mag0_bins, color = 'orange', label = 'Gain = 4',
             histtype = 'stepfilled')
    plt.xlim(xmin = -5, xmax = 5)    
    plt.ylim(ymin = 0, ymax = 40)
    plt.xlabel(x_label, fontsize = 8)
    plt.ylabel(y_label, fontsize = 8)
    plt.grid(True)
    plt.legend(loc = 'upper right', fontsize = 8)
    
    plt.tight_layout()
    pp.savefig(fig1)
    
    # normalized histograms and PDFs, model fits
    fig2 = plt.figure()
    fig2.suptitle(title, fontsize = 8, fontweight = 'bold', y = 0.99)
    
    plt.subplot(311)
    n2, mag0_bins2, patches2 = plt.hist(mag0, bins = 30, color = 'b',
                                        range = (-5., 5), label = 'Gain = 0',
                                        normed = True, histtype = 'stepfilled')
    mean_0 = numpy.mean(mag0)
    std_0 = numpy.std(mag0)
    
    x_plot = numpy.linspace(-5, 5, 10000)
    plt.plot(x_plot, stats.norm.pdf(x_plot, mean_0, std_0), 'k--', 
             label = 'norm fit')
    plt.grid(True)
    plt.ylim(ymin = 0, ymax = 1)
    plt.ylabel('Probability', fontsize = 8)
    plt.legend(loc = 'best', fontsize = 8)
    
    plt.subplot(312)
    mean_2 = numpy.mean(mag2)
    std_2 = numpy.std(mag2)
    
    plt.hist(mag2, bins = mag0_bins2, color = 'r', normed = True, 
             label = 'Gain = 2', histtype = 'stepfilled')    
    plt.plot(x_plot, stats.norm.pdf(x_plot, mean_2, std_2), 'k--', 
             label = 'norm fit')
    plt.grid(True)    
    plt.ylim(ymin = 0, ymax = 1)
    plt.ylabel('Probability', fontsize = 8)
    plt.legend(loc = 'best', fontsize = 8)
    
    plt.subplot(313)
    mean_4 = numpy.mean(mag4)
    std_4 = numpy.std(mag4)
    plt.hist(mag4, bins = mag0_bins2, color = 'g', normed = True, 
             label = 'Gain = 4', histtype = 'stepfilled')
    plt.plot(x_plot, stats.norm.pdf(x_plot, mean_4, std_4), 'k--', 
             label = 'norm fit')
    #plt.plot(x_plot, stats.gamma.pdf())
    plt.grid(True)
    plt.ylim(ymin = 0, ymax = 1)
    plt.xlabel('z-score', fontsize = 8)
    plt.ylabel('Probability', fontsize = 8)
    plt.legend(loc = 'best', fontsize = 8)
    
    plt.tight_layout()
    pp.savefig(fig2)
       
    
    # histograms combined
    fig3 = plt.figure()
    fig3.suptitle(title, fontsize = 8, fontweight = 'bold', y = 0.99)
    
#    n, mag0_bins, patches = plt.hist(mag0, bins = 50, color = 'b', 
#                                     label = 'Gain = 0')
#    #mag 2
#    plt.hist(mag2, bins = mag0_bins, color = 'r', label = 'Gain = 2',
#             alpha = 0.50)
#    #mag 4
#    plt.hist(mag4, bins = mag0_bins, color = 'orange', label = 'Gain = 4', 
#             alpha = 0.75)
    
    plt.hist([mag0,mag2,mag4], bins = 30, histtype = 'bar', range = (-5.,5.),
             label = ['Gain = 0', 'Gain = 2', 'Gain = 4'], fill = True,
             color = ['b', 'r', 'g'])
    
    # set up rest of figure & labels
    plt.grid(True)
    plt.xlabel(x_label, fontsize = 8)
    plt.ylabel(y_label, fontsize = 8)
    plt.xlim(xmin = -5, xmax = 5)
    plt.ylim(ymin = 0, ymax = 40)
    plt.legend(loc = 'upper right', fontsize = 8)
    
    plt.tight_layout()
    
    pp.savefig(fig3)
    plt.close('all')
    pp.close()
    
    return
    
    

###############################################################################
##  MAIN
###############################################################################

start = time.time()

print "\n\n*** MinContact_contactTimeGrapher.py - START *** \n\n"

os.chdir(dataFolder)

mag0_rawData = numpy.genfromtxt(dataFolder+mag0_contact, skip_header = 1)
mag2_rawData = numpy.genfromtxt(dataFolder+mag2_contact, skip_header = 1)
mag4_rawData = numpy.genfromtxt(dataFolder+mag4_contact, skip_header = 1)

mag0_zscore = mag0_rawData[:,8]
mag2_zscore = mag2_rawData[:,8]
mag4_zscore = mag4_rawData[:,8]

savePath_combined = rootFolder+pdfSaveName_combined
savePath_CDF = rootFolder+pdfSaveName_CDF
savePath_hist = rootFolder+pdfSaveName_hist

#plot CDF and hist together
plot_data(mag0_zscore, mag2_zscore, mag4_zscore, x_label, Hist_ylabel,
          CDF_ylabel, savePath_combined)

#plot CDFs side by side
plot_CDF(mag0_zscore, mag2_zscore, mag4_zscore, x_label, CDF_ylabel, 
         savePath_CDF)

#plot histograms side by side
plot_histograms(mag0_zscore, mag2_zscore, mag4_zscore, x_label, Hist_ylabel,
                savePath_hist)

print "\nAnalysis done!!"
print "Time elapsed: {0:03f} seconds".format(time.time() - start)
print "\n\n*** MinContact_contactTimeGrapher.py - END ***\n\n"
