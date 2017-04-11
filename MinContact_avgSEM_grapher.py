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
import matplotlib
import matplotlib.patches as patches
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import AutoMinorLocator
import matplotlib.pyplot as plt

# edit these constants to access figure
#rootFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts'
#dataFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts\data'
#outputFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts\output'
rootFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts'
dataFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts\data'
outputFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts\output'

data_mean = '\\MinContact_meanForce_avgSEM.txt'
data_SD = '\\MinContact_SDForce_avgSEM.txt'
data_pwr14 = '\\MinContact_pwr14_avgSEM.txt'
data_pwr47 = '\\MinContact_pwr47_avgSEM.txt'
data_pwr710 = '\\MinContact_pwr710_avgSEM.txt'
data_contact = '\\MinContact_contact_avgSEM.txt'

data_LFEmean_noFD = '\\LFE_meanForce_noFD.txt'
data_LFEmean_visFD = '\\LFE_meanForce_visFD.txt'
data_LFEstd_noFD = '\\LFE_stddev_noFD.txt'
data_LFEstd_visFD = '\LFE_stddev_visFD.txt'

data_LFEz_mag0 = '\\LFE_Zscore_mag0.txt'
data_LFEz_mag1 = '\\LFE_Zscore_mag1.txt'

pdfSaveName_mean = '\\MinContact_avgSEM_mean.pdf'
pdfSaveName_SD = '\\MinContact_avgSEM_SD.pdf'
pdfSaveName_pwr14 = '\\MinContact_avgSEM_pwr14.pdf'
pdfSaveName_pwr47 = '\\MinContact_avgSEM_pwr47.pdf'
pdfSaveName_pwr710 = '\\MinContact_avgSEM_pwr710.pdf'
pdfSaveName_power = '\\MinContact_avgSEM_powerCombined.pdf'
pdfSaveName_contact = '\\MinContact_avgSEM_contact.pdf'

pdfSaveName_LFEmean_noFD = '\\LFE_mean_noFD.pdf'
pdfSaveName_LFEmean_visFD = '\\LFE_mean_visFD.pdf'
pdfSaveName_LFEstd_noFD = '\\LFE_std_noFD.pdf'
pdfSaveName_LFEstd_visFD = '\\LFE_std_visFD.pdf'

pdfSaveName_LFEz_CDF = '\\LFE_windowTime_Zscore_CDF.pdf'

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
    
    N = 13
    
    # extract data
    if (title == 'Power (1 to 4 Hz)' or 
    title == 'Power (4 to 7 Hz)' or
    title == 'Power (7 to 10 Hz)'):
        # account for 
        mean = data[:,0]/5.
        SD = data[:,1]/5.
    else:
        mean = data[:,0]
        SD = data[:,1]
    
    push_mean = mean[0:3]
    push_SD = SD[0:3]
    push_SEM = push_SD/math.sqrt(N)
    
    pull_mean = mean[3:]
    pull_SD = SD[3:]
    pull_SEM = pull_SD/math.sqrt(N)
    
    gain = numpy.array([1, 2, 3, 4, 5], dtype = 'int64')
    data_gain = numpy.array([1,3,5],dtype='int64')
    ticks = ['0','1', '2','3', '4']
    width = 0.05
    
    pp = PdfPages(savePath)

    fig = plt.figure()
    #fig.suptitle(title, fontsize = 8, fontweight = 'bold', y = 0.99)   
    
    ax = plt.subplot(111)
    plt.errorbar(data_gain, push_mean, yerr = push_SEM, color = 'k', marker = 'o',
                 linestyle = '--', mec = 'k', label = r'$\mathrm{Push}$', 
                 ms = 6.5, mew = 1.1, capsize=0, fillstyle='none', dashes=(4,3))
    
    plt.errorbar(data_gain+width, pull_mean, yerr = pull_SEM, color = 'k',
                 marker = 'o', linestyle = '-', label = r'$\mathrm{Pull}$',
                 ms = 5.0, mec = 'k', mew = 1.1, capsize=0)
    
    # set up rest of figure & labels  
    plt.grid(which = 'major', axis = 'y')
    plt.ylim(ymin = 0)
    if title == 'Contact Time':
        plt.ylim(ymax = 5)
    if title == 'Power (1 to 4 Hz)':
        plt.ylim(ymax = 1.2e7)
    if (title == 'Power (4 to 7 Hz)' or title == 'Power (7 to 10 Hz)'):
        plt.ylim(ymax = 2.5e6)
        plt.ticklabel_format(style = 'sci', axis = 'y', scilimits = (0,0))
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    #plt.axvline(x = 0, linewidth = 1.0, color = 'k')
    
    # relabel x ticks & add minor labels; add x & y labels
    plt.xlim(xmin = 0, xmax = 6)
    plt.xticks(gain, ticks)
    minor_locator = AutoMinorLocator(2)
    ax.xaxis.set_minor_locator(minor_locator)
    ax.tick_params(axis = 'x', which='both', labelsize = tick_fontSize)
    ax.tick_params(axis = 'y', labelsize = tick_fontSize)
    ax.set_xticks([1,3], minor = True)
    plt.xlabel(r'$\mathrm{Magnification\ Factor}$', fontsize = label_fontSize)
    plt.ylabel(y_label, fontsize = label_fontSize)
    
    # update tick labels to correct font
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
        
    # remove error bars in legend, format labels for power spectra
    handles, labels = ax.get_legend_handles_labels()
    handles = [h[0] for h in handles]
    if title == 'Power (1 to 4 Hz)':
        labels = [r'$\mathrm{Push\ (1-4\ Hz)}$',r'$\mathrm{Pull\ (1-4\ Hz)}$']
    if title == 'Power (4 to 7 Hz)':
        labels = [r'$\mathrm{Push\ (4-7\ Hz)}$',r'$\mathrm{Pull\ (4-7\ Hz)}$']
    if title == 'Power (7 to 10 Hz)':
        labels = [r'$\mathrm{Push\ (7-10\ Hz)}$',r'$\mathrm{Pull\ (7-10\ Hz)}$']
    leg = ax.legend(handles, labels, loc = 'upper right', 
                    fontsize = legend_fontSize, numpoints = 2, handlelength=3)
    if title == 'Contact Time':
        leg = ax.legend(handles, labels, loc = 'lower right', 
                        fontsize = legend_fontSize, numpoints = 2, 
                        handlelength=3)
    leg.get_frame().set_linewidth(0.8)
    plt.tight_layout()
    
    pp.savefig(fig)
    plt.close('all')
    pp.close()    
    
    return

def plot_LFEdata(data, title, y_label, savePath):
    
    N = 4.0
    
    mean_magOff = data[:,1]
    mean_magOn = data[:,2]
    SEM_magOff = data[:,3]/math.sqrt(N)
    SEM_magOn = data[:,4]/math.sqrt(N)
    
    target = numpy.array([-4, -3, -2, -1, 1, 2, 3, 4], dtype = 'int64')
    ticks = ['-4','-3','-2','-1','1','2','3','4']
    width = 0.05
    
    pp = PdfPages(savePath)
    
    fig = plt.figure()
    ax = plt.subplot(111)
    
    plt.errorbar(target-width, mean_magOff, yerr = SEM_magOff, color = 'k', 
                 marker = 'o', linestyle = 'None', mec = 'k', 
                 label = r'$\mathrm{Magnification\ Off}$', ms = 6.5, mew = 1.1, 
                 capsize = 0, fillstyle = 'none')
    
    plt.errorbar(target+width, mean_magOn, yerr = SEM_magOn, color = 'k',
                 marker = 'o', linestyle = 'None', ms = 5.0, mec = 'k',
                 label = 'Magnification On', mew = 1.1, capsize = 0)   
    # remove error bars only
    handles_noErr, labels = ax.get_legend_handles_labels()
    handles_noErr = [h[0] for h in handles_noErr] 
    
    if (title == 'LFE mean noFD' or title == 'LFE mean visFD'):
        # y = x target line
        plt.plot(target, target, linestyle = '--', dashes = (4,3),
                 label = r'$\mathrm{Target\ Force}$', color = 'k')
        # y = x + 0.5 and y = x - 0.5 window lines
        #plt.plot(target, target+0.5, linestyle = '-', color = 'k')
        #plt.plot(target, target-0.5, linestyle = '-', color = 'k')
        plt.fill_between(target, target+0.5, target-0.5, alpha = 0.30, 
                         color = 'y')
        plt.axhline(y = 0, color = 'k')
        plt.text(-2.8, 18, r'$\mathrm{Push}$', fontsize = text_fontSize)
        plt.text(2.2, 18, r'$\mathrm{Pull}$', fontsize = text_fontSize)
        plt.ylim(ymin = -20, ymax = 20)
    
    if (title == 'LFE std noFD' or title == 'LFE std visFD'):
        plt.text(-2.8, 11.5, r'$\mathrm{Push}$', fontsize = text_fontSize)
        plt.text(2.2, 11.5, r'$\mathrm{Pull}$', fontsize = text_fontSize)
        plt.ylim(ymin = 0, ymax = 12)
    
    plt.axvline(x = 0, color = 'k')       
    plt.grid(which = 'major', axis = 'y')
    plt.xlim(xmin = -5, xmax = 5)

    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    #ax.spines['bottom'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    ax.tick_params(axis = 'both', labelsize = tick_fontSize)
    
    plt.xticks(target, ticks, fontsize = tick_fontSize)
    plt.xlabel(r'$\mathrm{Target\ Force}\ \mathit{(grams)}$', fontsize = label_fontSize)
    plt.ylabel(y_label, fontsize = label_fontSize)
    
    # update tick labels to correct font
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    
    # remove error bars only
    handles, labels = ax.get_legend_handles_labels()
    if (title == 'LFE mean noFD' or title == 'LFE mean visFD'):
        handles_noErr = numpy.append(handles_noErr, handles[0])
        leg = ax.legend(handles_noErr, 
                    [r'$\mathrm{Magnification\ Off}$',
                     r'$\mathrm{Magnification\ On}$', 
                     r'$\mathrm{Target\ Force}$'],
                    loc = 'lower right', fontsize = legend_fontSize,
                    numpoints = 2, handlelength = 2)
    else:
        handles = [h[0] for h in handles] 
        leg = ax.legend(handles, labels, loc = 'lower right', 
                        fontsize = legend_fontSize, numpoints = 2, 
                        handlelength = 2)
    leg.get_frame().set_linewidth(0.8)
    
    plt.tight_layout()
    pp.savefig(fig)
    plt.close('all')
    pp.close()
    
    return

def plot_CDF(z_magOff, z_magOn, y_label, savePath):
    
    visFD_magOff = z_magOff[:,0]
    visFD_magOn = z_magOn[:,0]
    noFD_magOff = z_magOff[:,1]
    noFD_magOn = z_magOn[:,1]    
    
    pp = PdfPages(savePath)
    
    fig = plt.figure()
    
    N = visFD_magOff.size
    sort_magOff_visFD = numpy.sort(visFD_magOff)
    sort_magOn_visFD = numpy.sort(visFD_magOn)
    sort_magOff_noFD = numpy.sort(noFD_magOff)
    sort_magOn_noFD = numpy.sort(noFD_magOn)
    y = numpy.array(range(N))/float(N)
    
    #plot all CDF together    
    plt.plot(sort_magOff_noFD, y, color = 'k', ls = '--', dashes = (4,3), 
             label = r'$\mathrm{No\ Feedback,\ Mag\ Off}$')
             
    plt.plot(sort_magOn_noFD, y, color = 'darkgrey', ls = '--', dashes = (4,3),
             label = r'$\mathrm{No\ Feedback,\ Mag\ On}$')
             
    plt.plot(sort_magOff_visFD, y, color = 'k', 
             label = r'$\mathrm{Feedback,\ Mag\ Off}$')
             
    plt.plot(sort_magOn_visFD, y, color = 'darkgrey', 
             label = r'$\mathrm{Feedback,\ Mag\ On}$')
    
    # set up rest of figure & labels  
    ax = fig.gca()
    ax.set_axisbelow(True)
    ax.yaxis.grid(linestyle = ':', zorder = 1)
    #plt.grid(which = 'major', axis = 'y')   
    plt.ylim(ymin = 0, ymax = 1)
    plt.xlim(xmax = 3)
    #xticks = numpy.array((ax.get_xticks() * 100), dtype= numpy.int)
    #ax.set_xticklabels(xticks, size = 8) 
    #ax.set_yticklabels(ax.get_yticks(), size = tick_fontSize)
    ax.tick_params(axis = 'both', labelsize = tick_fontSize)
    plt.legend(loc = 'lower right', fontsize = legend_fontSize-3, handlelength = 3)
    plt.xlabel(r'$\mathrm{\%\ Time\ in\ Window\ Z-Score}$', fontsize = label_fontSize)
    plt.ylabel(y_label, fontsize = label_fontSize)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.yaxis.set_ticks_position('left')
    ax.xaxis.set_ticks_position('bottom')
    
    # update tick labels to correct font
    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)
    
    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)
    
#==============================================================================
#     # create custom table legend
#     box_props = dict(boxstyle = 'square', lw = 1.0)
#     table_col = [r'$\mathrm{On}$', r'$\mathrm{Off}$']
#     table_row = [r'$\mathrm{Present}$', r'$\mathrm{Absent}$']
#     table_vals = [['',''],['','']]
# #    the_table = plt.table(rowLabels = table_row,
# #                          colLabels = table_col, cellText = table_vals,
# #                          bbox = [0.70, 0.03, 0.35, 0.30], zorder = 3)
#     the_table = plt.table(rowLabels = table_row,
#                           colLabels = table_col, cellText = table_vals,
#                           loc = 'lower right', zorder = 3)
#     table_props = the_table.properties()
#     table_cells = table_props['child_artists']
#     for cell in table_cells:
#         cell.set_height(0.10)
#         cell.set_width(0.15)
#         cell.set_linewidth(0)
#     
#     the_table.set_fontsize(12)
#     # use regular text to make right hand row labels
#     col_top = r'$\mathrm{Magnification}$'
#     row_left = r'$\mathrm{Visual}$'+'\n'+'$\mathrm{Feedback}$'
#     ax.text(1.68, 0.32, col_top, fontsize = legend_fontSize, zorder = 4)
#     ax.text(0.42, 0.18, row_left, fontsize = legend_fontSize, zorder = 4,
#             rotation = 90)
#     
#     # create box to block out grid lines
#     ax.add_patch(
#         patches.Rectangle(
#             (0.35, 0.028), # (x,y) start
#             2.58,          # width  
#             0.34,          # height
#             fill = True,
#             facecolor = 'white',
#             linewidth = 0,
#             zorder = 2
#         )
#     )
#     
#     # create box around table
#     ax.add_patch(
#         patches.Rectangle(
#             (0.35, 0.028), # (x,y) start
#             2.58,          # width  
#             0.34,          # height
#             fill = False,
#             linewidth = 1,
#             zorder = 4
#         )
#     )
#     
#     # add legend markers
#     plt.axhline(y = 0.17, xmin = 0.71, xmax = 0.80, color = 'darkgrey',
#                 zorder = 5)
#     plt.axhline(y = 0.17, xmin = 0.86, xmax = 0.95, color = 'k', zorder = 5)
#     plt.axhline(y = 0.07, xmin = 0.71, xmax = 0.80, color = 'darkgrey',
#                 ls = '--', zorder = 5, dashes=(4,3))
#     plt.axhline(y = 0.07, xmin = 0.86, xmax = 0.95, color = 'k', ls = '--',
#                 dashes = (4,3), zorder = 5)
#==============================================================================
    
        
    plt.tight_layout()
    pp.savefig(fig)
    
    # plot no visual feedback condition only
    fig2 = plt.figure()

    plt.plot(sort_magOff_noFD, y, color = 'k', 
             label = 'Magnification Off')
             
    plt.plot(sort_magOn_noFD, y, color = 'k', ls = '--', dashes = (4,3),
             label = 'Magnification On') 
    
    # set up rest of figure & labels  
    ax = fig.gca()
    plt.suptitle('No Visual Feedback')
    plt.grid(which = 'major', axis = 'y')  
    plt.ylim(ymin = 0, ymax = 1)
    plt.xlim(xmax = 3)
    #xticks = numpy.array((ax.get_xticks() * 100), dtype= numpy.int)
    #ax.set_xticklabels(xticks, size = 8) 
    ax.set_yticklabels(ax.get_yticks(), size = 8)
    plt.legend(loc = 'lower right', fontsize = legend_fontSize, handlelength = 3)
    plt.xlabel(r'$\mathrm{\%\ Time\ in\ Window\ Z-Score}$', fontsize = label_fontSize)
    plt.ylabel(y_label, fontsize = label_fontSize)
    plt.tight_layout()
    pp.savefig(fig2)
    
    # plot visual feedback condition only
    fig3 = plt.figure()

    plt.plot(sort_magOff_visFD, y, color = 'k', label = 'Magnification Off')
             
    plt.plot(sort_magOn_visFD, y, color = 'k', ls = '--', dashes = (4,3),
             label = 'Magnification On') 
    
    # set up rest of figure & labels  
    ax = fig.gca()
    plt.suptitle('Visual Feedback')
    plt.grid(which = 'major', axis = 'y')  
    plt.ylim(ymin = 0, ymax = 1)
    plt.xlim(xmax = 3)
    #xticks = numpy.array((ax.get_xticks() * 100), dtype= numpy.int)
    #ax.set_xticklabels(xticks, size = 8) 
    ax.set_yticklabels(ax.get_yticks(), size = 8)
    plt.legend(loc = 'lower right', fontsize = legend_fontSize, handlelength = 3)
    plt.xlabel(r'$\mathrm{\%\ Time\ in\ Window\ Z-Score}$', fontsize = label_fontSize)
    plt.ylabel(y_label, fontsize = label_fontSize)
    plt.tight_layout()
    pp.savefig(fig3)    
    
    plt.close('all')
    pp.close()
    
    #print K-S test stats for visFD and noFD
    D_visFD, p_visFD = stats.ks_2samp(sort_magOff_visFD, sort_magOn_visFD)
    D_noFD, p_noFD = stats.ks_2samp(sort_magOff_noFD, sort_magOn_noFD)
    D_mag1_visFDnoFD, p_mag1_visFDnoFD = stats.ks_2samp(sort_magOn_visFD, 
                                                      sort_magOn_noFD)
    D_mag0_visFDnoFD, p_mag0_visFDnoFD = stats.ks_2samp(sort_magOff_visFD, 
                                                      sort_magOff_noFD)                                                      
    
    print "Kolmogorov-Smirnov tests on magnification --"
    print "D_visFD = {};\tp_visFD = {}".format(D_visFD, p_visFD)
    print "D_noFD = {};\tp_noFD = {}".format(D_noFD, p_noFD)
    print "Kolmogorov-Smirnov tests on vision --"
    print "D_mag1_visFDnoFD = {};\tp_mag1_visFDnoFD = {}".format(D_mag1_visFDnoFD, p_mag1_visFDnoFD)
    print "D_mag0_visFDnoFD = {};\tp_mag0_visFDnoFD = {}".format(D_mag0_visFDnoFD, p_mag0_visFDnoFD)
    
    return

###############################################################################
##  MAIN
###############################################################################

start = time.time()

print "\n\n*** MinContact_avgSEM_grapher.py - START *** \n\n"

os.chdir(dataFolder)

mean_rawData = numpy.genfromtxt(dataFolder+data_mean, skip_header = 1)
SD_rawData = numpy.genfromtxt(dataFolder+data_SD, skip_header = 1)
pwr14_rawData = numpy.genfromtxt(dataFolder+data_pwr14, skip_header = 1)
pwr47_rawData = numpy.genfromtxt(dataFolder+data_pwr47, skip_header = 1)
pwr710_rawData = numpy.genfromtxt(dataFolder+data_pwr710, skip_header = 1)
contact_rawData = numpy.genfromtxt(dataFolder+data_contact, skip_header = 1)

LFEmean_noFD_rawData = numpy.genfromtxt(dataFolder+data_LFEmean_noFD, 
                                        skip_header = 1)
LFEmean_visFD_rawData = numpy.genfromtxt(dataFolder+data_LFEmean_visFD,
                                         skip_header = 1)
LFEstd_noFD_rawData = numpy.genfromtxt(dataFolder+data_LFEstd_noFD,
                                       skip_header = 1)
LFEstd_visFD_rawData = numpy.genfromtxt(dataFolder+data_LFEstd_visFD,
                                        skip_header = 1)   

LFEz_mag0 = numpy.genfromtxt(dataFolder+data_LFEz_mag0, skip_header = 1)
LFEz_mag1 = numpy.genfromtxt(dataFolder+data_LFEz_mag1, skip_header = 1)

savePath_mean = outputFolder + pdfSaveName_mean
savePath_SD = outputFolder + pdfSaveName_SD
savePath_pwr14 = outputFolder + pdfSaveName_pwr14
savePath_pwr47 = outputFolder + pdfSaveName_pwr47
savePath_pwr710 = outputFolder + pdfSaveName_pwr710
savePath_power = outputFolder + pdfSaveName_power
savePath_contact = outputFolder + pdfSaveName_contact

savePath_LFEmean_noFD = outputFolder + pdfSaveName_LFEmean_noFD
savePath_LFEmean_visFD = outputFolder + pdfSaveName_LFEmean_visFD
savePath_LFEstd_noFD = outputFolder + pdfSaveName_LFEstd_noFD
savePath_LFEstd_visFD = outputFolder + pdfSaveName_LFEstd_visFD

savePath_LFEz_CDF = outputFolder + pdfSaveName_LFEz_CDF

#plot data
plot_data(mean_rawData,' ', r'$\mathrm{Mean\ Contact\ Force}\ \mathit{(grams)}$',
          savePath_mean)

plot_data(SD_rawData,' ', r'$\mathrm{SD\ Contact\ Force}\ \mathit{(grams)}$',
          savePath_SD)

plot_data(pwr14_rawData, 'Power (1 to 4 Hz)', 
          r'$\mathrm{Power}\ \mathit{(grams^2 \; Hz)}$', savePath_pwr14)

plot_data(pwr47_rawData, 'Power (4 to 7 Hz)',
          r'$\mathrm{Power}\ \mathit{(grams^2 \; Hz)}$', savePath_pwr47)

plot_data(pwr710_rawData, 'Power (7 to 10 Hz)', 
          r'$\mathrm{Power}\ \mathit{(grams^2 \; Hz)}$', savePath_pwr710)

plot_data(contact_rawData, 'Contact Time', 
          r'$\mathrm{Contact\ \ Time}\ \mathit{(sec)}$', 
          savePath_contact)

plot_LFEdata(LFEmean_noFD_rawData, 'LFE mean noFD', 
             r'$\mathrm{Mean\ Applied\ Force}\ \mathit{(grams)}$', 
             savePath_LFEmean_noFD)

plot_LFEdata(LFEmean_visFD_rawData, 'LFE mean visFD',
             r'$\mathrm{Mean\ Applied\ Force}\ \mathit{(grams)}$', 
             savePath_LFEmean_visFD)

plot_LFEdata(LFEstd_noFD_rawData, 'LFE std noFD',
             r'$\mathrm{SD\ Applied\ Force}\ \mathit{(grams)}$', 
             savePath_LFEstd_noFD)

plot_LFEdata(LFEstd_visFD_rawData, 'LFE std visFD',
             r'$\mathrm{SD\ Applied\ Force}\ \mathit{(grams)}$',
             savePath_LFEstd_visFD)          

plot_CDF(LFEz_mag0, LFEz_mag1,
         r'$\mathrm{Cumulative\ Probability}$', savePath_LFEz_CDF)


print "\nGraphing done!!"
print "Time elapsed: {0:03f} seconds".format(time.time() - start)
print "\n\n*** MinContact_contactTimeGrapher.py - END ***\n\n"
