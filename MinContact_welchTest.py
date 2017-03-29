# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 18:11:02 2017

@author: Randy Lee
"""
##############################################################################
## Randy Lee
## VIA Lab
## Department of Bioengineering
## University of Pittsburgh
##
## Calculate Welch's t-test on column data
##
##############################################################################

import os
import numpy
import time
from scipy import stats
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# edit these constants to access figure
rootFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts'
dataFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts\data'
outputFolder = 'C:\Users\Randy Lee\Documents\VIA Lab\Python Scripts\output'
#rootFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts'
#dataFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts\data'
#outputFolder = 'D:\Randy Lee\Documents\VIA Lab\Python Scripts\output'

data_file = '\\MinContact_contactTime_all.txt'

def calc_stats(data):
    y1 = data[:,1]
    y2 = data[:,2]
    y3 = data[:,3]
    
    num1 = y1.shape[0]
    num2 = y2.shape[0]
    num3 = y3.shape[0]
    var1 = numpy.var(y1, ddof=1)
    var2 = numpy.var(y2, ddof=1)
    var3 = numpy.var(y3, ddof=1)
    print "var1 = {}; var2 = {}; var3 = {}\n".format(var1, var2, var3)
    
    #1 vs 2
    #curr_t, curr_p = stats.ttest_ind(y1, y2, equal_var=False)
    curr_t, curr_p = stats.ttest_rel(y1, y2)
    df1 = ((var1/num1 + var2/num2)**(2.0))/((var1/num1)**(2.0)/(num1-1) + (var2/num2)**(2.0)/(num2-1))
    print "**Testing y1 vs y2\nt = {}; p = {}".format(curr_t, curr_p)
    #print "Welch's df = {}".format(df1)
    print "df = {}".format(num1/2. - 1.)
    
    #1 vs 3
    #curr_t, curr_p = stats.ttest_ind(y1, y3, equal_var=False
    curr_t, curr_p = stats.ttest_rel(y1, y3)
    df2 = ((var1/num1 + var3/num3)**(2.0))/((var1/num1)**(2.0)/(num1-1) + (var3/num3)**(2.0)/(num3-1))
    print "\n\n**Testing y1 vs y3\nt = {}; p = {}".format(curr_t, curr_p)
    #print "Welch's df = {}".format(df2)
    print "df = {}".format(num1/2. - 1.)
    
    #2 vs 3
    #curr_t, curr_p = stats.ttest_ind(y2, y3, equal_var=False)
    curr_t, curr_p = stats.ttest_rel(y2, y3)
    df3 = ((var2/num2 + var3/num3)**(2.0))/((var2/num2)**(2.0)/(num2-1) + (var3/num3)**(2.0)/(num3-1))
    print "\n\n**Testing y2 vs y3\nt = {}; p = {}".format(curr_t, curr_p)
    #print "Welch's df = {}".format(df3)
    print "df = {}".format(num1/2. - 1.)
    
    return
    
##############################################################################
## MAIN
##############################################################################
if __name__ == '__main__':
    start = time.time()
    
    print "\n\n*** MinContact_welchTest.py - START ***\n\n"
    
    os.chdir(dataFolder)
    
    data = numpy.genfromtxt(dataFolder+data_file, skip_header=1)
    
    #savePath = outputFolder + save_file
    
    calc_stats(data)
    
    print "\n\n*** MinContact_welchTest.py - END ***"
    print "Time elapsed: {0:03f} seconds".format(time.time() - start)
    
