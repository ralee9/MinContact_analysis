# MinContact_analysis

A repository of assorted Python scripts to analyze voltage data from a psychophysics experiment.

Data are collected at 1 kHz and stored as rows, each voltage source as its own column.

Run MinContact _RawDataProcessing.py to split data and create graphs by trial
Run MinContact_DataStatistics.py to obtain mean, standard deviation, power, and contact time statistics

Run LowForceExp_RawDataProcessing.py to split Low Force Matching data and create graphs by trial 
Run LowForceExp_DataStatistics.py to obtain mean, standard deviation, power statistics

Scripts MinContact_avgSEM_grapher.py and MinContact_avgPower_grapher.py use matplotlib to create figures for 
publication
Script MinContact_welchTest.py calculates the Welch's t-test statistics for contact
