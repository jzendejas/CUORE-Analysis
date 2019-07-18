#!/usr/bin/env python
# coding: utf-8

# In[1]:


#import ROOT stuff
from ROOT import TF1, iostream, fstream, vector, TStyle, TMath, TH1D, TFile, gInterpreter
from ROOT import std, TBranch, TTree, TChain
#import numpy stuff
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy as sp
import scipy.signal as sig
import fnmatch
import os
import sys
#get_ipython().run_line_magic('matplotlib', 'inline')


# In[ ]:


def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]


def get_files(dataPath, key):
    print("The data path is: " + dataPath)
    pathdir = os.listdir(dataPath)
    file_array = []
    for file in pathdir:
        if (file.endswith(".root") and key in file):
			file_array.append(dataPath + file)
		   
    shortFilename = write_shortFilename(file_array[0])
    return file_array, shortFilename

    
def get_fullwave(filename, dataPath=None, num_of_files=None):
    filelist, shortFilename = get_files(dataPath, filename) 

    if num_of_files == None:
        print "Warning: Attempting to Merge all "+str(filename)+" files", '\n'
        print("Proceed with Caution")
        fulldata, sampling_rates = stitch_data(filelist)
        
    elif type(num_of_files) == int:
        fulldata, sampling_rates = stitch_data(filelist[0:num_of_files])
        
    elif type(num_of_files) == float:
        print("Warning: 'num_of_files' Given as float...", '\n', "..Interpreting as Integer")
        num_of_files = int(num_of_files)
        fulldata, sampling_rates = stitch_data(filelist[0:num_of_files])
        
    
    else:
        print("Invalid 'num_of_files' given")
        return

    return fulldata, sampling_rates, shortFilename
        

def get_plot(fulldata, sampling_rates, path, shortFilename, calculation=None):
    plt.figure()
    
    for j in range(fulldata["num_of_channels"]-1):
        plot_title, txt, my_file = plot_waves(fulldata, sampling_rates, j, path, shortFilename, calculation)
        
    plt.title(plot_title)
    plt.legend(loc=0)
    plt.figtext(0.5, -0.034, txt, wrap=True, horizontalalignment='center', fontsize=12)
    plt.grid()
    plt.savefig(my_file)
    plt.show()
    
    plot_mic(fulldata, sampling_rates, path, shortFilename, calculation)
        

def get_psd(fulldata, sampling_rates, windows_to_average=None, window_size=None):
    if windows_to_average == None and window_size == None:
		fulldata = psd_no_avg(fulldata, sampling_rates)
	
    else:
	    fulldata = psd_avg(fulldata, sampling_rates, windows_to_average, window_size)
	
    return fulldata
    

def grab_data(filename):
    myfile = TFile(str(filename), 'read')
    tree = myfile.Get("data_tree")
    tree.GetEvent()

    Branch_List = [tree.GetListOfBranches()[i].GetName() for i in range(len(tree.GetListOfBranches()))]
    nSamples = getattr(tree, Branch_List[0])
    dt = getattr(tree, Branch_List[3])

    nEntries  = tree.GetEntriesFast()
    nTotal = nSamples*nEntries
    sampling_rate = int(1/dt)

    data = {}
    data["times"] = np.linspace(0, (nTotal-1)*dt, nTotal)

    for i,j in enumerate(Branch_List[4::]): 
        print '\r',"Getting ",j,
        data["wave"+str(i)] = (np.asarray([[getattr(tree, j)[0::]] for event in tree])).flatten()
        
    return data, sampling_rate, dt
    
    
def plot_mic(fulldata, sampling_rates, path, shortFilename, calculation):
    plt.figure()
    
    if calculation == "ts":
        plt.plot(fulldata["times"], fulldata["wave"+str(fulldata["num_of_channels"]-1)], label="Mic")
        plt.title(shortFilename+"_TimeStream")
        txt = ""
        my_file2 = str(path) + str(shortFilename) + "_micTS.png"
    elif calculation == 'psd':
        plt.loglog(fulldata["PSD_f_wave"+str(fulldata["num_of_channels"]-1)], fulldata["PSD_wave"+str(fulldata["num_of_channels"]-1)], label="Mic")
        plt.title(shortFilename+"_PSD")
        try:
			txt = 'Average of '+str(fulldata["M"])+', '+str(fulldata["N"])+' point samples, Sampling Rate: '+str(sampling_rates[0]/1000)+' kHz'
        except KeyError:
            txt = ""
        my_file2 = str(path) + str(shortFilename) + "_micPSD.png"
    elif calculation == 'fft':
        plt.loglog(fulldata["FFT_f_wave"+str(fulldata["num_of_channels"]-1)], fulldata["FFT_wave"+str(fulldata["num_of_channels"]-1)], label="Mic")
        plt.title(shortFilename+"_FFT")
        txt = 'Average of '+str(fulldata["M"])+', '+str(fulldata["N"])+' point samples'
        my_file2 = str(path) + str(shortFilename) + "_micFFT.png"
    plt.figtext(0.5, -0.034, txt, wrap=True, horizontalalignment='center', fontsize=12)
    plt.legend(loc=0)
    plt.grid()
    plt.savefig(my_file2)
    plt.show()

    
def plot_waves(fulldata, sampling_rates, j, path, shortFilename, calculation):   
    if calculation == "ts":
        plt.plot(fulldata["times"], fulldata["wave"+str(j)], label = "wave "+str(j))
        plot_title = shortFilename+"_TimeStream"
        my_file = str(path) + str(shortFilename) + "_TS.png"
        txt = ""
    elif calculation == 'psd':
        plt.loglog(fulldata["PSD_f_wave"+str(j)], fulldata["PSD_wave"+str(j)], label = "wave "+str(j))
        plot_title = shortFilename+"_PSD"
        my_file = str(path) + str(shortFilename) + "_PSD.png"
        try:
			txt = 'Average of '+str(fulldata["M"])+', '+str(fulldata["N"])+' point samples, Sampling Rate: '+str(sampling_rates[0]/1000)+' kHz'
        except KeyError:
			txt = ""
    elif calculation == 'fft':
        plt.loglog(fulldata["times"], fulldata["wave"+str(j)], label = "wave "+str(j))
        plot_title = shortFilename+"_FFT"
        my_file = str(path) + str(shortFilename) + "_FFT.png"
        txt = 'Average of '+str(fulldata["M"])+', '+str(fulldata["N"])+' point samples'
    return plot_title, txt, my_file


def psd_avg(fulldata, sampling_rates, windows_to_average, window_size):
	
	assert (type(window_size+windows_to_average) == int), "Both window_size and windows_to_average must be integers"
		
	sums = 0
	samps = {}
	
	for i in range(fulldata["num_of_channels"]):
		for j in range(windows_to_average):
			samps[j] = fulldata["wave"+str(i)][j*window_size : ((j*window_size)+(window_size-1))]
			f, samps[j] = sig.periodogram(samps[j], fs=sampling_rates[0])
			
		for k in range(windows_to_average):
			sums += samps[k]
			
		PSD_output = sums/windows_to_average
		
		fulldata["PSD_f_wave"+str(i)], fulldata["PSD_wave"+str(i)] = f, PSD_output
		fulldata["M"], fulldata["N"] = windows_to_average, window_size
		return fulldata


def psd_no_avg(fulldata, sampling_rates):
	for i in range(fulldata["num_of_channels"]):
		fulldata["PSD_f_wave"+str(i)], fulldata["PSD_wave"+str(i)] = sig.periodogram(fulldata["wave"+str(i)], fs=sampling_rates[0])
	return fulldata	

	
def stitch_data(filelist):
    #used to check the number of channels being read in from the data file
    dummy_data, dummy_sampling_rate, dummy_dt = grab_data(filelist[0])
    num_of_channels = len(dummy_data)-1
    
    fulldata = {}
    fulldata["times"] = np.array([])
    for i in range(num_of_channels):
        fulldata["wave"+str(i)] = np.array([])
        
    sampling_rates = np.array([])
    
    for i in filelist:
        data, sampling_rate, dt = grab_data(i)
        for j in range(num_of_channels):
            fulldata["wave"+str(j)] = np.append(fulldata["wave"+str(j)], data["wave"+str(j)])
        fulldata["times"] = np.append(fulldata["times"], (data["times"] +  ( len(sampling_rates) * (data["times"][-1]+dt) )   ))
        sampling_rates = np.append(sampling_rates, sampling_rate)
        if len(sampling_rates) > 1:
            if sampling_rates[-1] != sampling_rates[-2]:
                print("Warning: Files of Varying Sampling Rates are being Merged")
				
				
    fulldata["num_of_channels"] = num_of_channels
    return fulldata, sampling_rates
    

def write_shortFilename(filename):
    shortFilename = filename.split("/")[-1]
    shortFilename_index = findOccurrences(shortFilename, '_')[-2]
    shortFilename = shortFilename[0:shortFilename_index]
    return shortFilename



def main(filename, num_of_files=None, loglog=True, xlim=None, ylim=None, windows_to_average=None, window_size=None, cutoff=None):
	repoPath = os.environ['ANALYSISREPO']
	dataPath = repoPath + "/data/root/"
	imgPath = repoPath + "/images/"
	
	fulldata, sampling_rates, shortFilename = get_fullwave(filename, dataPath, num_of_files)
	fulldata = get_psd(fulldata, sampling_rates, windows_to_average, window_size)
	get_plot(fulldata, sampling_rates, imgPath, shortFilename, calculation="ts")
	get_plot(fulldata, sampling_rates, imgPath, shortFilename, calculation="psd")
	#get_FFT(data, imgPath, window_size, windows_to_average, cutoff, xlim, ylim)
	
# ## Execute Main

argc = len(sys.argv)
if argc == 2:
    key = str(sys.argv[1])
    main(key)
elif argc == 3:
	key = str(sys.argv[1])
	num = int(sys.argv[2])
	main(key, num_of_files = num)
else: 
    print("Run as python Anal.py <key> ")

