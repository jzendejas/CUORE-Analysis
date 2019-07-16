#!/usr/bin/env python
# coding: utf-8

# In[228]:


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


# ## Helper Functions:

# ### Helper Functions: Getting Data

# #### Helper Functions: Getting Data/File Searching

# In[388]:

def findOccurrences(s, ch):
    return [i for i, letter in enumerate(s) if letter == ch]

def write_shortFilename(filename):
    shortFilename = filename.split("/")[-1]
    shortFilename_index = findOccurrences(shortFilename, '_')[-2]
    shortFilename = shortFilename[0:shortFilename_index]
    return shortFilename
	
	
def get_files(dataPath, key):
    print("The data path is: " + dataPath)
    pathdir = os.listdir(dataPath)
    file_array = []
    for file in pathdir:
        if (file.endswith(".root") and key in file):
			file_array.append(dataPath + file)
		   
    shortFilename = write_shortFilename(file_array[0])
    return file_array, shortFilename


# #### Helper Functions: Getting Data/Data Extraction

# In[230]:


def grab_data(filename):
    myfile = TFile(str(filename), 'read')
    tree = myfile.Get("data_tree")
    tree.GetEvent()
    nSamples = tree.NumberOfSamples
    nEntries  = tree.GetEntriesFast()
    nTotal = nSamples*nEntries
    dt = tree.SamplingWidth_s
    sampling_rate = int(1/dt)

    times = np.linspace(0, nTotal*dt-dt, nTotal)
    wave0 = (np.asarray([[tree.Waveform000[0:sampling_rate]] for event in tree])).flatten()
    wave1 = (np.asarray([[tree.Waveform001[0:sampling_rate]] for event in tree])).flatten()
    wave2 = (np.asarray([[tree.Waveform002[0:sampling_rate]] for event in tree])).flatten()
    
    return times, wave0, wave1, wave2, sampling_rate, dt

def stitch_data(filelist):
    fullwave0 = np.array([])
    fullwave1 = np.array([])
    fullwave2 = np.array([])
    fulltime = np.array([])
    sampling_rates = np.array([])
    
    for i in filelist:
        times, wave0, wave1, wave2, sampling_rate, dt = grab_data(i)
        fullwave0 = np.append(fullwave0, wave0)
        fullwave1 = np.append(fullwave1, wave1)
        fullwave2 = np.append(fullwave2, wave2)
        fulltime  = np.append(fulltime, (   times + ( len(sampling_rates) * (times[-1]+dt) )  )  )
        sampling_rates = np.append(sampling_rates, sampling_rate)
        if len(sampling_rates) > 1:
            if sampling_rates[-1] != sampling_rates[-2]:
                print("Warning: Files of Varying Sampling Rates are being Merged")
                #return
    return fulltime, fullwave0, fullwave1, fullwave2, sampling_rates

def get_fullwave(filename, dataPath, num_of_files=None):
    if num_of_files == None:
        print "Warning: Attempting to Merge all "+str(filename)+" files", '\n'
        print("Proceed with Caution")
        filelist, shortFilename = get_files(dataPath, filename)
        fulltime, fullwave0, fullwave1, fullwave2, sampling_rates = stitch_data(filelist)
        return fulltime, fullwave0, fullwave1, fullwave2, sampling_rates, shortFilename
        
    elif type(num_of_files) == int:
        filelist, shortFilename = get_files(dataPath, filename)
        fulltime, fullwave0, fullwave1, fullwave2, sampling_rates = stitch_data(filelist[0:num_of_files])
        return fulltime, fullwave0, fullwave1, fullwave2, sampling_rates, shortFilename
        
    elif type(num_of_files) == float:
        print("Warning: 'num_of_files' Given as float...", '\n', "..Interpreting as Integer")
        num_of_files = int(num_of_files)
        filelist, shortFilename = get_files(dataPath, filename)
        fulltime, fullwave0, fullwave1, fullwave2, sampling_rates = stitch_data(filelist[0:num_of_files])
        return fulltime, fullwave0, fullwave1, fullwave2, sampling_rates, shortFilename
    
    else:
        print("Invalid 'num_of_files' given")
        return
        
    


# ### Helper Functions: Plotting Assistance

# In[231]:



def find_max(hx, hy, hz, vx, vy, vz ):
    
    horizontal_axes = np.array([hx, hy, hz])
    vertical_axes = np.array([vx, vy, vz])
    
    h_maxes = np.array([max(hx), max(hy), max(hz)])
    h_max_ax = horizontal_axes[np.where(h_maxes == max(h_maxes))[0][0]]
    
    v_maxes = np.array([max(vx), max(vy), max(vz)])
    v_max_ax = vertical_axes[np.where(v_maxes == max(v_maxes))[0][0]]
    return h_max_ax, v_max_ax


def find_min(vx, vy, vz):
    vertical_axes = np.array([vx, vy, vz])
    v_mins = np.array([min(vx), min(vy), min(vz)])
    v_min_ax = vertical_axes[np.where(v_mins == max(v_mins))[0][0]]
    
    return v_min_ax

def axis_limits(x, y_max, y_min, xlim=None, ylim=None, psd=False):
    
    if xlim == None:
        xlim_index = np.where(y_max >= 0.00001*max(y_max))[0][-1]
        try:
            xlim = x[xlim_index+10]
            plt.xlim(0, xlim)
        except IndexError:
            xlim = x[xlim_index]
            plt.xlim(0, xlim)
    elif type(xlim) == int:
        plt.xlim(0,xlim)
    elif type(xlim) == float:
        plt.xlim(0,xlim)
    
    
    
    if ylim == None:
        if psd==False:
            ylim_max = max(y_max) + 0.15*max(y_max)
            ylim_min = min(y_min[0:xlim_index]) - 0.15*min(y_min[0:xlim_index])
            #ylim_min = y_min - 0.001*y_min
        if psd==True:
            ylim_max = max(y_max) + 0.25*max(y_max)
            ylim_min = min(y_min[1:xlim_index]) - 0.15*min(y_min[0:xlim_index])
            #ylim_min = y_min - 0.001*y_min
        plt.ylim(ylim_min, ylim_max)
    elif type(ylim) == int:
        plt.ylim(0,ylim)
    elif type(ylim) == float:
        plt.ylim(0,ylim)


# ### Helper Functions: Data Analysis:

# #### Helper Functions: Data Analysis/TimeStream

# In[385]:


def get_timestream(data, path):
    fulltime, fullwave0, fullwave1, fullwave2, sampling_rates, shortFilename = data
    
    my_file = str(path) + str(shortFilename) + "_Timestream.png"
    
    plt.figure()
    plt.plot(fulltime, fullwave0, color = 'r', label = "x")
    plt.plot(fulltime, fullwave1, color = 'b', label = "y")
    plt.plot(fulltime, fullwave2, color = 'g', label = "z")
    plt.xlabel("Time [Seconds]")
    plt.ylabel("Amplitude [mV ?]")
    plt.legend(loc=0)
    plt.title(str(shortFilename))
    plt.grid()
    
    plt.savefig(my_file)
    
    #plt.show()
    


# #### Helper Functions: Data Analysis/FFT

# In[383]:


def gauss_wind(N, sigma):
    M = N - 2 
    w = np.zeros(M+1)
    for i in range(M+1):
        w[i] = np.exp(-0.5*(((i - M)/2)/(sigma*M/2))*(((i - M)/2)/(sigma*M/2)))
    return w
        
def butter_lowpass(cutoff, fs, order=5):
    #fs sample rate in Hz
    #desired cutoff frequency of filter in Hz
    nyq = 0.5*fs
    normal_cutoff = cutoff / nyq
    b, a = sp.signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    #fs sample rate in Hz
    #desired cutoff frequency of filter in Hz
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = sp.signal.lfilter(b, a, data)
    return y

def plot_fft(FFT_Output, M, N, shortFilename, xlim, ylim, path):
    xf, yfx, yfy, yfz = FFT_Output
    
    my_file =  str(path) + str(shortFilename) + "_FFT.png"
	
    plt.figure()
    plt.loglog(xf, yfx, color = 'r', label = 'X; Gaussian Window' )
    plt.loglog(xf, yfy, color = 'b', label = 'Y; Gaussian Window' )
    plt.loglog(xf, yfz, color = 'g', label = 'Z; Gaussian Window' )
    plt.legend(loc=0)
    plt.title("FFT of File: "+str(shortFilename))
    
    xf_max, yf_max = find_max(xf, xf, xf, yfx, yfy, yfz)
    yf_min = find_min(yfx, yfy, yfz)
    axis_limits(xf_max, yf_max, yf_min, xlim, ylim)
    
    txt = 'Average of '+str(M)+', '+str(N)+' point samples'
    plt.figtext(0.5, -0.034, txt, wrap=True, horizontalalignment='center', fontsize=12)
    
    
    plt.grid()
    plt.savefig(my_file)
    
    #plt.show()


    
    
def calulate_FFT(windows_to_average, window_size, sampling_rate, cutoff, micG):
    sumsG = 0 
    sampsG = {}
    
    for j in range(windows_to_average):
        sampsG[j] = micG[j*window_size : ((j*window_size)+(window_size-1))]
        if type(cutoff) != type(None):
            sampsG[j] = butter_lowpass_filter(sampsG[j], cutoff, fs=sampling_rate, order=4)
        np.multiply(sampsG[j], gauss_wind(window_size, 0.5), out = sampsG[j], casting = 'unsafe')
        sampsG[j] = np.fft.fft(sampsG[j])

    
    for k in range(windows_to_average):
        sumsG += sampsG[k]

       
    FFT_output = sumsG/windows_to_average
    yfG = 2./window_size*np.absolute(FFT_output[1:window_size//2])
    
    return yfG
    
    

def get_FFT(data, path, window_size=None, windows_to_average=None, cutoff=None, xlim=None, ylim=None ):
    
    fulltime, fullwave0, fullwave1, fullwave2, sampling_rates, shortFilename = data
    sampling_rate = sampling_rates[0]
    
    if type(window_size) == type(None):
        window_size = len(fulltime)
    
    if type(windows_to_average) == type(None):
        windows_to_average = 1
    
    xf = np.linspace(0, 1.0/(2.0/sampling_rate), window_size//2-1)
    
    
    
    yx = np.copy(fullwave0)
    yy = np.copy(fullwave1)
    yz = np.copy(fullwave2)

    yfx = calulate_FFT(windows_to_average, window_size, sampling_rate, cutoff, yx)
    yfy = calulate_FFT(windows_to_average, window_size, sampling_rate, cutoff, yy)
    yfz = calulate_FFT(windows_to_average, window_size, sampling_rate, cutoff, yz)
    
    
    FFT_Output = (xf, yfx, yfy, yfz)
    
    plot_fft(FFT_Output, windows_to_average, window_size, shortFilename, xlim, ylim, path)

    

# #### Helper Functions: Data Analysis/PSD

# In[379]:


def calc_psd(data, windows_to_average, window_size ):
    fulltime, fullwave0, fullwave1, fullwave2, sampling_rates, shortFilename = data
    sampling_rate = sampling_rates[0]
    fx, px = psd_averaging(windows_to_average, window_size, sampling_rate, fullwave0)
    fy, py = psd_averaging(windows_to_average, window_size, sampling_rate, fullwave1)
    fz, pz = psd_averaging(windows_to_average, window_size, sampling_rate, fullwave2)
    
    return fx, fy, fz, px, py, pz, sampling_rate


def psd_averaging(windows_to_average, window_size, sampling_rate, micG):
    sumsG = 0 
    sampsG = {}
    for j in range(windows_to_average):
        sampsG[j] = micG[j*window_size : ((j*window_size)+(window_size-1))]
        f, sampsG[j] = sig.periodogram(sampsG[j], fs=sampling_rate)


    for k in range(windows_to_average):
        sumsG += sampsG[k]

        
    PSD_output = sumsG/windows_to_average
    
    
    return f, PSD_output


def plot_psd(psd_data, shortFilename, path, loglog, xlim, ylim, windows_to_average, window_size):
    fx, fy, fz, px, py, pz, sampling_rate = psd_data
    
    my_file = str(path) + str(shortFilename) + "_PSD.png"
    plt.figure()
    
    if loglog:
        plt.loglog(fx, px, color = 'r', label="x")
        plt.loglog(fy, py, color = 'b', label="y")
        plt.loglog(fz, pz, color = 'g', label="z")
    else: 
        plt.plot(fx, px, color = 'r', label="x")
        plt.plot(fy, py, color = 'b', label="y")
        plt.plot(fz, pz, color = 'g', label="z")
        
        
    plt.xlabel("Frequency [Hz]")
    plt.ylabel("Amplitude [V^2/Hz]")
    plt.title("PSD of: "+str(shortFilename))
    plt.legend(loc=0)
    plt.grid()

    f_max, p_max = find_max(fx,fy,fz,px,py,pz)
    p_min = find_min(px,py,pz)
    axis_limits(f_max, p_max, p_min, xlim, ylim, psd=True)
    
    
    text = """Average of %s, %s point samples \nSampling Rate: %s kHz"""
    plt.figtext(0.05,0.00, text %(windows_to_average, window_size, sampling_rate/1000), fontsize=8, va="top", ha="left")
    plt.savefig(my_file)
    #plt.show()
    
    
    
def get_psd(data, path, loglog=True, xlim=None, ylim=None, windows_to_average=None, window_size=None):
    #input file name as string (ex: "file")
    
    if type(window_size) == type(None):
        window_size = len(data[0])
    
    if type(windows_to_average) == type(None):
        windows_to_average = 1
    
    psd_data = calc_psd(data, windows_to_average, window_size)
    plot_psd(psd_data, data[-1], path, loglog, xlim, ylim, windows_to_average, window_size)


# ## Main Function

# In[389]:


def main(filename, num_of_files=None, loglog=True, xlim=None, ylim=None, windows_to_average=None, window_size=None, cutoff=None):
    repoPath = os.environ['ANALYSISREPO']
    dataPath = repoPath + "/data/root/"
    imgPath = repoPath + "/images/"
	
	
    data = get_fullwave(filename, dataPath, num_of_files)
    #paths = set_paths(path)
    get_timestream(data, path=imgPath)
    get_psd(data, imgPath, loglog, xlim, ylim, windows_to_average, window_size)
    get_FFT(data, imgPath, window_size, windows_to_average, cutoff, xlim, ylim)


# ## Execute Main

argc = len(sys.argv)
if argc == 2:
    key = str(sys.argv[1])
    main(key)

else: 
    print("Run as python Anal.py <key> ")



