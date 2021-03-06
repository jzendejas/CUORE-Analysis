#!/usr/bin/env python
# coding: utf-8
# 0:accx,1:accy,2:accz,3:mic1
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

###New: Plotting User Prompts
def request_plotting_setting():
    confirm = "n"
    while confirm != "y":
        ps = input("Would you like custom or automatic plotting? (custom/auto)" )
        if ps == "automatic":
            ps = "auto"
        confirm = input("Are you sure? (y/n)")
    return ps
	
	
def custom_plots_settings(ps):
    confirm = "n"
    while confirm != "y":
        print("Initializing custom plots..")
        num_of_plots = input("How many plots would you like to generate?")
        plot_channel_map = {}
        plot_title_map = {}
        assert type(num_of_plots) == int, "The Number of Plots to Generate must be an integer"
        for i in range(num_of_plots):
            plot_channel_map["plot"+str(i+1)] = input("Plot "+str(i+1)+" : Please enter devices to be plotted (accx,mic1,TES1,etc)")
            plot_title["plot"+str(i+1)] = input("Plot "+str(i+1)+" : Please enter general plot title (Accelerometers)")
        print("You have entered the following:", '\n')
        for i in range(num_of_plots):
            print("Plot "+str(i+1))
            print("Title")
            print(plot_title_map["plot"+str(i+1)])
            print("Devices")
            print(plot_channel_map["plot"+str(i+1)])
        confirm = input("Is this correct? (y/n)")
		
    plot_channel_map = separate_plot_map(plot_channel_map)
    return plot_channel_map, plot_title_map
	
def separate_plot_map(plot_channel_map):
	count = 0 
	map_2 = {}
	for j in plot_channel_map:
		commas = findOccurrences(plot_channel_map[j], ",")
		map_2[j] = []
		for i in commas:
			if count == 0:
				map_2[j].append(plot_channel_map[j][:i])
				count += 1
				i_prev = i
			else:
				map_2[j].append(plot_channel_map[j][i_prev+1:i])
				i_prev = i
				count += 1
			if i == commas[-1]:
				map_2[j].append(plot_channel_map[j][i+1::])
	return map_2


###Helper Functions: Used in user-requested Channel labeling and object instancing.
def print_channel_list(spaces, channel_list):
	count = 0
	for i in spaces:
		if count == 0:
			print("channel 0: "+str(channel_list[2:i-1])+" "+str(channel_list[i-1]))
			count +=1
			i_prev = i
		else:
			print("channel "+str(count)+": "+str(channel_list[i_prev+3:i-1])+" "+str(channel_list[i-1]))
			i_prev = i
			count += 1
			
		if i == spaces[-1]:
			print("channel "+str(count)+": "+str(channel_list[i+3:-1])+" "+str(channel_list[-1]))
	return
def confirm_channel_list():
	confirmed = False
	while confirmed == False:
		channel_list = input("Please enter channel number with associated device name (ex: 0:accx,1:mic1,2:TES1  etc.)")
		spaces = findOccurrences(channel_list, ",")
		print_channel_list(spaces, channel_list)
		confirm = input("Confirm? (y/n)")
		if confirm == "y":
			confirmed = True
		elif confirm == "n":
			pass
		else:
			print("Input not understood...")
			pass

	return spaces, channel_list
def write_channel_list(spaces, channel_list):
	count = 0
	channel_dict = {}
	
	for i in spaces:
		if count == 0:
			channel_dict[count] = (channel_list[2:i], channel_list[2:i-1])
			count +=1
			i_prev = i
		else:
			channel_dict[count] = (channel_list[i_prev+3:i], channel_list[i_prev+3:i-1])
			i_prev = i
			count += 1
			
		if i == spaces[-1]:
			channel_dict[count] = (channel_list[i+3::], channel_list[i+3:-1])
	print(channel_dict)
	return channel_dict

	
##Plotting Function (and its sub functions)
###Find best axis limits
def simple_search(a, max_or_min):
    if max_or_min == "max":
        return max(a)
    elif max_or_min == "min":
        return min(a)
    else:
        print("Please specify max or min.")
        return
def find_best_lims(channel_objs, calculation):
    uylims = []
    lylims = []
    uxlims = []
    lxlims = []
    for i in channel_objs:
        lims = i.axis_limits(calculation)
        uylims.append(lims[0])
        lylims.append(lims[1])
        uxlims.append(lims[2])
        lxlims.append(lims[3])
    Uylim = simple_search(uylims, max_or_min="max")
    Lylim = simple_search(lylims, max_or_min="min")
    Uxlim = simple_search(uxlims, max_or_min="max")
    Lxlim = simple_search(lxlims, max_or_min="min")
    Xlimits = [Lxlim, Uxlim]
    Ylimits = [Lylim, Uylim]
    return Xlimits, Ylimits


###Actual Plotting of the function
def hp_gen(calculation, channel_objs, shortFilename, path):
    if calculation == "ts":
        plotcalc = "plot"
    elif calculation == "fft" or calculation == "psd":
        plotcalc = "loglog"
    #getattr(plt, plotcalc)()
    plt.figure()
    plot_title = shortFilename+"_"+str(calculation)
    plt.title(plot_title)
    x_attribute = "x_"+str(calculation)
    y_attribute = "y_"+str(calculation)
    temp_dev_objs = []
    for i,j in enumerate(channel_objs):
        if i == 0:
            prev_dev = getattr(j, "device")
            getattr(plt, plotcalc)(getattr(j, x_attribute), getattr(j, y_attribute), label=str(getattr(j, "name")))
            temp_dev_objs.append(j)
        else: 
            dev = getattr(j, "device")
            if prev_dev == dev:
                getattr(plt, plotcalc)(getattr(j, x_attribute), getattr(j, y_attribute), label=str(getattr(j, "name")))
                prev_dev = getattr(j, "device")
                temp_dev_objs.append(j)
            else:
                plt.legend(loc=0)
                my_file = str(path)+str(shortFilename)+str(prev_dev)+"_"+str(calculation)
                Xlimits, Ylimits = find_best_lims(temp_dev_objs, calculation)
                plt.xlim(Xlimits)
                plt.ylim(Ylimits)
                plt.savefig(my_file)
                plt.show()
                temp_dev_objs = []
                plt.figure()
                getattr(plt, plotcalc)(getattr(j, x_attribute), getattr(j, y_attribute), label=str(getattr(j, "name")))
                prev_dev = getattr(j, "device")
                temp_dev_objs.append(j)
    plt.legend(loc=0)
    my_file = str(path)+str(shortFilename)+str(dev)+"_"+str(calculation)
    Xlimits, Ylimits = find_best_lims(temp_dev_objs, calculation)
    plt.xlim(Xlimits)
    plt.ylim(Ylimits)
    plt.savefig(my_file)
    plt.show()
    return	
def hp_custom(calculation, channel_objs, plot_channel_map, plot_title_map, num_of_plots, shortFilename, path):	
    if calculation == "ts":
        plotcalc = "plot"
    elif calculation == "fft" or calculation == "psd":
        plotcalc = "loglog"
		
    x_attribute = "x_"+str(calculation)
    y_attribute = "y_"+str(calculation)
	
    for i,l in enumerate(plot_title_map):
        plt.figure()
        plt.tile(plot_title_map[l])
        temp_plot_objs = []
        for j, k in enumerate(plot_channel_map["plot"+str(i+1)]):
            channel_name = k
            obje = compare_names(channel_objs, channel_name)
            temp_plot_objs.append(obje)
            getattr(plt, plotcalc)(getattr(obje, x_attribute), getattr(obje, y_attribute), label=str(k))
        plt.legend(loc=0)
        my_file = str(path)+str(shortFilename)+str(plot_title_map[l])+"_"+str(calculation)
        Xlimits, Ylimits = find_best_lims(temp_plot_objs, calculation)
        plt.xlim(Xlimits)
        plt.ylim(Ylimits)
        plt.savefig(my_file)
        plt.show()
		
		
def compare_names(channel_objs, channel_name):
	for i in channel_objs:
		if getattr(i, name) == channel_name:
			return i
		else:
			pass 
###Quoted functions have (hopefully) been replaced with hp_gen	
"""	
def hp_TS(calculation, channel_objs, shortFilename, path):
	plt.figure()
	plot_title = shortFilename+"_Timestream"
	plt.title(plot_title)
	for i,j in enumerate(channel_objs):
		if i == 0:
			prev_dev = getattr(j, device)
			plt.plot(getattr(j, time), getattr(j, wave), label=str(getattr(j, name)))
		else: 
			dev = getattr(j, device)
			if prev_dev == dev:
				plt.plot(getattr(j, time), getattr(j, wave), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
			else:
				plt.legend(loc=0)
				my_file = str(path)+str(shortFilename)+str(prev_dev)+"_TS"
				plt.savefig(my_file)
				plt.show()
				plt.figure()
				plt.plot(getattr(j, time), getattr(j, wave), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
	plt.legend(loc=0)
	my_file = str(path)+str(shortFilename)+str(dev)+"_TS"
	plt.savefig(my_file)
	plt.show()
	return	
def hp_FFT(calculation, channel_objs, shortFilename, path):
    plt.figure()
    plot_title = shortFilename+"_FFT"
    plt.title(plot_title)
    for i,j in enumerate(channel_objs):
        if i == 0:
			prev_dev = getattr(j, device)
			plt.plot(getattr(wave+str(i), freqsf), getattr(j, fft), label=str(getattr(j, name)))
        else: 
			dev = getattr(j, device)
            if prev_dev == dev:
				plt.plot(getattr(j, freqsf), getattr(j, fft), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
            else:
				plt.legend(loc=0)
				my_file = str(path)+str(shortFilename)+str(prev_dev)+"_FFT"
				plt.savefig(my_file)
				plt.show()
				plt.figure()
				plt.plot(getattr(j, freqsf), getattr(j, fft), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
    plt.legend(loc=0)
    my_file = str(path)+str(shortFilename)+str(dev)+"_FFT"
    plt.savefig(my_file)
    plt.show()
    return
def hp_PSD(calculation, channel_objs, shortFilename, path):
	plt.figure()
	plot_title = shortFilename+"_PSD"
	plt.title(plot_title)
	for i,j in enumerate(channel_objs):
		if i == 0:
			prev_dev = getattr(j, device)
			plt.plot(getattr(wave+str(i), freqp), getattr(j, psd), label=str(getattr(j, name)))
		else: 
			dev = getattr(j, device)
			if prev_dev == dev:
				plt.plot(getattr(j, freqp), getattr(j, psd), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
			else:
				plt.legend(loc=0)
				my_file = str(path)+str(shortFilename)+str(prev_dev)+"_PSD"
				plt.savefig(my_file)
				plt.show()
				plt.figure()
				plt.plot(getattr(j, freqp), getattr(j, psd), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
				
	plt.legend(loc=0)
	my_file = str(path)+str(shortFilename)+str(dev)+"_PSD"
	plt.savefig(my_file)
	plt.show()
	return
"""		
def hyper_plots(calculation, channel_objs, shortFilename, path):
    if calculation == "ts" or calculation == "fft" or calculation == "psd":
        hp_gen(calculation, channel_objs, shortFilename, path)
    else:
        print("Input not Understood: Please use Keys; 'ts', 'fft', 'psd'")
    return

	
###Simple Notch filter functions	
def notch_filter(yf, f0, fs, Q):
    b, a = sig.irrnotch(f0, Q, fs)
    freq, h = sig.freqz(b, a, fs=fs)
    yf *= h
    return yf
def harmonic_notch_filter(yf, f0, fs, Q, harmonic_multiple):
    for i in range(1,harmonic_multiple+1):
        yf = notch_filter(yf, f0*i, fs, Q)
    return yf
	

	
class Data_Channel:
	
	#Initializer
    def __init__(self, name, device, time, freqs_fft, freqs_psd, wave, fft, psd, sampling_rate):
        self.name = name
        self.device = device
        self.x_ts = time
        self.x_fft = freqs_fft
        self.x_psd = freqs_psd
        self.y_ts = wave
        self.y_fft = fft 
        self.y_psd = psd
        self.fs = sampling_rate 
		
		
		
	#Instance Method
    def axis_limits(self, calculation):
        if calculation == "ts":
            Uylim, Lylim = max(getattr(self, "y_ts")) + max(getattr(self, "y_ts"))*0.1, min(getattr(self, "y_ts")) - min(getattr(self, "y_ts"))*0.1
            Uxlim_index = np.where(self.y_ts>=0.1*max(self.y_ts))[0][-1] 
            Uxlim = self.x_ts[Uxlim_index]
            Uxlim += ((Uxlim - min(self.x_ts))*0.1)
            Uxlim = int(Uxlim)
        elif calculation == "fft":
            Uylim, Lylim = max(getattr(self, "y_fft")[1::]) + max(getattr(self, "y_fft")[1::])*0.1, min(getattr(self, "y_fft")) - min(getattr(self, "y_fft"))*0.1
            Uxlim_index = np.where(self.y_fft[1::]>=0.1*max(self.y_fft[1::]))[0][-1] 
            Uxlim = self.x_fft[Uxlim_index+1]
            Uxlim += ((Uxlim - min(self.x_fft))*0.1)
            Uxlim = int(Uxlim)
        elif calculation ==  "psd":
            Uylim, Lylim = max(getattr(self, "y_psd")[1::]) + max(getattr(self, "y_psd")[1::])*0.1, min(getattr(self, "y_psd")[1::]) - min(getattr(self, "y_psd")[1::])*0.1
            Uxlim_index = np.where(self.y_psd[1::]>=0.1*max(self.y_psd[1::]))[0][-1]
            Uxlim = self.x_psd[Uxlim_index+1]
            Uxlim += ((Uxlim - min(self.x_psd))*0.1)
            Uxlim = int(Uxlim)
        Lxlim = min(self.x_ts)
        lims = np.array([Uylim, Lylim, Uxlim, Lxlim])
        return lims




# In[ ]:
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


def findOccurrences(s, ch):
	return [i for i, letter in enumerate(s) if letter == ch]

	
def fft_avg(fulldata, sampling_rates, windows_to_average, window_size, cutoff):
	assert (type(window_size+windows_to_average) == int), "Both window_size and windows_to_average must be integers"

	sums = 0
	samps = {}
	
	for i in range(fulldata["nInputs"]):
		for j in range(windows_to_average):
			samps[j] = fulldata["wave"+str(i)][j*window_size : ((j*window_size)+(window_size-1))]
			if type(cutoff) != type(None):
				samps[j] = butter_lowpass_filter(samps[j], cutoff, fs=sampling_rates[0], order=4)
			np.multiply(samps[j], gauss_wind(window_size, 0.5), out = samps[j], casting = 'unsafe')
			samps[j] = np.fft.fft(samps[j])
		
		for k in range(windows_to_average):
			sums += samps[k]
		
		FFT_output = sums/windows_to_average
		yF = 2./window_size*np.absolute(FFT_output[1:window_size//2])
		
		fulldata["FFT_wave"+str(i)] = yF
	fulldata["M"], fulldata["N"] = windows_to_average, window_size
	return fulldata


def fft_no_avg(fulldata, sampling_rates, cutoff):
	for i in range(fulldata["nInputs"]):
		if type(cutoff) != type(None):
			fulldata["FFT_wave"+str(i)] = butter_lowpass_filter(fulldata["wave"+str(i)], cutoff, fs=sampling_rates[0], order=4)
			np.multiply(fulldata["FFT_wave"+str(i)], gauss_wind(len(fulldata["times"])+1, 0.5), out = fulldata["FFT_wave"+str(i)], casting = 'unsafe')
		else:
			fulldata["FFT_wave"+str(i)] = np.empty(len(fulldata["wave"+str(i)]))
			np.multiply(fulldata["wave"+str(i)], gauss_wind(len(fulldata["times"])+1, 0.5), out = fulldata["FFT_wave"+str(i)], casting = 'unsafe')
		fulldata["FFT_wave"+str(i)] = np.fft.fft(fulldata["FFT_wave"+str(i)])
		fulldata["FFT_wave"+str(i)] = 2./len(fulldata["times"])*np.absolute(fulldata["FFT_wave"+str(i)][1:len(fulldata["times"])//2])
	return fulldata


def gauss_wind(N, sigma):
	M = N - 2 
	w = np.zeros(M+1)
	for i in range(M+1):
		w[i] = np.exp(-0.5*(((i - M)/2)/(sigma*M/2))*(((i - M)/2)/(sigma*M/2)))
	return w

	
def get_files(dataPath, key):
	print("The data path is: " + dataPath)
	pathdir = os.listdir(dataPath)
	file_array = []
	for file in pathdir:
		if (file.endswith(".root") and key in file):
			file_array.append(dataPath + file)
		   
	shortFilename = write_shortFilename(file_array[0])
	return file_array, shortFilename

	
def get_fft(fulldata, sampling_rates, windows_to_average, window_size, cutoff):
	if windows_to_average == None and window_size == None:
		fulldata = fft_no_avg(fulldata, sampling_rates, cutoff)
	
	else:
		fulldata = fft_avg(fulldata, sampling_rates, windows_to_average, window_size)
		
	for i in range(fulldata["nInputs"]):	
		try:	
			fulldata["FFT_f_wave"+str(i)] = np.linspace(0, 1.0/(2.0/sampling_rates[0]), window_size//2-1)
		except TypeError:
			fulldata["FFT_f_wave"+str(i)] = np.linspace(0, 1.0/(2.0/sampling_rates[0]), len(fulldata["times"])//2-1)
		
	return fulldata

	
def get_fullwave(filename, dataPath=None, num_of_files=None):
	filelist, shortFilename = get_files(dataPath, filename) 

	if num_of_files == None:
		print("Warning: Attempting to Merge all "+str(filename)+" files")
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
	
	for j in range(fulldata["nInputs"]):
		plot_title, txt, my_file = plot_waves(fulldata, sampling_rates, j, path, shortFilename, calculation)
		
	plt.title(plot_title)
	plt.legend(loc=0)
	plt.figtext(0.5, -0.034, txt, wrap=True, horizontalalignment='center', fontsize=12)
	plt.grid()
	plt.savefig(my_file)
	plt.show()
	
	#plot_mic(fulldata, sampling_rates, path, shortFilename, calculation)
		

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
	nInputs = len(Branch_List) - 4
	nSamples = getattr(tree, Branch_List[0])
	dt = getattr(tree, Branch_List[3])

	nEntries  = tree.GetEntriesFast()
	nTotal = nSamples*nEntries
	sampling_rate = int(1/dt)

	data = {}
	data["times"] = np.linspace(0, (nTotal-1)*dt, nTotal)

	for i,j in enumerate(Branch_List[4::]):
		data["wave"+str(i)] = (np.asarray([[getattr(tree, j)[0::]] for event in tree])).flatten()
		
	return data, sampling_rate, dt, nInputs
	
	
def plot_mic(fulldata, sampling_rates, path, shortFilename, calculation):
	plt.figure()
	
	if calculation == "ts":
		plt.plot(fulldata["times"], fulldata["wave"+str(fulldata["nInputs"]-1)], label="Mic")
		plt.title(shortFilename+"_TimeStream")
		txt = ""
		my_file2 = str(path) + str(shortFilename) + "_micTS.png"
	elif calculation == 'psd':
		plt.loglog(fulldata["PSD_f_wave"+str(fulldata["nInputs"]-1)], fulldata["PSD_wave"+str(fulldata["nInputs"]-1)], label="Mic")
		plt.title(shortFilename+"_PSD")
		try:
			txt = 'Average of '+str(fulldata["M"])+', '+str(fulldata["N"])+' point samples, Sampling Rate: '+str(sampling_rates[0]/1000)+' kHz'
		except KeyError:
			txt = ""
		my_file2 = str(path) + str(shortFilename) + "_micPSD.png"
	elif calculation == 'fft':
		plt.loglog(fulldata["FFT_f_wave"+str(fulldata["nInputs"]-1)], fulldata["FFT_wave"+str(fulldata["nInputs"]-1)], label="Mic")
		plt.title(shortFilename+"_FFT")
		try:
			txt = 'Average of '+str(fulldata["M"])+', '+str(fulldata["N"])+' point samples, Sampling Rate: '+str(sampling_rates[0]/1000)+' kHz'
		except KeyError:
			txt = ""
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
		plt.loglog(fulldata["FFT_f_wave"+str(j)], fulldata["FFT_wave"+str(j)], label = "wave "+str(j))
		plot_title = shortFilename+"_FFT"
		my_file = str(path) + str(shortFilename) + "_FFT.png"
		try:
			txt = 'Average of '+str(fulldata["M"])+', '+str(fulldata["N"])+' point samples, Sampling Rate: '+str(sampling_rates[0]/1000)+' kHz'
		except KeyError:
			txt = ""
	return plot_title, txt, my_file


def psd_avg(fulldata, sampling_rates, windows_to_average, window_size):
	
	assert (type(window_size+windows_to_average) == int), "Both window_size and windows_to_average must be integers"
		
	sums = 0
	samps = {}
	
	for i in range(fulldata["nInputs"]):
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
	for i in range(fulldata["nInputs"]):
		fulldata["PSD_f_wave"+str(i)], fulldata["PSD_wave"+str(i)] = sig.periodogram(fulldata["wave"+str(i)], fs=sampling_rates[0])
	return fulldata	

	
def stitch_data(filelist):
    #no longer necessary to do this since grab_data now returns the number of input channels
	#dummy_data, dummy_sampling_rate, dummy_dt = grab_data(filelist[0])
	#nInputs = len(dummy_data)-1
	
	fulldata = {}
	sampling_rates = np.array([])
	
	for i in range(len(filelist)):
		data, sampling_rate, dt, nInputs = grab_data(filelist[i])
		if i == 0:
			fulldata["times"] = np.array([])
			for k in range(nInputs):
				fulldata["wave"+str(k)] = np.array([])
		for j in range(nInputs):
			fulldata["wave"+str(j)] = np.append(fulldata["wave"+str(j)], data["wave"+str(j)])
		fulldata["times"] = np.append(fulldata["times"], (data["times"] +  ( len(sampling_rates) * (data["times"][-1]+dt) )   ))
		sampling_rates = np.append(sampling_rates, sampling_rate)
		if len(sampling_rates) > 1:
			if sampling_rates[-1] != sampling_rates[-2]:
				print("Warning: Files of Varying Sampling Rates are being Merged")
				
				
	fulldata["nInputs"] = nInputs
	return fulldata, sampling_rates
	

def write_shortFilename(filename):
	shortFilename = filename.split("/")[-1]
	shortFilename_index = findOccurrences(shortFilename, '_')[-2]
	shortFilename = shortFilename[0:shortFilename_index]
	return shortFilename



def main(filename, num_of_files=None, loglog=True, xlim=None, ylim=None, windows_to_average=None, window_size=None, cutoff=None):
    repoPath = os.environ['ANALYSISREPO']
    #repoPath = "/global/u2/z/zendejas/CUORE/Raw_Data/root_files/Analysis"
    dataPath = repoPath + "/data/root/"
    imgPath = repoPath + "/images/"
	
	#Request Channel Info from User		
    spaces, channel_list = confirm_channel_list()
    channel_dict = write_channel_list(spaces, channel_list)
    
	#Collect Object Values
    print("Collecting full datastream...")	
    fulldata, sampling_rates, shortFilename = get_fullwave(filename, dataPath, num_of_files)
    
    print("Getting PSD...")
    fulldata = get_psd(fulldata, sampling_rates, windows_to_average, window_size)
    
    print("Getting FFT...")
    fulldata = get_fft(fulldata, sampling_rates, windows_to_average, window_size, cutoff)
	
	#Create obj instances from fulldata
    channel_objs = []
    for i in range(fulldata["nInputs"]):
	    channel_objs.append("channel"+str(i))
    print(channel_objs)
    for i, j in enumerate(channel_objs):
        print(fulldata["wave"+str(i)])
        j = Data_Channel(channel_dict[i][0], channel_dict[i][1], fulldata["times"], fulldata["FFT_f_wave"+str(i)], fulldata["PSD_f_wave"+str(i)], fulldata["wave"+str(i)], fulldata["FFT_wave"+str(i)], fulldata["PSD_wave"+str(i)], sampling_rates[0])
        channel_objs[i] = j
	
    print(channel_objs)
    plotss = ["ts", "fft", "psd"]
    for i in plotss:
        print(i)
        hyper_plots(str(i), channel_objs, shortFilename, imgPath)
	
	#Hyper_Plts should handle plotting better
    #get_plot(fulldata, sampling_rates, imgPath, shortFilename, calculation="ts")
    #get_plot(PSDdata, sampling_rates, imgPath, shortFilename, calculation="psd")
	#get_plot(FFTdata, sampling_rates, imgPath, shortFilename, calculation="fft")
	
    plt.show()

#Execute Main

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


	