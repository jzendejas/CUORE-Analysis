#import ROOT stuff
from ROOT import TF1, iostream, fstream, vector, TStyle, TMath, TH1D, TFile, gInterpreter
from ROOT import std, TBranch, TTree, TChain
#import numpy stuff
import numpy as np
import matplotlib
#matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import scipy as sp
import scipy.signal as sig
import fnmatch
import os
import sys
#get_ipython().run_line_magic('matplotlib', 'inline')

###GENERAL GLOBAL FUNCTIONS
def assignmets(fulldata, channel_dict, sampling_rates, channel_objs, data):
    device, name, time, wave, sampling_rate, channel_number = data
    try:
        if device == "acc":
            channel_objs[device][name] = accelerometer(name, time, wave, sampling_rate, channel_number)
        elif device == "mic":
            channel_objs[device][name] = microphone(name, time, wave, sampling_rate, channel_number)
        elif device == "TES":
            channel_objs[device][name] = TES(name, time, wave, sampling_rate, channel_number)
        
    except KeyError:
        channel_objs[device] = {}
        
        if device == "acc":
            channel_objs[device][name] = accelerometer(name, time, wave, sampling_rate, channel_number)
        elif device == "mic":
            channel_objs[device][name] = microphone(name, time, wave, sampling_rate, channel_number)
        elif device == "TES":
            channel_objs[device][name] = TES(name, time, wave, sampling_rate, channel_number)
            
    return channel_objs

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

def create_objects(fulldata, channel_dict, sampling_rates):
    channel_objs = {}
    
    for i in range(fulldata["nInputs"]):
        device = channel_dict[i][1]
        name = channel_dict[i][0]
        time = fulldata["times"]
        wave = fulldata["wave"+str(i)]
        sampling_rate = sampling_rates[0]
        channel_number = i
        data = (device, name, time, wave, sampling_rate, channel_number)
        channel_objs = assignmets(fulldata, channel_dict, sampling_rates, channel_objs, data)
            
    return channel_objs

def custom_plotter(channel_objs, plot_objs, calc, partial_myfile):
    organize = []
    for l in channel_objs:
        for m in channel_objs[l]:
            organize.append(channel_objs[l][m])
    print(organize)
    for i in range(len(plot_objs)):
        plt.figure()
        for j in getattr(plot_objs[i], "csv"):
            current_channel = next((x for x in organize if x.channel_number == int(j)), None)
            print(current_channel)
            getattr(plt, getattr(plot_objs[i], "scaling"))(getattr(current_channel, "x_"+str(calc)), getattr(current_channel, "y_"+str(calc)), label=current_channel.name)
            peaks, x1, y1 = current_channel.find_peaks(calc)
            getattr(plt, getattr(plot_objs[i], "scaling"))(x1[peaks], y1[peaks], "x")
        plt.title(str(getattr(plot_objs[i], "title"))+"_"+str(calc))
        plt.legend()
        for k in ["x_lim", "y_lim"]:
            try_options(plot_objs[i], k)
		
        try:
            plt.savefig(str(partial_myfile)+str(getattr(plot_objs[i], "savefig"))+str(calc))
        except AttributeError:
            pass

        plt.show()
    return
	
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

def initialze_plot_object(channel_dict):
    num_of_plots = int(input("How many plots would you like?"))
    plot_objs = []
    print("These are the channel assignments you gave.")
    print(channel_dict)
    print("Please write all listed inputs as comma separated lists without spaces")
    print("eg: 1,2,3")
    for i in range(num_of_plots):
        print("List the channels you want on this plot")
        print("e.g.: 1,3,7")
        csv_list = list(input("Plot "+str(i+1)+":"))
        print("Enter the plot title")
        plot_title = input("Plot "+str(i+1)+" title:")
        scaling = input("Would you like linear or loglog scaling? (linear/loglog)")
        temp_plot = plot(i, csv_list, plot_title, scaling)
        temp_plot.sort("csv")
        save = input("Would you like to set a filename for img saving?(y/n)")
        if save == "y" or save == "Y":
            save_title = input("Enter title for saving")
            temp_plot.assign_save(save_title)
        limits = input("Would you like to set axis limits? (y/n)")
        if limits == "y" or save == "Y":
            axes = input("which axes limits would you like to assign? (x/y/both)")
            if axes == "x":
                x_lim = list(input("List x axis limits like so: 0,100"))
                temp_plot.assign_x_lim(x_lim)
            elif axes == "y":
                y_lim = list(input("List y axis limits like so: 0,100"))
                temp_plot.assign_y_lim(y_lim)
            elif axes == "both":
                x_lim = list(input("List x axis limits like so: 0,100"))
                y_lim = list(input("List y axis limits like so: 0,100"))
                temp_plot.assign_x_lim(x_lim)
                temp_plot.assign_y_lim(y_lim)
        plot_objs.append(temp_plot)
        
    return plot_objs

def plotter(channel_objs, calculation, partial_myfile):
    
    for i in channel_objs:
        plt.figure()
        for j in channel_objs[i]:
            getattr(channel_objs[i][j], "plot")(calculation)

        my_file = str(partial_myfile)+str(type(channel_objs[i][j]))[17:20]+"_"+str(calculation)
        plt.title(str(type(channel_objs[i][j]))[17:20]+"_"+str(calculation))
        plt.savefig(my_file)
        plt.legend()
        plt.show()
		
def plot_decision_map(channel_dict, channel_objs, partial_myfile):
    plotting_method = request()
    if plotting_method == "auto":
        for k in ["ts", "fft", "psd"]:
            plotter(channel_objs, k, partial_myfile)
    elif plotting_method == "custom":
        plot_objs = initialze_plot_object(channel_dict)
        for k in ["ts", "fft", "psd"]:
            custom_plotter(channel_objs, plot_objs, k, partial_myfile)

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

def request():
    plotting_method = input("Would you like 'auto' or 'custom' plotting?")
    return plotting_method

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
	
def try_options(obj, option):  
    try:
        getattr(plt, option)(getattr(obj, option))
    except AttributeError:
        pass
    return

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
	
def write_shortFilename(filename):
	shortFilename = filename.split("/")[-1]
	shortFilename_index = findOccurrences(shortFilename, '_')[-2]
	shortFilename = shortFilename[0:shortFilename_index]
	return shortFilename

###GENERAL CLASS METHODS

def fft_avg_in(wave, times, sampling_rate, windows_to_average, window_size, cutoff):
    assert (type(window_size+windows_to_average) == int), "Both window_size and windows_to_average must be integers"

    sums = 0
    samps = {}

    for j in range(windows_to_average):
        samps[j] = wave[j*window_size : ((j*window_size)+(window_size-1))]
        if type(cutoff) != type(None):
            samps[j] = butter_lowpass_filter(samps[j], cutoff, fs=sampling_rate, order=4)
        np.multiply(samps[j], sig.hann(len(times)), out = samps[j], casting = 'unsafe')
        samps[j] = np.fft.fft(samps[j])

    for k in range(windows_to_average):
        sums += samps[k]

    FFT_output = sums/windows_to_average
    yF = 2./window_size*np.absolute(FFT_output[1:window_size//2])

    return yF

def fft_bins(times, sampling_rate, windows_to_average, window_size):
    try:
        freqs = np.linspace(0, 1.0/(2.0/sampling_rate), window_size//2-1)
    except TypeError:
        freqs = np.linspace(0, 1.0/(2.0/sampling_rate), len(times)//2-1)

    return freqs

def fft_no_avg_in(wave, times, sampling_rate, cutoff):

    if type(cutoff) != type(None):
        FFT_wave = butter_lowpass_filter(wave, cutoff, fs=sampling_rate, order=4)
        np.multiply(FFT_wave, sig.hann(len(times)), out = FFT_wave, casting = 'unsafe')
    else:
        FFT_wave = np.empty(len(wave))
        np.multiply(wave, sig.hann(len(times)), out = FFT_wave, casting = 'unsafe')
    FFT_wave = np.fft.fft(FFT_wave)
    FFT_wave = 2./len(times)*np.absolute(FFT_wave[1:len(times)//2])
    return FFT_wave

def psd_avg(wave, sampling_rate, windows_to_average, window_size):

    assert (type(window_size+windows_to_average) == int), "Both window_size and windows_to_average must be integers"

    sums = 0
    samps = {}


    for j in range(windows_to_average):
        samps[j] = wave[j*window_size : ((j*window_size)+(window_size-1))]
        f, samps[j] = sig.welch(samps[j], fs=sampling_rate)

    for k in range(windows_to_average):
        sums += samps[k]

    PSD = sums/windows_to_average

    return f, PSD

def psd_no_avg(wave, sampling_rate):
    f, PSD = sig.welch(wave, fs=sampling_rate)
    return f, PSD

###CLASS DEFINITIONS

class Data_Channel:
	
	#Initializer
    def __init__(self):
        return
		
		
		
	#Instance Methods
    def fft(self, windows_to_average=None, window_size=None, cutoff=None):
        if windows_to_average == None:
            y_output = fft_no_avg_in(self.y_ts, self.x_ts, self.fs, cutoff)
        elif windows_to_average != None:
            y_output = fft_avg_in(self.y_ts, self.x_ts, self.fs, windows_to_average, window_size, cutoff)
        
        x_output = fft_bins(self.x_ts, self.fs, windows_to_average, window_size)
        
        self.x_fft, self.y_fft = x_output, y_output
        
        return x_output, y_output
    
    def psd(self, windows_to_average=None, window_size=None):
        if windows_to_average == None:
            f, PSD = psd_no_avg(self.y_ts, self.fs)
        elif windows_to_average != None:
            f, PSD = psd_avg(self.y_ts, self.fs, windows_to_average, window_size)
            
        self.x_psd, self.y_psd = f, PSD
        
        return f, PSD
        
    
    
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
    
class plotting:
    #Initializer
    def __init__(self):
        return
    
    #Instance Methods
    def find_peaks(self, calculation):
        ext_y = "y_"+str(calculation)
        ext_x = "x_"+str(calculation)
        y = getattr(self, ext_y)
        x = getattr(self, ext_x)
        peaks, _ = sig.find_peaks(y, threshold = np.std(y))
        n = 0
        while len(peaks) > 3:
            n += 1
            peaks, _ = sig.find_peaks(y, threshold = 3*np.std(y), prominence = n*np.std(y))
        return peaks, x, y
    
    def plot(self, calc):
        if calc == "ts":
            plt.plot( getattr(self, "x_"+str(calc)), getattr(self, "y_"+str(calc)), label=self.name )
        else:
            try:
                peaks, x, y = self.find_peaks(calculation=calc)
                plt.loglog( getattr(self, "x_"+str(calc)), getattr(self, "y_"+str(calc)), label=self.name )
                plt.loglog( x[peaks], y[peaks], "x")
            except AttributeError:
                print("Warning: "+str(calc)+" called without being computed.")
                print("Call: device."+str(calc)+"() prior to plotting")

class plot:
    def __init__(self, plot_number, csv_list, plot_title, scaling):
        self.num = plot_number
        self.csv = csv_list
        self.title = plot_title
        if scaling == "linear":
            self.scaling = "plot"
        else:
            self.scaling = scaling
		
        
    def sort(self, att): 
        out = findOccurrences(getattr(self, att), ",")
        for i in range(len(getattr(self, att))):
            try:
                getattr(self, att).remove(",")
            except ValueError:
                break
        return 

    def assign_x_lim(self, x_lim):
        self.x_lim = x_lim
        self.sort("x_lim")
        return
        
    def assign_y_lim(self, y_lim):
        self.y_lim = y_lim
        self.sort("y_lim")
        return
    
    def assign_save(self, save_title):
        self.savefig = save_title
        return

class accelerometer(Data_Channel, plotting):
    def __init__(self, name, time, wave, sampling_rate, channel_number):
        self.name = name
        self.x_ts = time
        self.y_ts = wave
        self.fs = sampling_rate
        self.channel_number = channel_number
        return
    
class microphone(Data_Channel, plotting):
    def __init__(self, name, time, wave, sampling_rate, channel_number):
        self.name = name
        self.x_ts = time
        self.y_ts = wave
        self.fs = sampling_rate
        self.channel_number = channel_number
        return

class TES(Data_Channel, plotting):
    def __init__(self, name, time, wave, sampling_rate, channel_number):
        self.name = name
        self.x_ts = time
        self.y_ts = wave
        self.fs = sampling_rate
        self.channel_number = channel_number
        return

###	MAIN FUNCTION

def main(filename, num_of_files, windows_to_average=None, window_size=None, cutoff=None):
    repoPath = os.environ['ANALYSISREPO']
    dataPath = repoPath + "/data/root/"
    imgPath = repoPath + "/images/"

    #Request Channel Info from User		
    spaces, channel_list = confirm_channel_list()
    channel_dict = write_channel_list(spaces, channel_list)
    
    #Collect Object Values
    print("Loading in Data")
    fulldata, sampling_rates, shortFilename = get_fullwave(filename, dataPath, num_of_files)
    print("Creating Objects from Data")
    channel_objs = create_objects(fulldata, channel_dict, sampling_rates)
    partial_myfile = str(imgPath)+str(shortFilename)
    for i in channel_objs:
        for j in channel_objs[i]:
            print("Getting "+str(getattr(channel_objs[i][j], "name"))+" fft")
            getattr(channel_objs[i][j], "fft")(windows_to_average, window_size, cutoff)
            print("Getting "+str(getattr(channel_objs[i][j], "name"))+" psd")
            getattr(channel_objs[i][j], "psd")(windows_to_average, window_size)

    plot_decision_map(channel_dict, channel_objs, partial_myfile)
	

    



### EXECUTE MAIN

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
	


