import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.signal as sig
import fnmatch
import os
import csv
from pathlib import Path
import sys

def axis_limits(x, y_max, y_min, xlim, ylim, psd=False):
    
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

def calc_fs(time):
    step_size = time[5]-time[4]
    fs = int(1/step_size)
    print("Sampling Rate = "+str(fs)+" Hz")
    return fs

def calc_psd(wave0, wave1, wave2, sampling_rate, windows_to_average, window_size ):
    fx, px = psd_averaging(windows_to_average, window_size, sampling_rate, wave0)
    fy, py = psd_averaging(windows_to_average, window_size, sampling_rate, wave1)
    fz, pz = psd_averaging(windows_to_average, window_size, sampling_rate, wave2)
    
    return fx, fy, fz, px, py, pz

def calculate_FFT(windows_to_average, window_size, sampling_rate, cutoff, micG):
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
    yfG = 2/window_size*np.abs(FFT_output[1:window_size//2])
    
    return yfG

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

def gauss_wind(N, sigma):
    M = N - 2 
    w = np.zeros(M+1)
    for i in range(M+1):
        w[i] = np.exp(-0.5*(((i - M)/2)/(sigma*M/2))*(((i - M)/2)/(sigma*M/2)))
    return w

def get_FFT(filename, path, window_size=None, windows_to_average=None, cutoff=None, xlim=None, ylim=None ):
    
    my_path = path
    time, wave0, wave1, wave2, filename = get_file_data(filename)
    sampling_rate = calc_fs(time)
    
    if type(window_size) == type(None):
        window_size = len(time)
    
    if type(windows_to_average) == type(None):
        windows_to_average = 1
    
    xf = np.linspace(0, 1.0/(2.0/sampling_rate), window_size//2-1)
          
    yx = np.copy(wave0)
    yy = np.copy(wave1)
    yz = np.copy(wave2)

    yfx = calculate_FFT(windows_to_average, window_size, sampling_rate, cutoff, yx)
    yfy = calculate_FFT(windows_to_average, window_size, sampling_rate, cutoff, yy)
    yfz = calculate_FFT(windows_to_average, window_size, sampling_rate, cutoff, yz)
    
    plot_fft(xf, yfx, yfy, yfz, windows_to_average, window_size, filename, xlim, ylim, my_path)

def get_file_data(filename):
    try: 
        err = "Error: 0-D Array, skipping "
        count = 0
        time =  []
        wave0 = []
        wave1 = []
        wave2 = []


        with open(filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                time.append(float(row[0]))
                wave0.append(float(row[1]))
                wave1.append(float(row[2]))
                wave2.append(float(row[3]))

    except FileNotFoundError:
        #Error Handling: If a filenotfound error occurs during the data processing
        time = np.linspace(0,50,500000)
        wave = np.ones(500000)
        print("FILE NOT FOUND:")
        print("Plots petaining to: "+str(filename)+" should be ignored")
        filename = "FILE NOT FOUND: PLEASE IGNORE"
        
    return time, wave0, wave1, wave2, filename

def get_files(dataPath, key, sort = True):
    print("The data path is: " + dataPath)
    pathdir = os.listdir(dataPath)
    file_array = []
    for file in pathdir:
        if (file.endswith(".csv") and key in file):
           file_array.append(dataPath + file)
    return file_array

def get_psd(filename, imgPath, loglog=True, xlim=None, ylim=None, windows_to_average=None, window_size=None):
    #input file name as string (ex: "file")
    time, wave0, wave1, wave2, filename = get_file_data(filename)
    sampling_rate = calc_fs(time)
    
    if type(window_size) == type(None):
        window_size = len(time)
    
    if type(windows_to_average) == type(None):
        windows_to_average = 1
    
    fx, fy, fz, px, py, pz = calc_psd(wave0, wave1, wave2, sampling_rate, windows_to_average, window_size)
    plot_psd(fx, fy, fz, px, py, pz, filename, imgPath, loglog, xlim, ylim, windows_to_average, window_size)

def plot_psd(fx, fy, fz, px, py, pz, filename, imgPath, loglog, xlim, ylim, windows_to_average, window_size):
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
    plt.title("PSD of File: "+str(filename)) 
    plt.grid()

    f_max, p_max = find_max(fx,fy,fz,px,py,pz)
    p_min = find_min(px,py,pz)
    axis_limits(f_max, p_max, p_min, xlim, ylim, psd=True)
    
    txt = 'Average of '+str(windows_to_average)+', '+str(window_size)+' point samples'
    plt.figtext(0.5, -0.034, txt, wrap=True, horizontalalignment='center', fontsize=12)
    imgPath = filename.split(".")[0] + "_psd.png" 
    plt.savefig(imgPath)
    #plt.show()

def plot_fft(xf, yfx, yfy, yfz, M, N, filename, xlim, ylim, imgPath):
        
    plt.figure()
    plt.loglog(xf, yfx, label = 'X; Gaussian Window' )
    plt.loglog(xf, yfy, label = 'Y; Gaussian Window' )
    plt.loglog(xf, yfz, label = 'Z; Gaussian Window' )
    plt.legend()
    plt.title("FFT of File: "+str(filename) + "( M = "+str(M)+', N= '+str(N)+")")
    
    xf_max, yf_max = find_max(xf, xf, xf, yfx, yfy, yfz)
    yf_min = find_min(yfx, yfy, yfz)
    axis_limits(xf_max, yf_max, yf_min, xlim, ylim)
    
    txt = 'Average of '+str(M)+', '+str(N)+' point samples'
    plt.figtext(0.5, -0.034, txt, wrap=True, horizontalalignment='center', fontsize=12)
    print("filename is : " + filename) 
    imgPath = filename.split(".")[0] +"_fft.png"
    plt.savefig(imgPath)
    #plt.show()  

def plotter(filename, imgPath):
    #input file name as string (ex: "file")
    time, wave0, wave1, wave2, filename = get_file_data(filename)
    
    plt.plot(time, wave0, color = 'r', label = "x")
    plt.plot(time, wave1, color = 'b', label = "y")
    plt.plot(time, wave2, color = 'g', label = "z")
    plt.xlabel("Time [Seconds]")
    plt.ylabel("Amplitude [mV ?]")
    plt.legend()
    plt.grid()
    imgPath = imgPath + filename.split('.')[0].split('/')[-1] + "_timestream.png"
    print(imgPath + "Img path")
    title = input("Enter timestream title: ")
    plt.savefig(imgPath)
    #plt.show()

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

def main(key, Unix, sort, window_size=None, windows_to_average=None, cutoff=None, xlim=None, ylim=None, loglog=True):
    repoPath = os.environ['ANALYSISREPO']
    dataPath = repoPath + "/data/csv/"
    imgPath = repoPath + "/images/"

    if Unix == False:
        dataPath = PureWindowsPath(dataPath)
        imgPath = PureWindowsPath(imgPath)
    
    #print("Calling get_files on: " + str(dataPath) + " with key " +  str(key))
    file_array = get_files(dataPath, key, sort=sort)
    for i in file_array:
        plotter(i, imgPath)
        #get_psd(i, imgPath, loglog, xlim, ylim, windows_to_average, window_size)
        #get_FFT(i, imgPath, window_size, windows_to_average, cutoff, xlim, ylim)

    plt.show()

argc = len(sys.argv)
if argc == 2:
    key = str(sys.argv[1])
    main(key,Unix = True, sort = False)
elif (argc > 2 and argc <= 4):
    Unix = sys.argv[2]
    sort = sys.argv[3]	
    if (Unix == True or Unix == False and (sort == True or sort == False)):
        main(key, Unix, sort)
else: 
    print("Run as python Accel_Analysis_csv.py <key> <Unix = True or False> <sort = True or False> ")
    
    
