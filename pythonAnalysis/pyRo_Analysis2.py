def findOccurrences(s, ch):
	return [i for i, letter in enumerate(s) if letter == ch]
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
		channel_list = input("Please enter channel number with associated device name (ex: 0:accx 1:mic1 2:TES1  etc.)")
		spaces = findOccurrences(channel_list, " ")
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


spaces, channel_list = confirm_channel_list()
channel_dict = write_channel_list(spaces, channel_list)


for i in channel_dict:
	print(channel_dict[i][:-1])


#Create obj class
class Data_Channel:
	
	#Initializer
	def __init__(self, name, device, time, freqs_fft, freqs_psd, wave, fft, psd, sampling_rate):
		self.name = name
		self.device = device
		self.time = time
		self.freqsf = freqs_fft
		self.freqp = freqs_psd
		self.wave = wave
		self.fft = fft 
		self.psd = psd
		self.fs = sampling_rate
		
#Create obj instances from fulldata
channel_objs = []
for i in range(fulldata[nInputs]):
	channel_objs.append("channel"+str(i))
for i,j in enumerate(channel_objs):
	j = Data_Channel(channel_dict[i], channel_dict[i][0], channel_dict[i][1], fulldata["times"], fulldata["FFT_f_wave"+str(i)], fulldata["PSD_f_wave"+str(i)], i, fulldata["FFT_wave"+str(i)], fulldata["PSD_wave"+str(i)], sampling_rates[0])
	channel_objs.append(wave+str(i))
	

def hp_TS(channel_objs):
	plt.figure()
	plt.title("TimeStream")
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
				plt.show()
				plt.figure()
				plt.plot(getattr(j, time), getattr(j, wave), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
				
	
def hp_FFT(channel_objs):
	plt.figure()
	plt.title("FFT")
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
				plt.show()
				plt.figure()
				plt.plot(getattr(j, freqsf), getattr(j, fft), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
				
				
def hp_PSD(channel_objs):
	plt.figure()
	plt.title("PSD")
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
				plt.show()
				plt.figure()
				plt.plot(getattr(j, freqp), getattr(j, psd), label=str(getattr(j, name)))
				prev_dev = getattr(j, device)
				
				
def hyper_plots(calculation):
	if calculation == "ts":
		hp_TS(channel_objs)
		
	elif calculation == "fft":
		hp_FFT(channel_objs)
		
	elif calculation == "psd":
		hp_PSD(channel_objs)
	else:
		print("Input not Understood: Please use Keys; 'ts', 'fft', 'psd'")
		
	
	
	
	
	
	