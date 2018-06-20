from scipy import signal
import csv
import numpy as np
import matplotlib.pyplot as plt
import pywt
from hybrid_tdmg_FE_new import hybrid_tdmg_FE_new

##############################################################
##############################################################

def unique(list1):
 
    unique_list = []
    for x in list1:
        # check if exists in unique_list or not
        if x not in unique_list:
            unique_list.append(x)
    # print list
    return unique_list


def filtering(x):
	f = 122.0/2.0
	b,a = signal.butter(4,0.05/f,btype='high')
	x_filtrato = signal.filtfilt(b,a,x)
	b_,a_ = signal.butter(4,40.0/f,btype='low')
	x_filtrato2 = signal.filtfilt(b,a,x_filtrato)
	return x_filtrato2


def indices(a, func):
    return [i for (i, val) in enumerate(a) if func(val)]

########################################################
#######################################################


def Boundary_Detection(ECG,name):

	xnew1 = filtering(np.asarray(ECG))
	#xnew2 = xnew1
	#xnew3 = xnew1 + xnew1

	Fs = 5000 
	xnew4 = signal.resample(xnew1,Fs)

	fs = 1000;
	y2_1 = signal.medfilt(xnew4,fs/5 + 1)
	y2_2 = signal.medfilt(y2_1,1 + (3*fs)/5)
	y2 = xnew4-y2_2;

	ecg_signal = y2

   	################ Wavelet Transforms #############
	cA_l1, cD_l1 = pywt.dwt(ecg_signal,'haar');
	cA_l2, cD_l2 = pywt.dwt(cA_l1,'haar');
	cA_l3, cD_l3 = pywt.dwt(cA_l2,'haar');
	cA_l4, cD_l4 = pywt.dwt(cA_l3,'haar');
	cA_l5, cD_l5 = pywt.dwt(cA_l4,'haar');

	#######  Finding ThreshHold #####################
	m = 1 - 1;
	n = 175 - 1;

	max_chunk = [];
	for i in range(0,len(cD_l3)/175):
		max_chunk.append(np.amax(max(cD_l3[m:n])))
		m = n + 1
		n = n + 175
	min_chunk = np.amin(max_chunk)
	T = min_chunk*0.45;


	cD_l3_logic = np.zeros(len(cD_l3))

	###### Storing logic 0 and 1 in the memory_______________________________________________________Start_____%%%%%
	for i in range(0,len(cD_l3)):
		if(abs(cD_l3[i]) >= T):
			cD_l3_logic[i] = 1
		else:
			cD_l3_logic[i] = 0

	cD_l3_logic[i] = cD_l3_logic[i]
	 
	for i in range(0,len(cD_l3)):
		if(i < len(cD_l3)-2):
			if(cD_l3_logic[i] == 1):
				if(cD_l3_logic[i] == 1 and cD_l3_logic[i+1] == 1):
					cD_l3_logic[i] = 0           
				else:
					cD_l3_logic[i] = cD_l3_logic[i]
			else:
				cD_l3_logic[i] = cD_l3_logic[i]

	j = 0
	k = 0

	max_index = []
	min_index = []
	value = []
	R_peak_value = []
	R_peak_index = []

	###############     Finding R_peak_index

	for i in range(0,len(cD_l3)):
		if(cD_l3_logic[i] == 1 and (i < (len(cD_l3) - 10)) and i > 10):		
			max_index.append(i)
			value.append(np.amin(cD_l3[max_index[j] - 10 : max_index[j] + 10]))
			min_index.append(np.argmin(cD_l3[max_index[j] - 10 : max_index[j] + 10]))
			if(min_index[k] > 10):
				min_index[k] = max_index[j] + (min_index[k] - 10)
			else:
				min_index[k] = max_index[j] - 11 + min_index[k]

			if(max_index[j] < min_index[k]):
				R_peak_index.append(np.argmin(ecg_signal[max_index[j]*8:min_index[k]*8]))
				R_peak_value.append(np.amin(ecg_signal[max_index[j]*8:min_index[k]*8]))
				R_peak_index[j] = R_peak_index[j] + max_index[j]*8
			else:
				R_peak_index.append(np.argmax(ecg_signal[min_index[j]*8:max_index[k]*8]))
				R_peak_value.append(np.amax(ecg_signal[min_index[j]*8:max_index[k]*8]))
				R_peak_index[j] = R_peak_index[j] + min_index[j]*8
			j = j + 1
			k = k + 1



	R_peak_index = np.unique(R_peak_index)


	############################################

	#########    Eliminating peaks with relative distance < 290 

	i = 0
	while (i < len(R_peak_index) - 1):
		if ((R_peak_index[i + 1] - R_peak_index[i]) <= 290):
			if(abs(ecg_signal[R_peak_index[i]]) > (abs(ecg_signal[R_peak_index[i + 1]]))):
				R_peak_index[i + 1] = 0
				i = i + 1
			elif(abs(ecg_signal[R_peak_index[i]]) < (abs(ecg_signal[R_peak_index[i + 1]]))):
				R_peak_index[i] = 0
				i = i + 1
		else:
			i = i + 1

	R_peak_index_new = []
	for i in R_peak_index:
		if(i != 0):
			R_peak_index_new.append(i)

	R_peak_index = R_peak_index_new

	############################################################


	#######################

	for i in range(0,len(R_peak_index)):
		if ((R_peak_index[i] > 56) and (R_peak_index[i] < 5000 -56)):
			index = np.argmax(np.absolute(ecg_signal[R_peak_index[i] - 56 : R_peak_index[i] + 56]))
			if(index >= 56):
				index = index + R_peak_index[i] - 56
			else:
				index = R_peak_index[i] - (56 - index)
			R_peak_index[i] = index
			R_peak_value[i] = ecg_signal[index]



	while (i < len(R_peak_index) - 1):
		if ((R_peak_index[i + 1] - R_peak_index[i]) <= 290):
			if(abs(ecg_signal[R_peak_index[i]]) > (abs(ecg_signal[R_peak_index[i + 1]]))):
				R_peak_index[i + 1] = 0
				i = i + 1
			elif(abs(ecg_signal[R_peak_index[i]]) < (abs(ecg_signal[R_peak_index[i + 1]]))):
				R_peak_index[i] = 0
				i = i + 1
		else:
			i = i + 1

	R_peak_index_new = []
	for i in R_peak_index:
		if(i != 0):
			R_peak_index_new.append(i)

	R_peak_index = R_peak_index_new
	R_peak_index = np.sort(R_peak_index)

	#################

	max_R_peak = max(ecg_signal[R_peak_index])

	################## Finding Missing Peaks within 1200 radius and 0.37 of max R Height

	index9 = [0]
	s = 0
	for i in range(0,len(R_peak_index) - 1):
		if(R_peak_index[i + 1] - R_peak_index[i] > 1200):
			for k in range(R_peak_index[i] + 50, R_peak_index[i + 1] - 50):
				if(ecg_signal[k] >= (0.37*max_R_peak)):
					index9[0] = (R_peak_index[i])
					index9.append(k)
					s = 1
					break;

	if(s):
		index8 = index9[1:]
		j = np.argmax(ecg_signal[index8])
		if(j - index9[0] > 400):
			R_peak_index.append(j)
			R_peak_index = np.sort(R_peak_index)

	index9 = [0]
	s = 0
	for i in range(0,len(R_peak_index) - 1):
		if(R_peak_index[i + 1] - R_peak_index[i] > 1200):
			for k in range(R_peak_index[i] + 50, R_peak_index[i + 1] - 50):
				if(ecg_signal[k] >= (0.37*max_R_peak)):
					index9[0] = (R_peak_index[i])
					index9.append(k)
					s = 1
					break;

	if(s):
		index8 = index9[1:]
		j = np.argmax(ecg_signal[index8])
		if(j - index9[0] > 400):
			R_peak_index.append(j)
			R_peak_index = np.sort(R_peak_index)


	###################

	max_R_peak = max(ecg_signal[R_peak_index])



	################### Finding R Peaks again with 0.5 times Max R height

	index9 = [0]
	s = 0
	for i in range(0,len(R_peak_index) - 2):
		if(R_peak_index[i + 2] - R_peak_index[i + 1] >=  1.8*(R_peak_index[i + 1] - R_peak_index[i] )):
			for k in range(R_peak_index[i + 1] + 50, R_peak_index[i + 2] - 50):
				if(ecg_signal[k] >= (0.5*max_R_peak)):
					index9[0] = (R_peak_index[i + 1])
					index9.append(k)
					s = 1
					break;

	if(s):
		index8 = index9[1:]
		j = np.argmax(ecg_signal[index8])
		if(j - index9[0] > 100):
			R_peak_index.append(j)
			R_peak_index = np.sort(R_peak_index)


	i = 0
	while(i < (len(R_peak_index)) - 1):
		if(R_peak_index[i + 1] - R_peak_index[i] < 250):
			R_peak_index[i + 1] = 0
			if (i < len(R_peak_index) - 2):
				i = i + 2
			else:
				i = i + 1
		elif(i > len(R_peak_index) - 2):
			i = i + 1
		else:
			i = i + 1

	R_peak_index_new = []
	for i in R_peak_index:
		if(i != 0):
			R_peak_index_new.append(i)

	R_peak_index = R_peak_index_new
	R_peak_index = np.sort(R_peak_index)


	################### More Filtering ??

	for i in range(0,len(R_peak_index)):
		if ((R_peak_index[i] > 56) and (R_peak_index[i] < 5000 -56)):
			index = np.argmax(np.absolute(ecg_signal[R_peak_index[i] - 56 : R_peak_index[i] + 56]))
			if(index >= 56):
				index = index + R_peak_index[i] - 56
			else:
				index = R_peak_index[i] - (56 - index)
			R_peak_index[i] = index
			R_peak_value[i] = ecg_signal[index]

	for i in range(0,len(R_peak_index)):
		if ((R_peak_index[i] > 56) and (R_peak_index[i] < 5000 -56)):
			index = np.argmax(np.absolute(ecg_signal[R_peak_index[i] - 56 : R_peak_index[i] + 56]))
			if(index >= 56):
				index = index + R_peak_index[i] - 56
			else:
				index = R_peak_index[i] - (56 - index)
			R_peak_index[i] = index
			R_peak_value[i] = ecg_signal[index]


	##################

	i = 0

	while (i < len(R_peak_index) - 1):
		if(abs(ecg_signal[R_peak_index[i + 1]]) < 0.2*abs(ecg_signal[R_peak_index[i]])):
			R_peak_index[i + 1] = 0
			i = i + 2
		i = i + 1

	R_peak_index_new = []
	for i in R_peak_index:
		if(i != 0):
			R_peak_index_new.append(i)

	R_peak_index = R_peak_index_new
	R_peak_index = np.sort(R_peak_index)

	#################

	index12 = []

	for i in range(0,len(R_peak_index) - 2):
		if (abs(ecg_signal[R_peak_index[i + 1]]) < (0.67*abs(ecg_signal[R_peak_index[i]]))):
			index = np.argmax(np.absolute(ecg_signal[R_peak_index[i + 1] - 10 : R_peak_index[i + 1] + 150]))
			if(index > 10):
				index = R_peak_index[i + 1] + index - 10
			else:
				index = R_peak_index[i + 1] - index
			index12.append(index)
			R_1 = R_peak_index[0:i + 1]
			R_2 = R_peak_index[i + 3:]
			R_peak_index = R_1
			R_peak_index.append(index9)		
			R_peak_index = R_peak_index + R_2

	R_peak_index = np.sort(R_peak_index)


	################### More Filtering ??

	for i in range(0,len(R_peak_index)):
		if ((R_peak_index[i] > 56) and (R_peak_index[i] < 5000 -56)):
			index = np.argmax(np.absolute(ecg_signal[R_peak_index[i] - 56 : R_peak_index[i] + 56]))
			if(index >= 56):
				index = index + R_peak_index[i] - 56
			else:
				index = R_peak_index[i] - (56 - index)
			R_peak_index[i] = index
			R_peak_value[i] = ecg_signal[index]

	for i in range(0,len(R_peak_index)):
		if ((R_peak_index[i] > 56) and (R_peak_index[i] < 5000 -56)):
			index = np.argmax(np.absolute(ecg_signal[R_peak_index[i] - 56 : R_peak_index[i] + 56]))
			if(index >= 56):
				index = index + R_peak_index[i] - 56
			else:
				index = R_peak_index[i] - (56 - index)
			R_peak_index[i] = index
			R_peak_value[i] = ecg_signal[index]

	####################

	################### Something to do with negative R-index values

	R_peak_prev = R_peak_index
    
   	for i in range(0,len(R_peak_index)):
		if((ecg_signal[R_peak_index[i]] < 0) and (i == 0)):
			if(R_peak_index[0] > 89):
				index = np.argmax((ecg_signal[R_peak_index[i] - 89 : R_peak_index[i]]))
				R_peak_index[i] = R_peak_index[i] - 89 + index
			else:
				index = np.argmax((ecg_signal[1 : R_peak_index[i]]))
				R_peak_index[i] = index
		elif(ecg_signal[R_peak_index[i]] < 0):
			index = np.argmax((ecg_signal[R_peak_index[i] - 89 : R_peak_index[i]]))
			R_peak_index[i] = R_peak_index[i] - 89 + index

	R_peak_value = ecg_signal[R_peak_index]


	################## Boundary Calculations

	boundary = []

	already_in = 0;

	for i in range(0,len(R_peak_index) - 1):
		if((R_peak_index[i + 1] - R_peak_index[i] < 658) and ((((R_peak_index[i + 1] - R_peak_index[i])/2) + R_peak_index[i] + 190) < R_peak_index[i + 1])):
			array = np.arange(int(round((R_peak_index[i] + R_peak_index[i + 1])/2)),int(round((R_peak_index[i] + R_peak_index[i + 1])/2) + 191)).tolist()
			for k in array:
				if(((ecg_signal[k]*ecg_signal[k + 1]) < 0) or (ecg_signal[k] == 0)):
					if (already_in == 0):
						boundary.append(k)
						already_in = 1
					else:
						boundary[i] = k
					already_in = 0
					break
				else:
					if(already_in == 0):
						boundary.append((R_peak_index[i] + R_peak_index[i + 1])/2)
						already_in = 1
					else:
						boundary[i] = (R_peak_index[i] + R_peak_index[i + 1])/2
		elif(R_peak_index[i + 1] - R_peak_index[i] < 190):
			boundary.append((R_peak_index[i] + R_peak_index[i + 1])/2) 
		else:
			boundary.append(((R_peak_index[i] + R_peak_index[i + 1])/2) + 40)
		already_in = 0

	################## R peak correction using S pattern Morphology

	if (R_peak_index[0] < boundary[0]) :
		boundary_dummy = boundary
		R_peak_index_dummy = R_peak_index[1:]
		R_peak_prev_dummy = R_peak_prev[1:]
	else:
		boundary_dummy = boundary
		R_peak_index_dummy = R_peak_index
		R_peak_prev_dummy = R_peak_prev

	grad = np.diff(ecg_signal)
	
	for i in range(0,len(R_peak_index_dummy)):
		if(R_peak_index_dummy[i] - boundary_dummy[i] < 142) and (ecg_signal[R_peak_index_dummy[i]] < 0.19):
			pos = np.argmax(grad[R_peak_index_dummy[i] : R_peak_prev_dummy[i] - 14])
			pos = pos + R_peak_index_dummy[i] - 10
			pos1 = np.argmax(ecg_signal[pos: pos + 18])
			pos1 = pos1 + pos
			R_peak_index_dummy[i] = int(pos1)


	if (len(R_peak_index) > len(R_peak_index_dummy)):
		R_peak_index[1:] = R_peak_index_dummy
	else:
		R_peak_index = R_peak_index_dummy        

	################## HRV Calculations

	HRV = []

	for i in xrange(0,len(R_peak_index) - 1):
		HRV.append(R_peak_index[i + 1] - R_peak_index[i])

	################## More Boundary Calculations using HRV

	for i in range(0,len(HRV) - 1):
		if((abs(HRV[i+1] - HRV[i]) >= 100) and ((((R_peak_index[i + 1] - R_peak_index[i])/2) + R_peak_index[i] + 190) < R_peak_index[i + 1])):
			array = np.arange(int(round((R_peak_index[i] + R_peak_index[i + 1])/2)),int(round((R_peak_index[i] + R_peak_index[i + 1])/2) + 191)).tolist()
			for k in array:
				if(((ecg_signal[k]*ecg_signal[k + 1]) < 0) or (ecg_signal[k] == 0)):
					boundary[i] = (k)
					break
		else:
			boundary[i] = boundary[i]


	################## Boundary End and Start

	if((R_peak_index[0] >= (R_peak_index[1] - boundary[0])) and (R_peak_index[1] - R_peak_index[0] > 500)):
		boundary_start = (R_peak_index[1] - boundary[0])
		boundary_start = (R_peak_index[0] - boundary_start)
	else:
		boundary_start = boundary[0]

	if((len(ecg_signal) - R_peak_index[-1]) >= (boundary[-1] - R_peak_index[-2])):
		boundary_end = boundary[-1] - R_peak_index[-2]
		boundary_end = boundary_end + R_peak_index[-1]
	else:
		boundary_end = boundary[-1]

	boundary.append(boundary_start)
	boundary.append(boundary_end)
	boundary = np.sort(boundary)
	boundary = unique(boundary)

	##################	Setting Proper Boundaries if value > 0.03

	for k in range(0,len(boundary)):
		if ((ecg_signal[boundary[k]] > 0.03) and k < (len(boundary))):
			for kk in range(boundary[k],boundary[k] + 74):
				if((kk < len(ecg_signal)) and ((ecg_signal[kk]*ecg_signal[kk + 1] < 0) or (ecg_signal[kk] < 0.015))):
					boundary[k] = kk
		elif ((k == len(boundary)) and (len(ecg_signal) - boundary[-1] <= 100) and (abs(ecg_signal[boundary[k]]) > 0.03)):
			for kk in range(boundary[k],len(ecg_signal)):
				if((kk < len(ecg_signal)) and ((ecg_signal[kk]*ecg_signal[kk + 1] < 0) or (ecg_signal[kk] < 0.02) or (ecg_signal[kk + 1] < 0.02))):
					boundary[k] = kk
		elif ((k == len(boundary)) and (len(ecg_signal) - boundary[-1] > 100) and (abs(ecg_signal[boundary[k]]) > 0.03) and (ecg_signal[boundary[-1]] + 1 > ecg_signal[boundary[-1]])):
			for kk in range((R_peak_index[-1] + boundary[-1])/2,boundary[-1] + 80):
				if((kk < len(ecg_signal)) and ((ecg_signal[kk]*ecg_signal[kk + 1] < 0) or (ecg_signal[kk] < 0.015))):
					boundary[k] = kk
		

	if(ecg_signal[boundary[-1]] > 1):
		if((boundary[-1] < len(ecg_signal)) and (ecg_signal[boundary[-1] + 1] > ecg_signal(boundary[-1])) and (boundary[-1] - R_peak_index[-1] > 400)):
			for k in range(boundary[-1] - 125,boundary[-1] + 1):
				if((kk < len(ecg_signal)) and ((ecg_signal[kk]*ecg_signal[kk + 1] < 0) or (ecg_signal[kk] < 0.015))):
					boundary[-1] = kk

	

	##################


	max_R_peak = max(ecg_signal[R_peak_index])

	boundary_line = ecg_signal[boundary]

	R_peak_value = ecg_signal[R_peak_index]

	plt.plot(boundary,boundary_line,'ro')
	plt.plot(R_peak_index,R_peak_value,'bo')
	plt.plot(ecg_signal,'b')
	plt.savefig(File_name + ".png")
	plt.show()

	for j in range(0,len(boundary) - 1):
		hybrid_tdmg_FE_new(ecg_signal[boundary[j] : boundary[j + 1]],range(boundary[j],boundary[j + 1]), R_peak_index[j] - boundary[j])



		

File_name = "mitdb_113_22_27"

with open(File_name + ".csv") as f:
    readCSV = csv.reader(f, delimiter=" ")
    ECG_str = []
    for row in readCSV:
    	data = row[0]
    	ECG_str.append(data)

ECG_signal = []

for i in range(0,len(ECG_str)):
	ECG_signal.append(float(ECG_str[i]))


Boundary_Detection(ECG_signal,File_name)

