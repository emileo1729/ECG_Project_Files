# -*- coding: utf-8 -*-
"""
Created on Mon Jun 18 21:41:12 2018

@author: DELL
"""
import numpy as np
from notching import notching
from zerocross import zerocross

def morphologies(ecg_signal,start_qrs,end_qrs):
    
    y = []
    crossings = []
    
    for i in range(start_qrs + 10,end_qrs + 10):
        if((ecg_signal[i] > 0 and ecg_signal[i + 1] < 0) or (ecg_signal[i] < 0 and ecg_signal[i + 1] > 0)):
            crossings.append(i)
            
    crossings = [start_qrs] + crossings + [end_qrs]
    
    for k in range(0,len(crossings) - 1):
        window = ecg_signal[crossings[k] : crossings[k + 1]]
        m = 1
        grad = np.diff(window)
        y1 = []
        for m in range(0,len(grad) - 1):
            if((grad[m] > 0 and grad[m + 1] < 0) or (grad[m] < 0 and grad[m + 1] > 0)):
                y1.append(m - crossings[k])
        if(len(y1) > 1):
            pos = np.argmax(np.absolute(ecg_signal[y1]))
            y1 = y1[pos]
        y.append(y1)
    
    y2 = []
    
    if(len(y) > 4):
        if(y[1] < y[0] + 20 and y[-1] > end_qrs - 10):
            y2 = y[1 : len(y) - 1]
        elif(y[1]  < y[0] + 20):
            y2 = y[1 : ]
        elif(y[-1] > end_qrs - 10):
            y2 = y[: len(y) - 1]
        else:
            y2 = y
    else:
        y2 = y
    
    
    y3 = []
    
    if(bool(y2)):
        if(y2[0] <= start_qrs + 5 and y2[-1] >= end_qrs - 10):
            y3 = y2[1 : len(y2) - 1]
        elif(y2[0] <= start_qrs + 5):
            y3 = y2[1 : ]
        elif(y2[-1] >= end_qrs - 10):
            y3 = y2[: len(y) - 1]
        else:
            y3 = y
            
    # line 59
    
    y4 = []
    
    if(len(y3) > 4):
        values = ecg_signal[y3]
        
        max_1_pos = np.argmax(values)
        values[max_1_pos] = 0            
        max_1_pos = y3[max_1_pos]
        max_2_pos = np.amax(values)
        values[max_2_pos] = 0
        max_2_pos = y3[max_2_pos]
        
        min_1_pos = np.argmin(values)
        values[min_1_pos] = 0
        min_1_pos = y3[min_1_pos]
        min_2_pos = np.argmin(values)
        min_2_pos = y3[min_2_pos]
        
        y4 = [max_1_pos] + [max_2_pos] + [min_1_pos] + [min_2_pos]
        y4 = np.sort(y4).tolist()
    else:
        y4 = y3
        
    ## line 90
    
    y5 = []
    
    if(len(y4) > 2):
        if(y4[1] - y4[0] < 15):
            y5 = y4[1:]
        else:
            y5 = y4
    else:
        y5 = y4
    
    
    if(len(y5) == 3 and ecg_signal[y5[0]] > 0):
        if(y5[1] - y5[0] < 20):
            y5 = y5[1:]
    
    if(len(y5) == 1):
        if(ecg_signal[y5[0]] > 0):
            not2 = notching(ecg_signal, start_qrs + 5,y5)
            if(bool(not2)):
                val = np.amin(ecg_signal[not2])
                posm = np.argmin(ecg_signal[not2])
                if(val < 0):
                    posm = not2[posm]
                    y5 = [posm] + y5
            elif(not not2):
                not3 = notching(ecg_signal,y5,end_qrs)
                if(bool(not3)):
                    val = np.amin(ecg_signal[not3])
                    posm = np.argmin(ecg_signal[not3])
                    if(val < 0):
                        posm = not3[posm]
                        y5 = y5 + [posm]
    
    
    #line 130
    
    if((len(y5) == 3) and (ecg_signal[y5[0]] > 0 and ecg_signal[y5[1]] < 0 and ecg_signal[y5[2]] > 0)):
        values = [y5[0]] + [y5[2]]
        
        val1 = np.amax(ecg_signal[values])
        pos1 = np.argmax(ecg_signal[values])
        pos1 = values[pos1]
        
        val2 = np.amin(ecg_signal[values])
        pos2 = np.argmin(ecg_signal[values])
        pos2 = values[pos2]
        
        pos3 = y5[2]
        
        val3 = abs(ecg_signal[pos3])
        
        values1 = [pos1] + [pos3]
        
        val4 = np.amax(ecg_signal[values1])
        pos4 = np.argmax(ecg_signal[values1])
        pos4 = values1[pos4]
        
        val5 = np.amin(ecg_signal[values1])
        pos5 = np.argmin(ecg_signal[values1])
        pos5 = values1[pos5]
        
        if((val2 < 0.08*val1) or (val5 < 0.1*val4)):
            y5 = [y5[0]] + [y5[1]]
    
    # line 149
        
    if(len(y5) == 3 and ecg_signal[y5[0]] > 0 and ecg_signal[y5[1]] < 0 and ecg_signal[y5[2]] > 0):
        val = np.amin(ecg_signal[y5[2]:end_qrs])
        pos = np.argmin(ecg_signal[y5[2]:end_qrs])
        
        pos = pos + y5[2] - 1
        if(abs(val) > 0.75*abs(ecg_signal[y5[1]])):
            y5 = y5 + [pos]
    
    elif(len(y5) == 3 and ecg_signal[y5[0]] > 0 and ecg_signal[y5[1]] < 0 and ecg_signal[y5[2]] < 0):
        val = np.amax(ecg_signal[y5[1]:y5[2]])
        pos = np.argmax(ecg_signal[y5[2]:y5[2]])
        pos = pos + y5[1] - 1
        y5 = [y5[0]] + [pos] + [y5[2]]
    
    # line 164
    
    if(len(y5) == 3 and ecg_signal[y5[0]] > 0 and ecg_signal[y5[1]] > 0 and ecg_signal[y5[2]] < 0):
        val = np.amin(ecg_signal[y5[0]:y5[1]])
        pos = np.argmin(ecg_signal[y5[2]:y5[1]])
        
        pos = pos + y5[0] - 1
        if(abs(val) < 0.4*abs(ecg_signal[y5[2]])):
            y5 = [y5[0]] + [pos] + y5[1:]
    
    # line 174
    if(len(y5) == 4 and ecg_signal[y5[0]] > 0 and ecg_signal[y5[1]] < 0 and ecg_signal[y5[2]] > 0 and ecg_signal[y5[3]] < 0):
        values = ecg_signal[y5]
        
        max_1_pos = np.argmax(values)
        values[max_1_pos] = 0            
        max_1_pos = y5[max_1_pos]
        
        min_1_pos = np.amin(values)
        values[min_1_pos] = 0
        min_1_pos = y5[min_1_pos]
        
        min_2_pos = np.argmin(values)
        min_2_pos = y5[min_2_pos]
        
        absmax = max(abs(ecg_signal[max_1_pos]),abs(ecg_signal[min_1_pos]))
        absmin = min(abs(ecg_signal[max_1_pos]),abs(ecg_signal[min_1_pos]))
        
        absmin1 = max(abs(ecg_signal[min_2_pos]),abs(ecg_signal[min_1_pos]))
        absmin2 = min(abs(ecg_signal[min_2_pos]),abs(ecg_signal[min_1_pos]))
        
        if(absmin > 0.09*absmax and absmin2 > 0.24*absmin1):
            y5 = y5
        elif(ecg_signal[y5[0]] > 0 and ecg_signal[y5[1]] > ecg_signal[y5[3]]):
            y5 = [y5[0]] + [y5[2]] + [y5[3]]
        else:
            y5 = y5[0:3]
    
    if(len(y5) == 3 and ecg_signal[y5[0]] > 0 and ecg_signal[y5[1]] > 0 and ecg_signal[y5[2]] < 0 and ecg_signal[y5[1]] < 0.05):
        y5 = [y5[0]]
    
    if(len(y5) == 1 and ecg_signal[y5[0]] > 0):
        
        pos = zerocross(ecg_signal,start_qrs,y5)
        
        
    
        