# -*- coding: utf-8 -*-
"""
Created on Fri Jun 15 17:49:59 2018

@author: DELL
"""

import numpy as np

def notching (ecg_signal,on,off):
    cD = np.diff(ecg_signal)
    notch = []
    for i in range(on + 1, off - 2):
        if(((cD[i] > 0) and (cD[i + 1] < 0)) or ((cD[i] < 0) and (cD[i + 1] > 0))):
            if((ecg_signal[i] > 0.0015) or (ecg_signal[i] < -0.0015)):
                notch.append(i)
                
    notch1 = []
    
    for i in range(0,len(notch)):
        if(notch[i] < off - 20):
            notch1.append(notch[i])
    
    return notch1
    
    