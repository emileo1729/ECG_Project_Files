# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 13:46:12 2018

@author: DELL
"""
import numpy as np

def zerocross (ecg_signal,on,off):
    grad1 = np.diff(ecg_signal[on:off + 1])
    pos = []
    for i in range(0,len(grad1) - 1):
        if(((grad1[i] > 0) and (grad1[i + 1] < 0)) or ((grad1[i] < 0) and (grad1[i + 1] > 0))):
            pos.append(on + i-1)
    return pos