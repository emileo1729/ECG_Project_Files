# -*- coding: utf-8 -*-
"""
Created on Sat Jun 16 14:19:44 2018

@author: DELL
"""

import numpy as np
from movingslope import movingslope
import matplotlib.pyplot as plt
import pywt

def smooth(a,WSZ):
    # a: NumPy 1-D array containing the data to be smoothed
    # WSZ: smoothing window size needs, which must be odd number,
    # as in the original MATLAB implementation
    out0 = np.convolve(a,np.ones(WSZ,dtype=int),'valid')/WSZ    
    r = np.arange(1,WSZ-1,2)
    start = np.cumsum(a[:WSZ-1])[::2]/r
    stop = (np.cumsum(a[:-WSZ:-1])[::2]/r)[::-1]
    return np.concatenate((  start , out0, stop  ))


def hybrid_tdmg_FE_new(ecg_signal,time,R_peak_index):

    fs = 1000.0
    fsb = 1000.0
    
    rat = fs/fsb
    
    cA_l1, cD_l1 = pywt.dwt(ecg_signal,'haar');
    cA_l2, cD_l2 = pywt.dwt(cA_l1,'haar');
    cA_l3, cD_l3 = pywt.dwt(cA_l2,'haar');
    cA_l4, cD_l4 = pywt.dwt(cA_l3,'haar');
    cA_l5, cD_l5 = pywt.dwt(cA_l4,'haar');
    
    resol_l5 = 2**5
    
            
    min_val = np.amin(ecg_signal)
    min_pos = np.argmin(ecg_signal)
    max_val = np.amax(ecg_signal)
    max_pos = np.argmax(ecg_signal)
    
    window = 13
    
    grad1 = movingslope(ecg_signal,5,1)
    #plt.plot(grad1,'b')
    grad1 = smooth(grad1,window)
    #plt.plot(grad1,'g')
    grad2 = movingslope(grad1,5,1)
    #plt.plot(grad2,'r')
    grad2 = smooth(grad2,window)
    #plt.plot(grad2,'black')
    #plt.show()
    
    feat = 1.3*abs(grad1) + 1.1*abs(grad2)
    feat = np.power(feat,2)
    feat = smooth(feat,window)
    
    #plt.plot(feat)
    #plt.show()
    
    feat_max = np.amax(feat)
    feat_max_time = np.argmax(feat)
    
    if(feat_max < 200):
        limit_one_feat = 0.01*feat_max;
        limit_two_feat = 0.01*feat_max
    elif(feat_max < 850):
        limit_one_feat = 0.015*feat_max;
        limit_two_feat = 0.01*feat_max
    elif(feat_max < 1400):
        limit_one_feat = 0.009*feat_max;
        limit_two_feat = 0.009*feat_max
    elif(feat_max < 3000):
        limit_one_feat = 0.003*feat_max;
        limit_two_feat = 0.003*feat_max
    elif(feat_max < 10000):
        limit_one_feat = 0.0009*feat_max;
        limit_two_feat = 0.0009*feat_max
    else:
        limit_one_feat = 0.001*feat_max;
        limit_two_feat = 0.001*feat_max
    
    
    R_peak_pos_ref = R_peak_index
    R_peak_ref = ecg_signal[R_peak_index]
    
    if(feat_max_time > 200*rat):
        qw_start = round(200*rat)
    else:
        qw_start = feat_max_time - 1
    
    
    xa = np.where(feat[int(feat_max_time - qw_start) : int(feat_max_time - round(30*rat))] < limit_one_feat)[-1]
    if(xa.size):
        xa = xa[-1]
        start_qrs = feat_max_time - qw_start + xa - 1
    else:
        xa = -1
        start_qrs = -1
    
    k = 1
    limit_one_feat1 = limit_one_feat
    while(((xa == -1) and k < 10) or ((start_qrs != -1) and R_peak_pos_ref - start_qrs > 120*rat and k < 5)):
        limit_one_feat1 = limit_one_feat1 + 0.001*k*np.amax(feat)
        xa = np.where(feat[int(feat_max_time - qw_start) : int(feat_max_time - round(30*rat))] < limit_one_feat1)[-1]
        if(not xa.size):
            xa = -1
            start_qrs = -1
        else:
            xa = xa[-1]
            start_qrs = feat_max_time - qw_start + xa            
        k = k + 1
        
    
    k = 1
    start_qrs_er1 = start_qrs
        
    if(start_qrs_er1 - round(50*rat) > 0):
        xa_ind = np.where(feat[int(start_qrs_er1 - round(50*rat)) : int(start_qrs_er1)] > 20*limit_one_feat1)[-1]
        if(not xa_ind.size):
            xa_er = -1
        else:
            xa_er = xa_ind[-1]
    else:
        xa_ind = np.where(feat[int(1) : int(start_qrs_er1)] > 20*limit_one_feat1)[-1]
        if(not xa_ind.size):
            xa_er = -1
        else:
            xa_er = xa_ind[-1]
        
        
        
    while(xa_er != -1 and k < 3):
        if(start_qrs_er1 - round(50*rat) > 0):
            xa_ind = np.where(feat[int(feat_max_time - round(qw_start*rat)) : int(start_qrs_er1 - round(50*rat) + xa_er - round(15*rat))] < limit_one_feat1)[-1]
            if(not xa_ind.size):
                xa = -1
            else:
                xa = xa_ind[-1]
        
        else:
            xa_ind = np.where(feat[int(feat_max_time - round(qw_start*rat)) : int(1 + xa_er - round(15*rat))] < limit_one_feat1)[-1]
            if(not xa_ind.size):
                xa = -1
            else:
                xa = xa_ind[-1]
        
        if(xa != -1):
            start_qrs_er1 = feat_max_time - round(qw_start*rat) + xa
        
        k = k + 1
        
        if(start_qrs_er1 - round(50*rat) > 0):
            xa_ind = np.where(feat[int(start_qrs_er1 - round(50*rat)) : int(start_qrs_er1)] > 20*limit_one_feat1)[-1]
            if(not xa_ind.size):
                xa_er = -1
            else:
                xa_er = xa_ind[-1]
        
        else:
            xa_ind = np.where(feat[int(1) : int(start_qrs_er1)] > 20*limit_one_feat1)[-1]
            if(not xa_ind.size):
                xa_er = -1
            else:
                xa_er = xa_ind[-1]
        
            
        if(R_peak_pos_ref - start_qrs_er1 < 160*rat):
            start_qrs = start_qrs_er1
    
    if(start_qrs > R_peak_pos_ref):
        start_qrs = R_peak_pos_ref - 2
        
    if(not start_qrs):
        start_qrs = 1
        
    ######## QRS refinement based on amplitude
    
    zcf = 0
    if(ecg_signal[int(start_qrs)] < -0.1 and start_qrs - 25 > 0):
        
        for i in  range(start_qrs - 25,start_qrs):
            if(ecg_signal[i] > 0 and ecg_signal[i + 1] < 0):
                start_qrs = i
                zcf = 1
        
        if(zcf == 0):
            start_qrs = start_qrs - 15
    
    
    zcf1 = 0
    if(ecg_signal[int(start_qrs)] > 1.1 and start_qrs - 90 > 0):
    
        for i in range(start_qrs - 90,start_qrs):
            if(ecg_signal[i] < 0 and ecg_signal[i + 1] > 0):
                start_qrs = i
                zcf1 = 1
        
        if(zcf1 == 0):
            start_qrs = start_qrs - 85
    
    elif(ecg_signal[int(start_qrs)] > 0.14 and start_qrs - 25 > 0):
        for i in range(start_qrs - 25,start_qrs):
            if(ecg_signal[i] < 0 and ecg_signal[i + 1] > 0):
                start_qrs = i
                zcf1 = 1
        
        if(zcf1 == 0):
            start_qrs = start_qrs - 25
    
    ###############################################
    ################## End Of QRS #################
    
    if(feat_max_time + (200*rat) < len(ecg_signal)):
        qw_end = int(round(200*rat))
    else:
        qw_end = len(ecg_signal) - feat_max_time - 1
    
    
    xb = np.where(feat[int(feat_max_time + round(30*rat)) : int(feat_max_time + qw_end)] < limit_two_feat)[-1]
    if(xb.size):
        xb = xb[0]
        end_qrs = feat_max_time + xb + round(30*rat)
    else:
        xb = -1
        end_qrs = -1   
    
    print feat_max_time
    
    k = 1
    # line 226
    
    limit_two_feat1 = limit_two_feat
    
    while(((xb == -1) and (k < 15)) or ((end_qrs != -1) and (end_qrs - R_peak_pos_ref > 160*rat) and (k < 10))):
        
        limit_two_feat1 = limit_two_feat1 - 0.001*k*np.amax(feat)
        xb = np.where(feat[int(feat_max_time + round(30*rat)) : int(feat_max_time + qw_end)] < limit_two_feat)[-1]
        if(xb.size):
            xb = xb[0]
            end_qrs = feat_max_time + xb + round(30*rat)
        else:
            xb = -1
            end_qrs = -1
        k = k + 1
    
    
    #line276
    
    end_qrs_er1 = end_qrs
    
    k = 1
    
    if(end_qrs_er1 + round(50*rat) < len(feat)):
        xb_er = np.where(feat[int(end_qrs_er1) : int(end_qrs_er1 + round(50*rat))] > limit_two_feat1)[-1]
        if(xb_er.size):
            xb_er = xb_er[0]
        else:
            xb_er = -1
    else:
        xb_er = np.where(feat[int(end_qrs_er1) : int(len(feat))] > limit_two_feat1)[-1]
        if(xb_er.size):
            xb_er = xb_er[0]
        else:
            xb_er = -1
        
    #line 287
    
    while((xb_er != -1) and (k < 3)):
        xb = np.where(feat[int(end_qrs_er1 + xb_er + round(15*rat)) : int(feat_max_time + round(qw_end*rat))] > limit_two_feat1)[-1]
        if(xb.size):
            xb = xb[0]
        else:
            xb = -1
        
        if(xb != -1):
            end_qrs_er1 = end_qrs_er1 + xb_er + round(15*rat) + xb
        
        k = k + 1
        
        if(end_qrs_er1 + round(50*rat) < len(feat)):
            xb_er = np.where(feat[int(end_qrs_er1) : int(end_qrs_er1 + round(50*rat))] > limit_two_feat1)[-1]
            if(xb_er.size):
                xb_er = xb_er[0]
            else:
                xb_er = -1
        else:
            xb_er = np.where(feat[int(end_qrs_er1) : int(len(feat))] > limit_two_feat1)[-1]
            if(xb_er.size):
                xb_er = xb_er[0]
            else:
                xb_er = -1
            
        if(end_qrs_er1 - R_peak_pos_ref < round(180*rat)):
            end_qrs = end_qrs_er1

    
    if(end_qrs == -1):
        end_qrs = 1
        
    # line 316
    
    # QRS end refinement based on amplitudes
    
    zf = 0
    
    if(ecg_signal[int(end_qrs)] < -0.9 and ecg_signal[R_peak_index] < 0.5):
        for i in range(end_qrs, end_qrs + 101):
            if(ecg_signal[i] < 0 and ecg_signal[i + 1]> 0):
                end_qrs = i
                zf = 1
        if(zf == 0):
            end_qrs = end_qrs + 100
    elif(ecg_signal[int(end_qrs)] < -0.55 and ecg_signal[int(R_peak_index)] < 0.5):
        for i in range(end_qrs,end_qrs + 45):
            if(ecg_signal[i] < 0 and ecg_signal[i + 1]> 0):
                end_qrs = i
                zf = 1
        if(zf == 0):
            end_qrs = end_qrs + 45
    
    zf2 = 0
    
    if(ecg_signal[int(end_qrs)] > 0.45):
        for i in range(end_qrs - 85,end_qrs):
            if(ecg_signal[i] < 0 and ecg_signal[i + 1] > 0):
                end_qrs = i
                zf2 = 1
        if(zf2 == 0):
            end_qrs = end_qrs - 85
    
    #######################################
    
    x1 = []
    y1 = []
    
    x1.append(start_qrs)
    y1.append(feat[int(start_qrs)])
    
    x2 = []
    y2 = []
    
    x2.append(end_qrs)
    y2.append(feat[int(end_qrs)])
    
    plt.plot(feat)
    plt.plot(x1,y1,'ro')
    plt.plot(x2,y2,'go')
    plt.show()
    
    
    print "start_qrs = ",
    print start_qrs
    
    print "end_qrs = ",
    print end_qrs