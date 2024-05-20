# -*- coding: utf-8 -*-
"""
Created on Tue Feb 6 11:36:57 2024

@authors: Eleanor Church & Kate Sewart

06/05/24

Takes the volunteer motion data and tests on first 7s then tests on rest of data to find
resp gating for end of exhale
DOES have option to flip axis
DOES have end of pause
USES other volunteer data & plots it

finds length of pause + averages
"""

#%% IMPORTS

import numpy as np
import matplotlib.pyplot as plt
from datetime import date
import matplotlib.style as mplstyle
mplstyle.use('fast')

#%% OPEN DATA FILES

volunteernum = 20
folder = 'C:/Users/HP EliteBook/Documents/4th year/Motion in upright MRI/DATA_findexhale/Volunteer'+str(volunteernum)+'/'

# open the three translations, four quarternion parameters and time data
TX=np.loadtxt(folder+"TX.csv")
TY=np.loadtxt(folder+"TY.csv")
TZ=np.loadtxt(folder+"TZ.csv")
QX=np.loadtxt(folder+"QX.csv")
QY=np.loadtxt(folder+"QY.csv")
QZ=np.loadtxt(folder+"QZ.csv")
QW=np.loadtxt(folder+"QW.csv")
time=np.loadtxt(folder+"Time.csv")

#%% DEFINE FUNCTIONS

def AverageVsMedian(DATA_NORM):
    if np.median(DATA_NORM) > np.mean(DATA_NORM):
        Multiplication_factor = -1
    else:
        Multiplication_factor = +1
    return Multiplication_factor
# compares the mean and median values in the normalised data.  
# if the median is greater than the mean, the greatest number of plateaus
# in the data occur closer to the maximum than the minimum so the data is inverted

def ZeroesVsOnes(DATA_INV_NORM_1_ROUND):
    if len(DATA_INV_NORM_1_ROUND)-np.sum(DATA_INV_NORM_1_ROUND) > np.sum(DATA_INV_NORM_1_ROUND):
        Multiplication_factor = 1
    else:
        Multiplication_factor = -1
    return Multiplication_factor
# compares the numbers of zeros and ones in the normalised data.
# the number of ones being greater than the number of zeros corresponds to
# the plateau being at the top of the plot so the data is inverted in this case to
# bring the plateau to the bottom

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    print(idx, array[idx])
    return idx
# finds the index of an array nearest to a value given
# For this it will be 7s, but could be changed to work on different lengths of cycle

#%% NORMALISE DATA AND CHECK IF IN/OUT OF PHASE

def data_norm(data):
    data_norm = (data-np.min(data))/(np.max(data)-np.min(data))
    data_norm_round = np.round(data_norm, decimals=0)
    if np.argmax([data_norm_round]) == 0:
        Multiplication_factor = -1
    else:
        Multiplication_factor = 1
    #invert data
    data = data * Multiplication_factor
    return data
 # data are normalised between 0 and 1 and inverted if
 # the index of the maximum value in the array is zero

TX = data_norm(TX)
TY = data_norm(TY)
TZ = data_norm(TZ)
QX = data_norm(QX)
QY = data_norm(QY)
QZ = data_norm(QZ)
QW = data_norm(QW)

def CharacteriseResp(x, T):
    T = np.transpose(T)  
    t2 = np.concatenate((T, T + (T[-1] - T[0]), T + 2 * (T[-1] - T[0]), T + 3 * (T[-1] - T[0])))
    X = x - np.mean(x, axis=0)
    x2 = np.zeros((len(t2)))
    x2[:X.shape[0]] = X * np.exp(-np.transpose(T) / 20)
   
    Y = np.fft.fft(x2.squeeze())
    Y = np.abs(Y)
   
    nd = len(t2)
    f = np.arange(nd) / nd
   
    Ymax = np.max(Y[np.where((f >= 0.1) & (f <= 0.4))])
    if Ymax == 0 or np.where(Y == Ymax)[0][0] < 0.2 or np.where(Y == Ymax)[0][0] >= 0.4:
        print('Max of the spectra not found OR smaller than 0.2 Hz OR greater than 0.4 Hz: f_resp will be the default one.')
        f_resp = 0.3
    else:
        f_resp = f[np.where(Y == Ymax)][0]
    # searches for a respiration frequecny between 0.2Hz and 0.4Hz
    # uses 0.3Hz as default if none found
    
    respCycle_seconds = 1 / f_resp
    respCycle_frames = np.where(T >= (T[0] + respCycle_seconds))[0][0]
   
    return respCycle_frames
# takes the motion data, x, and the time data, T, to obtain an expected number
# of points for one respiration cycle.


#%% TRAINING

def training(TX,TY,TZ,QX,QY,QZ,QW,y):
    # take first 7 seconds of each parameter
    idx = find_nearest(time,7)
    TX_sample = TX[:idx]
    TY_sample = TY[:idx]
    TZ_sample = TZ[:idx]
    QX_sample = QX[:idx]
    QY_sample = QY[:idx]
    QZ_sample = QZ[:idx]
    QW_sample = QW[:idx]

    # find the standard deviation of each parameter
    TX_std = np.std(TX_sample)
    TY_std = np.std(TY_sample)
    TZ_std = np.std(TZ_sample)
    QX_std = np.std(QX_sample)
    QY_std = np.std(QY_sample)
    QZ_std = np.std(QZ_sample)
    QW_std = np.std(QW_sample)
   
    # make an array of standard deviations
    sd_array = np.array([TX_std, TY_std, TZ_std, QX_std, QY_std, QZ_std, QW_std])
    motion_parameters = np.array(['TX', 'TY', 'TZ', 'QX', 'QY', 'QZ', 'QW'])

    # put the values in a new array in descending order
    sorted_array_desc = np.sort(sd_array)[::-1]
   
    # take the first two values
    best_vals = sorted_array_desc[0:2]
    print(best_vals)
   
    # find associated motion parameters and sort from smallest to largest
    sorted_motionparams = [x for _,x in sorted(zip(sd_array,motion_parameters))]
    best_params = sorted_motionparams[5:7]
    print(best_params)
   
    # find associated norms with highest standard deviations
    sorted_norms = [x for _,x in sorted(zip(sd_array,np.array([TX,TY,TZ,QX,QY,QZ,QW])))]
    best_norms = sorted_norms[5:7]
   
    # comparison between mean and median to get multiplication factor
    mult_1 = AverageVsMedian(best_norms[0])
    mult_2 = AverageVsMedian(best_norms[1])
   
    # invert data if necessary
    inv_1 = mult_1*best_norms[0]
    inv_2 = mult_2*best_norms[1]
   
    # renormalise data
    inv_1a = data_norm(inv_1)
    inv_2a = data_norm(inv_2)


    # find if in/out of phase
    if max(inv_1a+inv_2a)<1.6:
        print('Out of phase')
       
        # approximate to nearest integer
        inv_norm_round_1 = np.round(inv_1a)
        inv_norm_round_2 = np.round(inv_2a)
       
        ## as before ##
        mult_1 = ZeroesVsOnes(inv_norm_round_1[0])
        mult_2 = ZeroesVsOnes(inv_norm_round_2[1])
   
        inv_1 = mult_1*best_norms[0]
        inv_2 = mult_2*best_norms[1]
   
        inv_1b = data_norm(inv_1)
        inv_2b = data_norm(inv_2)
        ## ##
       
        # redefines inverted data to be parsed in later code
        inv_1a = inv_1b
        inv_2a = inv_2b
   
    # find number of points in resp cycle
    data_inv_array = ([inv_1, inv_2])
    sum_data_inv = np.sum(data_inv_array, axis=0) # adds two most significant parameters together
    time_subset = time
    RespCycle = CharacteriseResp(sum_data_inv, time_subset) # num points in one resp cycle
    RespCycle_points = round(0.9*RespCycle) # rough minimum of cycle length
    motion_parameters = np.array([TX, TY, TZ, QX, QY, QZ, QW])
   
    sorted_motionparams = [x for _,x in sorted(zip(sd_array,motion_parameters))]
    data = sorted_motionparams[5:7]
   
    # invert data
    data[0] = mult_1*data[0]
    data[1] = mult_2*data[1]
   
    # sum the signals along the second axis - summing each row
    data = np.sum(data, axis=0)
   
    if y == 'y':
        fig= plt.figure(figsize=(10,5))
        ax1 = fig.gca()
        ax1.plot(time,data)
        plt.waitforbuttonpress()
        flip_y = input("Do you want to flip the y-axis? (yes/no): ")
    else:
        flip_y='no'
    if flip_y.lower() == 'yes':
        data = -1*data
        plt.cla()
        ax1.plot(time,data)
    # checks by human eye if axis needs flipping again
    # so that plateau is on the bottom
   
    return RespCycle, RespCycle_points, data



#%% TESTING

def testing(RespCycle, RespCycle_points, data, grad_std):
    trigger = 0 # counts the number of triggers
    dt = np.mean(np.diff(time))
    RespCycle_sec =2
    flag_top =0
    timesince=1000
    istart_plot = 4*RespCycle
    grdlimL = np.round(RespCycle_sec/4/dt) # length of negative gradient to be assumed exhale
    grdlimS = (grdlimL*grad_std) # length of gradient uptick to be considered end of exhale
   
    # defines arrays to be filled with triggers and data to plot
    trigger_sec = np.zeros((data[:, 0], data[:, 1]))
    Trigger_Value = np.zeros((data[:, 0], data[:, 1]))
    data_1 = np.zeros((data[:, 0], data[:, 1]))
    data_2 = data_1
   
    # empty arrays to update running min, max and mean
    running_min = np.zeros((data[:, 0], data[:, 1]))
    running_mean = np.zeros((data[:, 0], data[:, 1]))
    running_max = np.zeros((data[:, 0], data[:, 1]))
   
    # empty arrays to update local min, max and mean
    local_min = np.zeros((data[:, 0], data[:, 1]))
    local_mean = np.zeros((data[:, 0], data[:, 1]))
    local_max = np.zeros((data[:, 0], data[:, 1]))


    for j in range(RespCycle+1,time.size): # iterates over each respiration cycle
        # performs drift removal
        data_1[j] = np.sum(data[j-8:j])/9
        local_mean[j] = np.mean(data_1[j-RespCycle:j])
        data_2[j] = data_1[j]-local_mean[j]
       
        # updates running min, max and mean
        if j>istart_plot:
            running_min[j] = data_2[istart_plot:j].min()
            running_mean[j] = np.mean(data_2[istart_plot:j])
            running_max[j] = data_2[istart_plot:j].max()
       
        # updates local min, max and mean
        local_min[j] = data_2[j-RespCycle:j,0].min()
        local_mean[j] = np.mean(data_2[j-RespCycle:j,0])
        local_max[j] = data_2[j-RespCycle:j,0].max()
       
        # if statement triggers if data is at a peak
        if data_2[j,0]>local_mean[j,0]:
            flag_top=1
       
        # resets grdlimL and grdlimS if j is too small
        if j <= grdlimL:
            grdlimL = RespCycle
        if j <= grdlimS:
            grdlimS = RespCycle
   
   
        gradLong = data[j]- data[int(j - grdlimL)]
        gradShort = data[j] - data[int(j - grdlimS)]
        # gradient must be negative and uptick recently for a trigger.
       
        # calculates time since last trigger
        if trigger > 0:
            timesince = (time[j] - trigger_sec[trigger])[0]

        # defines all conditions needed for there to be a trigger
        if (flag_top == 1 and gradLong < 0 and gradShort > gradLong * grad_std  and timesince > RespCycle_sec and data_2[j,0] < local_mean[j,0] and abs(data_2[j,0] - local_mean[j,0]) > 0.4 * abs(local_min[j,0] - local_mean[j,0])):
            trigger += 1  # increases trigger count
            trigger_sec[trigger] = time[j]
            Trigger_Value[trigger] = data_2[j]
            flag_top = 0

    print(trigger)
    return trigger_sec


#%% CALLING FUNCTIONS

# triggering start of pause first
RespCycle, RespCycle_points, data = training(TX,TY,TZ,QX,QY,QZ,QW,'y')

# finds sampling frequency for dataset
sampling_frequency = len(data)/len(time)

# finds the number of data points within interval of 0.1 seconds
data_points_in_interval = int(0.1 * sampling_frequency)

# groups the data and calculate the gradients across each group
data_groups = data.reshape(-1, data_points_in_interval)
gradients = np.diff(data_groups, axis=1)

# makes a 1D array of gradients
flattened_gradients = gradients.flatten()

# calculates the standard deviation of gradients over intervals of 0.1 seconds
std_dev_gradients = np.std(flattened_gradients)

print("Standard deviation of gradients over 0.1 second groups:", std_dev_gradients)


# applying triggers to remaining data
trigger_sec=testing(RespCycle, RespCycle_points, data, std_dev_gradients)

# triggering ends of pauses 
TX_flip = np.flip(TX)
TY_flip = np.flip(TY)
TZ_flip = np.flip(TZ)
QX_flip = np.flip(QX)
QY_flip = np.flip(QY)
QZ_flip = np.flip(QZ)
QW_flip = np.flip(QW)

RespCycle_end, RespCycle_points_end, data_end = training(TX_flip,TY_flip,TZ_flip,QX_flip,QY_flip,QZ_flip,QW_flip,'n')
trigger_sec_end = max(time)-testing(RespCycle_end, RespCycle_points_end, data_end,std_dev_gradients)


#%% PLOTTING THE DATA AND TRIGGERS

# sets up figure
fig= plt.figure(figsize=(20,10))
ax1 = plt.subplot(2,1,1)
ax2 = plt.subplot(2,1,2)

# plots for the start of pause triggers
ax1.plot(time,data)
ax1.vlines(x=trigger_sec, ymin=min(data), ymax=max(data), color='r')
# plots for the end of pause triggers
ax2.plot(time,data)
ax2.vlines(x=trigger_sec_end, ymin=min(data), ymax=max(data), color='g')

# automatically saves figure to file (change filedir to your own)
date = date.today()
filedir = 'C:/Users/HP EliteBook/Documents/4th year/Motion in upright MRI/Graphs/'
plt.savefig(filedir+str(date)+' Volunteer'+str(volunteernum)+' with end of exhale and stats')

#%% FINDING AVERAGE EXHALATION LENGTH FOR THE VOLUNTEER

# defines variables and array to be used in iteratiom
exhale_length = []
i_prev = 0 # initial previous start of pause
j_prev = 0.0000001 # initial previous end of pause

# redefining trigger arrays to just the locations of triggers
begin_exhale = trigger_sec[:,0][trigger_sec[:,0]!=0]
end_exhale = np.flip(trigger_sec_end[:,0][trigger_sec_end[:,0]!=max(time)],0)

# iterates through trigger arrays to calculate average
for i in begin_exhale:
    for j in end_exhale:
        if i<j and i_prev<j_prev: # prevents two triggers in same cycle
            print(i,j)
            exhale_length += [j-i] # all exhalation lengths
            break
        j_prev = j
    i_prev = i

# finds mean exhalation length
avg_exhale = np.mean(exhale_length)

print("lengths of pause periods" + exhale_length)
print("average length of pause period" + avg_exhale)