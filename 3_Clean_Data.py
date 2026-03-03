import numpy as np
import os
import csv
import matplotlib.pyplot as plt
from tkinter import filedialog, messagebox
import tkinter as tk

def graphic_thres_mark(values,t_up,t_dw,low):
    '''

    :param values: Array of value to clean
    :param t_up: Rising time to be cleaned before, once the value is above threshold
    :param t_dw: Falling time to be cleaned after, once the value is below threshold
    :param low: if True, the threshold work as a lower threshold
    :return: removal mark (0 -> keep it, 1 -> put NaN instead)
    '''
    global Select_data,idt

    t_up = t_up/3600.
    t_dw = t_dw/3600.

    N = len(values)
    answer = False
    while not answer:
        plt.plot(Select_data[:,idt],values)
        if not low: plt.title('Clean spikes by drawing an upper threshold value')
        else: plt.title('Clean spikes by drawing a lower threshold value')

        for t in t_break:
            plt.axvline(x=t, ymin=0, ymax = 1, linewidth=1, color='green')
        for t in t_change:
            plt.axvline(x=t, ymin=0, ymax = 1, linewidth=1, color='red')
        plt.ylim([np.nanmean(values)-2*np.nanstd(values),np.nanmean(values)+4*np.nanstd(values)])
        pts = plt.ginput(n=-1,timeout=-1)
        plt.show()
        X = [np.min(Select_data[:,idt])]
        Y = [pts[0][1]]
        for (x,y) in pts:
            if x> np.min(Select_data[:,idt]) and x<np.max(Select_data[:,idt]):
                X.append(x)
                Y.append(y)
        X.append(np.max(Select_data[:,idt]))
        Y.append(pts[-1][1])

        Y = np.array(Y)[np.argsort(X)]
        X = np.array(X)[np.argsort(X)]

        thres = np.interp(Select_data[:,idt],X,Y)

        plt.plot(Select_data[:,idt],values)
        plt.title('Here is the drawn threshold')
        plt.plot(Select_data[:,idt],thres,c='black')
        plt.show()

        application_window = tk.Tk()
        answer = messagebox.askyesno("Question","Happy with the threshold?", parent = application_window)
        application_window.destroy()

    if low==True:
        values *=-1
        thres *=-1

    rmv_mark = np.zeros(N)
    for i in range(N):
        if values[i]>thres[i]:
            rmv_mark[i] = 1
        if values[i-1]<thres[i-1] and values[i]>thres[i]:
            for j in range(0,i):
                if abs(Select_data[j,idt]-Select_data[i,idt])<t_up:
                    rmv_mark[j] = 1
        if values[i-1]>thres[i-1] and values[i]<thres[i]:
            for j in range(i,N):
                if abs(Select_data[j,idt]-Select_data[i,idt])<t_dw:
                    rmv_mark[j] = 1

    if low==True:
        values *= -1
    return rmv_mark

def HK_clean(idc):
    '''

    :param idc: Index of the compound in Select_data
    :return: removal mark (0 -> keep it, 1 -> put NaN instead)
    '''
    global Select_data,headers_raw
    print('Cleaning using HK parameters')

    idp = headers_raw.index('P_from_Ctrl')
    idf = headers_raw.index('F_from_Ctrl')
    idpcfa = headers_raw.index('P_CFA')
#    idmelt = headers_raw.index('MeltRate')
    p0,p1 = 19.8,20.2
    pcfa0,pcfa1 = 200.,800.
    f0,f1 = 0.65,2.
#    melt0, melt1 = -1, 8.

    rmv_mark = np.zeros(len(Select_data[:,idc]))
    # Clean using SARA Pressure
    rmv_mark = np.where(np.logical_or(Select_data[:,idp]<p0,Select_data[:,idp]>p1),1,0) + rmv_mark
    # Clean using Module Pressure
    rmv_mark = np.where(np.logical_or(Select_data[:, idpcfa] < pcfa0, Select_data[:, idpcfa] > pcfa1), 1, 0) + rmv_mark
    # Clean using Flow
    rmv_mark = np.where(np.logical_or(Select_data[:, idf] < f0, Select_data[:, idf] > f1), 1, 0) + rmv_mark
    # Clean Melt Rate
    # rmv_mark = np.where(np.logical_or(Select_data[:, idmelt] < melt0, Select_data[:, idmelt] > melt1), 1, 0) + rmv_mark

    rmv_mark = np.clip(rmv_mark,0,1)

    return rmv_mark

def Kero_clean():

    global Select_data,idt

    #file_kero = easygui.fileopenbox(msg='Choose a Kerosene file')
    file_kero = '/media/fourteau/KevinF/Data/CFA-vostok/Kero-infil-0222.txt'
    if file_kero is None: raise Exception('No file selected...')
    starts,stops = np.loadtxt(file_kero,unpack=True)
    T_shift = starts[0]
    starts = starts[1:] - T_shift/3600.
    stops = stops[1:] - T_shift/3600.
    rmv_mark = np.zeros(len(Select_data[:,idt]),dtype=bool)
    for start,stop in zip(starts,stops):
        rmv_mark = np.where(np.logical_and(Select_data[:,idt]>=start,Select_data[:,idt]<=stop),1,rmv_mark)
        wh = np.where(np.logical_and(Select_data[:, idt] >= start, Select_data[:, idt] <= stop))
    return rmv_mark

def Der_clean(idc,t_width=30,t_up=20,t_dw=20,low=False):
    '''
    Clean using spikes in the derivative of the compound
    :param idc: Index of the compound in Select_data
    :param t_width: Temporal width to compute the rmsd
    :param t_up: Rising time to be cleaned before, once the value is above threshold
    :param t_dw: Falling time to be cleaned after, once the value is below threshold
    :param low: if True, the threshold work as a lower threshold
    :return: removal mark (0 -> keep it, 1 -> put NaN instead)
    '''
    t_width=t_width/3600.

    global Select_data, idt
    print('Cleaning using derivative')

    N = len(Select_data[:,idt])
    der = (Select_data[1:,idc]-Select_data[0:-1,idc])/(Select_data[1:,idt]-Select_data[0:-1,idt])
    der = np.append(der,0)

    for i in range(N):
        if Select_data[i,idt]-Select_data[0,idt]> t_width:
            i_st = i
            break

    rmsd = np.zeros(N)
    jlast = 0
    for i in range(i_st,N):
        tmp = np.array([])
        pts = np.array([])
        for j in range(jlast,i):
            if abs(Select_data[i,idt]-Select_data[j,idt])<=t_width:
                tmp = np.append(tmp,der[j])
                pts = np.append(pts,j)
        if len(pts)==0:
            plop = 1
            rmsd[i] = np.NaN
        if len(pts)>0:
            jlast = int(np.min(pts))
            rmsd[i] = np.std(tmp)
    rmsd[0:i_st] = rmsd[i_st:2*i_st]

    rmv_mark = graphic_thres_mark(rmsd,t_up,t_dw,low)

    return rmv_mark

def Thres_clean(idc,t_up=20,t_dw=20,low=False):
    '''
    Clean using spikes in the concentration of the compound
    :param idc: Index of the compound in Select_data
    :param t_up: Rising time to be cleaned before, once the value is above threshold
    :param t_dw: Falling time to be cleaned after, once the value is below threshold
    :param low: if True, the threshold work as a lower threshold
    :return: removal mark (0 -> keep it, 1 -> put NaN instead)
    '''

    print('Cleaning using concentrations')

    rmv_mark = graphic_thres_mark(Select_data[:,idc],t_up,t_dw,low)
    return rmv_mark

def Depth_clean(idc,dmin,dmax):
    '''
    Clean using a min and max depth
    :param idc: Index of the compound in Select_data
    :param dmin: Min depth
    :param dmax: Max depth
    :return: removal mark (0 -> keep it, 1 -> put NaN instead)
    '''
    global Depth
    print('Cleaning using min and max depths')

    rmv_mark = np.logical_or(Depth<dmin,Depth>dmax)
    return rmv_mark

def Local_clean(idc):

    global Select_data,idt
    t_win = 600/60. #largeur plot successif pour nettoyage data
    t0 = min(Select_data[:,idt])-10/60.
    t = t0
    rmv_mark = np.zeros(len(Select_data[:,idt]))
    while not t>max(Select_data[:,idt]):
        good = False
        wh_time = np.where(np.logical_and(t<=Select_data[:,idt], Select_data[:,idt]<=t+t_win))
        if np.all(np.isnan(Select_data[wh_time,idc])):
            t += t_win
            good = True
        while not good:
            print(t)
        
            plt.plot(Select_data[:,idt],Select_data[:,idc])
            for tb in t_break:
                plt.axvline(x=tb, ymin=0, ymax = 1, linewidth=1, color='green')
            for tb in t_change:
                plt.axvline(x=tb, ymin=0, ymax = 1, linewidth=1, color='red')
            plt.xlim([t,t+t_win])
            wh = np.where(np.logical_and(Select_data[:,idt]>=t,Select_data[:,idt]<=t+t_win))

            plt.ylim([min(Select_data[wh,idc][0]),min(Select_data[wh,idc][0])+200])
            pts = plt.ginput(n=-1,timeout=-1)
            plt.show()
            X = []
            for (x,y) in pts:
                X.append(x)
            if len(X)%2 !=0 :
                print('The number of points should be even')
                continue
            if len(X) !=0 :
                X = np.array(X)[np.argsort(X)]
                for i in range(0,len(X),2):
                    x0,x1 = X[i], X[i+1]
                    rmv_mark[np.where(np.logical_and(Select_data[:,idt]>=x0,Select_data[:,idt]<=x1))] = 1
            t += t_win
            good = True
    plt.plot(rmv_mark)
    plt.show()
    return rmv_mark

def get_rmv_mark(idc_key,hk=True,der=False,thres=False,depth=False):
    '''

    :param idc: Index of the compound in Select_data
    :param hk: if true use HK_clean
    :param der: if true use Der_clean
    :param thres: if true use Thres_clean
    :param depth: if true use Depth_clean
    :return: Nothing but update Select_data
    '''
    global Select_data,Depth,idt
    # Create a removal mark (if mark=one -> value=NaN)
    rmv_mark = np.zeros(len(Select_data[:,idt]))
    rmv_mark = np.clip(rmv_mark + Local_clean(idc_key), 0, 1)
    if hk : rmv_mark = np.clip(rmv_mark+HK_clean(idc_key),0,1)
    if der :rmv_mark = np.clip(rmv_mark+Der_clean(idc_key),0,1)
    if thres : rmv_mark = np.clip(rmv_mark+Thres_clean(idc_key),0,1)
    if depth : rmv_mark = np.clip(rmv_mark+Depth_clean(idc_key,85000,95000),0,1)
    return rmv_mark

def clean_main(d_compounds,wd=None):
    '''
    Clean data
    :param compounds: Array of the compounds to be cleaned (ex : ['CH4','CO'])
    :param wd: Working Directory. If None, one is chosen with gui
    :return: nothing, but save cleaned data
    '''

    # Get working directory and Date
    if wd is None:
        application_window = tk.Tk()
        wd = filedialog.askdirectory(title='Select day directory', parent=application_window)
        application_window.destroy()
    if not wd: raise Exception('No directory selected...')

    Date = os.path.basename(wd)
    Year = Date[0:4]
    Month = Date[4:6]
    Day = Date[6:8]

    # Get file
   #file = os.path.join(wd,'SARA_cleaned_break_no_calib_'+Date+'.txt')
    file = os.path.join(wd,'SARA_cleaning_'+Date+'.txt')

    # Load Spreadsheets
    Daily_cfa_Eventlog_File = os.path.join(wd,'CFA_Event_Log_'+Day+'_'+Month+'_'+Year+'.xls')
    Depth_File=os.path.join(wd,'RUN-'+Year+'-'+Month+'-'+Day+'.txt')

    # Put Data and idt in global scope
    global Select_data, Depth, idt,t_break,t_change,headers_raw

    # Get Depth and Time scale from Run file
    #Time_run, Depth_run, Melt_run = np.loadtxt(Depth_File,unpack=True)
    Time_str, Depth_run, Melt_run = np.loadtxt(Depth_File, dtype=str, unpack=True,skiprows=1)
    Time_run = np.array([int(t.split(':')[0]) +
                         int(t.split(':')[1])/60 +
                         int(t.split(':')[2])/3600
                         for t in Time_str], dtype=float)

    Depth_run = Depth_run.astype(float)
    Melt_run = Melt_run.astype(float)

    # Get min and max time of Run file
    mint = np.min(Time_run)
    maxt = np.max(Time_run)

    # Get break and core change Melter times
    t_break = np.array([])
    t_change = np.array([])
    with open(Daily_cfa_Eventlog_File, 'r') as csvEvent:
        csvreader = csv.reader(csvEvent, delimiter='\t')
        for r in csvreader:
            if str(r[3]) == 'Core Break':
                t_break = np.append(t_break, \
                            float(r[2].split(':')[0]) + float(r[2].split(':')[1])/60. + float(r[2].split(':')[2])/3600.)
            if str(r[3]) == 'Change of core':
                t_change = np.append(t_change,\
                           float(r[2].split(':')[0]) + float(r[2].split(':')[1])/60. + float(r[2].split(':')[2])/3600.)

    # Load data and read header
    Raw_data = np.loadtxt(file,delimiter='\t')
    f=open(file)
    headers_raw = f.readline().replace('# ','').rstrip().split('\t')
    f.close()

    # Get Index of interest using header
    idt = headers_raw.index('Time')

    # Select data that are only in the good range of time
    Select_data = Raw_data[np.where(np.logical_and(Raw_data[:,idt]>mint,Raw_data[:,idt]<maxt))]
    # Interpolate depth of Run file
    Depth = np.interp(Select_data[:,idt],Time_run,Depth_run)
    # Smooth depth scale
    window = 10
    # Depth=Depth.rolling(windows=10).mean()
    # Depth = np.convolve(Depth, np.ones(window)/window, mode='same')
    #Depth = rolling_mean(Depth, window=10, min_periods=1, center=True)

    # Clean concentrations
    for c_key in d_compounds.keys():
        idc_key = headers_raw.index(c_key)
        rmv_mark = get_rmv_mark(idc_key,hk=False,der=False,thres=False,depth=False)
        for c in d_compounds[c_key]:
            idc = headers_raw.index(c)
            Select_data[np.where(rmv_mark==1),idc] = np.nan


    # Create the output Array
    Out_array = np.column_stack((Depth,Select_data))
    wh = np.where(~np.logical_or(np.append(np.diff(Depth),1)==0, np.isnan(Depth)))
    Out_array = Out_array[wh,:]
    Out_array = Out_array[0]

    # Save data using numpy.savetxt
    outfile = os.path.join(wd,'SARA_cleaned_break_no_calib_' + Date + '.txt')
    header = 'Depth'
    for h in headers_raw:
        header = header + '\t' + h
    np.savetxt(outfile,Out_array,delimiter='\t',header=header,fmt='%.6e')

    # Plot cleaned data
    already = set()
    for group in d_compounds.values():
        for c in group:
            if c in already : continue
            already.add(c)
            idc = headers_raw.index(c)
            plt.plot(Select_data[:,idt],Select_data[:,idc])
            plt.title(c+' vs Time of melt')
            plt.show()
            plt.plot(Depth,Select_data[:,idc])
            plt.title(c+' vs Depth')
            plt.show()

if __name__ == '__main__':
    clean_main({'CH4':['CH4','CO']})
