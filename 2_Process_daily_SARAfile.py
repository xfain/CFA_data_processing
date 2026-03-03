import numpy as np
import os
from openpyxl import load_workbook
import xlrd
import matplotlib.pyplot as plt
import matplotlib as mpl
from tkinter import filedialog
import tkinter as tk
import csv


def process_main(wd=None):
    '''

    :param wd: Working Directory. If None, one is chosen with gui
    :return: nothing but save selected data
    '''
    # Get the data file name and Date

    if wd is None:
        application_window = tk.Tk()
        wd = filedialog.askdirectory(title='Select day directory', parent=application_window)
        application_window.destroy()
    if not wd: raise Exception('No directory selected...')

    Date = os.path.basename(wd)
    Year = Date[0:4]
    Month = Date[4:6]
    Day = Date[6:8]
    file = os.path.join(wd,'SARA_MergedFile_'+Date+'.txt')

    # Load Gas and Event Spreadsheets
    Daily_cfa_Eventlog_File = os.path.join(wd,'CFA_Event_Log_'+Day+'_'+Month+'_'+Year+'.xls')
    Daily_cfa_Gaslog_File = os.path.join(wd,'CFA_GasLog_'+Date+'.xls')
    Gas_log = xlrd.open_workbook(Daily_cfa_Gaslog_File)

    # Delays
    seam_internal = 0./3600
    seam_full = 0./3600

    # Get Time shift of the day between Melter and Instrument (and convert in hours)
    sheet = Gas_log.sheet_by_name('Timeshift')
    T_shift = float(sheet.cell(1,1).value)/3600.

    # Load data and read header
    Raw_data = np.loadtxt(file)
    f=open(file)
    headers_raw = f.readline().replace('# ','').rstrip().split('\t')
    f.close()

    # Get Index of interest using header
    idt = headers_raw.index('Time')
    idch4 = headers_raw.index('CH4')
    idco = headers_raw.index('CO')

    # Put instrument time apart
    Time = Raw_data[:,idt]

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

    # Get Instrument time of switches between loops and sample
    t_start = np.array([])
    t_end = np.array([])
    sheet = Gas_log.sheet_by_name('Log')
    for irow in range(1,sheet.nrows):
        if int(sheet.row(irow)[7].value)==1:
            row = sheet.row(irow)
            next_row = sheet.row(irow + 1)
            if int(sheet.row(irow-1)[4].value)==1:
                t_start = np.append(t_start, float(row[0].value) + float(row[1].value)/60. + \
                                    float(row[2].value) / 3600. + seam_internal)
                t_end = np.append(t_end, float(next_row[0].value) + float(next_row[1].value)/60. + \
                                    float(next_row[2].value) / 3600. )
            else:
                t_start = np.append(t_start, float(row[0].value) + float(row[1].value) / 60. + \
                                    float(row[2].value) / 3600. + seam_full)
                t_end = np.append(t_end, float(next_row[0].value) + float(next_row[1].value) / 60. + \
                                  float(next_row[2].value) / 3600.)

    # Keep data only where Time instrument >tstart and <tend
    whr = np.zeros(len(Time),dtype=bool)
    for i,(tst,tnd) in enumerate(zip(t_start,t_end)):
        whr = np.logical_or(np.logical_and(Time>=tst,Time<=tnd),whr)
    whr = np.logical_or(whr,np.append(np.diff(Time),1)<=0)

    # Select only good data in Raw_data
    Out_array = Raw_data[np.where(whr),:][0]

    # --------- !!!! Change time in Time Melter from now on !!!! -------------
    Out_array[:,idt] = Out_array[:,idt] - T_shift

    # Save data using numpy.savetxt
    outfile = os.path.join(wd,'SARA_cleaning_' + Date + '.txt')
    header = None
    for h in headers_raw:
        if not header is None:
            header = header + '\t' + h
        else:
            header = h

    np.savetxt(outfile,Out_array,delimiter='\t',header=header,fmt='%.6e')

    # Plot all that shit
    fig = plt.figure(facecolor='#ffffff')
    ax1 = plt.subplot(211)
    # ax1.plot(Raw_data[:,idt] - T_shift,Raw_data[:,idch4],c='blue')
    ax1.plot(Raw_data[:,idt] - T_shift ,Raw_data[:,idch4],c='k')
    ax1.plot(Out_array[:,idt],Out_array[:,idch4],c='#2C98CA')

    for t in t_break:
        ax1.axvline(x=t, ymin=0, ymax = 1, linewidth=1, color='green')
    for t in t_change:
        ax1.axvline(x=t, ymin=0, ymax = 1, linewidth=1, color='red')
    ax1.set_ylabel('CH4 (ppb)')
    ax1.set_xlim(0,24)

    ax2 = plt.subplot(212,sharex=ax1)
    #ax2.set_axis_bgcolor('#ffffff')
    # ax2.plot(Raw_data[:,idt] - T_shift,Raw_data[:,idco],c='blue')
    ax2.plot(Out_array[:,idt],Out_array[:,idco],c='#cd950c')
    for t in t_break:
        ax2.axvline(x=t, ymin=0, ymax = 1, linewidth=1, color='green')
    for t in t_change:
        ax2.axvline(x=t, ymin=0, ymax = 1, linewidth=1, color='red')
    ax2.set_ylabel('CO (ppb)')
    ax2.set_xlabel('Temps (h)')
    ax2.set_xlim(0,24)

    plt.show()

if __name__ == '__main__':
    process_main()
