import numpy as np
import os
from math import floor
from tkinter import filedialog
import tkinter as tk
import matplotlib.pyplot as plt

def merge_main(wd=None):
    '''

    :param wd: Working Directory. If None, one is chosen with gui
    :return: nothing but save file with merged HK and AvC files
    '''

    # -------------------- Get working directory and prepare parameters -----------------------------------
    if wd is None:
        application_window = tk.Tk()
        wd = filedialog.askdirectory(title='Select day directory', parent=application_window)
        application_window.destroy()
    if not wd: raise Exception('No directory selected...')

    DateP=os.path.basename(wd)
    renames_Avc = {'CH4(ppb)':'CH4', 'CO(ppb)':'CO', 'H2O(%)':'H2O'} # {'Old Name':'New Name'}
    dump_Avc = []
    renames_HK = {'inP_from_Ctrl':'P_CFA'}
    dump_HK = ['Avg_sym','Bsln1','Bsln0']

    # ------------------- MERGE AVC FILES -----------------------------------------------------------------
    print('Merging AvC files')
    store = []
    for f in os.listdir(wd): # Browse all files in wd
        # Only if end with AvC.dat
        if f.endswith('AvC.dat'):
            store.append(f)
    store.sort()
    for f in store:
        print(f)
        if 'Avc_raw' in locals(): # Avc_raw already filled with data -> Concatenate
            Avc_raw = np.vstack((Avc_raw,np.loadtxt(os.path.join(wd,f),skiprows=1)))
        else: #Avc_raw to be created
            Avc_raw = np.loadtxt(os.path.join(wd,f),skiprows=1)
            # Open file and read first line : Header
            f = open(os.path.join(wd,f),mode='r')
            headers_Avc = f.readline().rstrip().split('\t')
            f.close()

    # Remove points which are not from the good day (i.e. not the good floor)
    idt = headers_Avc.index('Time')
    T_floor = int(floor(Avc_raw[int(len(Avc_raw)/2),idt]))
    DataAvc = Avc_raw[np.where(np.floor(Avc_raw[:,idt])==T_floor)] # Keep only points where floor(time) == T_floor
    DataAvc[:,idt] = (DataAvc[:,idt] - T_floor)*24. # Scale in hours


    # -------------------- MERGE HK FILES -----------------------------------------------------------------------

    # Do exactly the same as the above with AvC
    print('Merging HK files')
    store = []
    for f in os.listdir(wd):
        if f.endswith('HK.dat'):
            store.append(f)
    store.sort()
    for f in store:
        print(f)
        if 'HK_raw' in locals():
            try : HK_raw = np.vstack((HK_raw,np.loadtxt(os.path.join(wd,f),skiprows=1)))
            except: print(f)
        else:
            HK_raw = np.loadtxt(os.path.join(wd,f),skiprows=1)
            f = open(os.path.join(wd,f),mode='r')
            headers_HK = f.readline().rstrip().split('\t')
            f.close()

    # Remove points which are not from the good day (i.e. not the good floor)
    idthk = headers_HK.index('Time')
    T_floor = int(floor(HK_raw[int(len(HK_raw)/2),idthk]))
    DataHK = HK_raw[np.where(np.floor(HK_raw[:,idthk])==T_floor)]
    DataHK[:,idthk] = (DataHK[:,idthk] - T_floor)*24.

    # Interpolate HK data on Instrument time scale
    DataHK_rs = np.zeros((len(DataAvc[:,0]),len(DataHK[0,:])))
    for i in range(len(DataHK[0,:])):
        DataHK_rs[:,i]=np.interp(DataAvc[:,idt],DataHK[:,idthk],DataHK[:,i])



    # ------------------------ Add Melt rate -------------------------------------- #placer ce bloc en commentaire si pas de fichier RUN
#    Date = os.path.basename(wd)
#    Year = Date[0:4]
#    Month = Date[4:6]
#    Day = Date[6:8]
#    file_melt = os.path.join(wd,'RUN-'+Year+'-'+Month+'-'+Day+'.txt')
#    time_melt, depth, melt = np.loadtxt(file_melt,unpack=True)
#    melt = rolling_median(melt, window=20)
#    melt= np.interp(DataAvc[:,idt], time_melt, melt)

    # ---------------- SAVE DATA ---------------------------------------------------------
    outfile = os.path.join(wd,'SARA_MergedFile_'+DateP+'.txt')

    if os.path.exists(outfile): os.remove(outfile) # If the file already exists, delete and rewrite.

    # Build header and output array keeping only column of interest and with right name
    header = None
    Array_out  = None
    for i,head in enumerate(headers_Avc):
        if head in dump_Avc: continue
        if head in renames_Avc: head =renames_Avc[head]
        if not header is None:
            header = header+ '\t' +head
        else:
            header = head
        if not Array_out is None:
            Array_out = np.column_stack((Array_out,DataAvc[:,i]))
        else:
            Array_out = DataAvc[:,i]
    for i,head in enumerate(headers_HK):
        if head in dump_HK: continue
        if head in renames_HK: head =renames_HK[head]
        if not header is None:
            header = header + '\t' + head
        else:
            header = head
        if not Array_out is None:
            Array_out = np.column_stack((Array_out,DataHK_rs[:,i]))
        else:
            Array_out = DataHK_rs[:,i]
#    header += '\tMeltRate' #placer ces deux lignes en commentaire si pas de fichier RUN
#    Array_out = np.column_stack((Array_out, melt))

    # Save data using numpy.savetxt
    np.savetxt(outfile,Array_out,delimiter='\t',header=header,fmt='%.6e')

if __name__ == '__main__':
    merge_main()
