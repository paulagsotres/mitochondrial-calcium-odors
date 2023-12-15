# -*- coding: utf-8 -*-
"""
Created on Tue May  9 15:58:03 2023

@author: pgomez
"""

from pathlib import Path
import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy import signal as ss
from scipy.signal import filtfilt
from scipy.signal import savgol_filter
from matplotlib.backends.backend_pdf import PdfPages
from configuration import dir_path, filename_analysis, filename_pdf, filter_window, framerate, photobleach, type_smoothing, window_length


save_path_analysis= Path (dir_path.joinpath(filename_analysis))
save_path_pdf = Path (dir_path.joinpath(filename_pdf))


#AFTER THIS, YOU SHOULDN'T HAVE TO CHANGE ANYTHING
def open_data (file_path):
    """
    Opens the signal from a matlab file that records two channels
    """
    with h5py.File (file_path) as f:
            signal = f['sig'][()][:-30, :]
    return signal #specific to open keys from a h5py

def take_all_files (dir_path):
    """
    Function to take all the files that contain the .mat file type in aa folder and extract all signal
    Returns a dictionary with the mouse number and the signal extracted for each file

    """
    # Create an empty dictionary to store the results
    signal_dict = {}
    info_experiment_list = []
    info_experiment = pd.DataFrame(columns=['mouse #', 'sensor', 'date', 'condition'])
    # Iterate over all files in the directory
    for file_path in Path(dir_path).glob("*.mat"):
        # Extract the mouse number from the file name
        mouse_num = file_path.name.split("_")[0]
        # Call the extract_behaviors_file function to get a dictionary of variables for this file
        variables_dict = open_data(str(file_path))
        # Check if the mouse number is already a key in the results dictionary
        if mouse_num in signal_dict:
            # If it is, append the variables dictionary to the list of dictionaries for that mouse
            signal_dict[mouse_num].append(variables_dict)
        else:
            # If it's not, create a new list with the variables dictionary and store it in the results dictionary
            signal_dict[mouse_num] = [variables_dict]
        #extract info from name to save it in a dataframe
        mouse_number = file_path.name.split("_")[0]
        sensor = file_path.name.split("_")[1]
        date = file_path.name.split("_")[2]
        condition = file_path.name.split("_")[3]
        # Add the data to the list
        info_experiment_list.append([mouse_number, sensor, date, condition])
        
    info_experiment = pd.DataFrame(info_experiment_list, columns=['mouse #', 'sensor', 'date', 'condition'])
    return signal_dict,info_experiment

all_signal, all_info = take_all_files (dir_path)


#--------------generate timestamp--------------

def timestamp (signal, framerate, photobleach):
    """
    This function generates two columns that we might find useful later, one with the timestamp in frames
    and the other one with the timestamp equivalent in minutes,  without counting the photobleaching period

    """
    timestamp_list = []
    for i in range(len(signal)):
        timestamp_list.append (i)
    timestamp_list = timestamp_list [photobleach:]
    timestamp_frames = np.array(timestamp_list)
    timestamp_seconds =  timestamp_frames/framerate
    return timestamp_frames, timestamp_seconds


#-------------remove period of photobleaching at the beginning
def remove_bleaching (signal, photobleach):
    """
    removes fast photobleaching period stated at the beginning of the code and excludes it from the signal processing
    and subsequent excel
    """
    isosbestic_new = signal[:, 0][photobleach:]
    rawsignal_new = signal[:, 1][photobleach:]
    return isosbestic_new, rawsignal_new


#------- smoothing of the signal

def smoothing (isosbestic, rawsignal, filter_window):
    """
    This function makes the average of several values stated by the filter window in case you want less resolution of the data
    """
    b = np.divide(np.ones((filter_window,)), filter_window)
    a = 1
    control_smooth = ss.filtfilt(b, a, isosbestic)
    signal_smooth = ss.filtfilt(b, a, rawsignal)
    return control_smooth, signal_smooth

#----------------
def savitzky_golay_filter(isosbestic, rawsignal, window_length, polyorder = 4):
    # Apply Savitzky-Golay filter to the signal
    control_smooth = savgol_filter(isosbestic, window_length, polyorder)
    signal_smooth = savgol_filter(rawsignal, window_length, polyorder)
    return control_smooth, signal_smooth


#---------------------fitting

# Fit the control channel to the signal channel using a least squares polynomial fit of degree 1.
def controlfit(isosbesticsmooth, rawsignalsmooth):
    """
    This functions adjusts the isosbestic signal to the sensor signal using a linear fit (polynomial of degree 1)
    """
    p = np.polyfit(isosbesticsmooth, rawsignalsmooth, 1)
    fitted_isosbestic = np.polyval(p, isosbesticsmooth)
    return fitted_isosbestic


#------------------calculate dff after fitting


def compute_delta_f_over_f(fitted_control_channel, signal_channel):
    """
    Computes ΔF/F by subtracting the fitted control channel from the signal channel and dividing by the fitted control channel.

    Parameters:
    control_channel (numpy array): The control channel signal.
    signal_channel (numpy array): The signal channel signal.

    Returns:
    delta_f_over_f (numpy array): The ΔF/F signal.
    """

    # Compute ΔF/F by subtracting the fitted control channel from the signal channel, and dividing by the fitted control channel.
    delta_f_over_f = (signal_channel - fitted_control_channel) / fitted_control_channel
    
    return delta_f_over_f


#--------------------calculate zscore
def z_score(signal):
    """
    Calculates the z-score of a signal by subtracting the mean and dividing by the standard deviation of
    the signal.
    """
    mean= np.mean(signal)
    std = np.std(signal)
    z_score = (signal - mean) / std

    return z_score

#----------plots the signal per mouse
def plot_processed_signal (timestamp, control_channel, fitted_control_channel, signal_channel, signal_smooth, deltaff):
    fig, axs = plt.subplots(4, 1, figsize=(10, 8), sharex = True)

    x = timestamp
    axs[0].plot(x, signal_channel, c = 'green', label = 'control')
    axs[0].set_ylabel('Intensity(mV)')
    
    axs[1].plot(x, signal_channel, c = 'green', label = 'signal')
    axs[1].plot(x, control_channel, c = 'violet', label = 'control')
    axs[1].set_ylabel('Intensity(mV)')
    
    axs[2].plot(x, signal_smooth, c = 'green', label = 'signal smoothed')
    axs[2].plot(x, fitted_control_channel, c = 'violet', label = 'Isosbestic fitted')
    axs[2].set_ylabel('Intensity (mV)')
    
    
    axs[3].plot(x, deltaff, c = 'black')
    axs[3].axhline (y=0, color='black', linestyle='-', alpha =0.1)
    axs[3].set_ylabel('%ΔF/F')
    

    axs[2].set_xlabel("Time (s)")
    plt.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    
    plt.subplots_adjust(top=0.92)
    
    return fig, axs
#
#
#
#Functions that compile all the other ones

def photometry_preprocessing(input_dict, photobleach, framerate, filter_window, type_smoothing):
    output_dict = {}
    # Loop through the input dictionary
    for key, value in input_dict.items():
        # Loop through the arrays for this key
        time_frames, time_seconds = timestamp (value[0], framerate, photobleach)
        isosbestic, raw_signal = remove_bleaching (value[0], photobleach)
        
        if type_smoothing is "moving_average":
            iso_smooth, signal_smooth = smoothing (isosbestic, raw_signal, filter_window)
        elif type_smoothing is "savitzky":
            iso_smooth, signal_smooth = savitzky_golay_filter(isosbestic, raw_signal, window_length, polyorder = 4)
        else:
            iso_smooth, signal_smooth = isosbestic, raw_signal
        fitted_iso = controlfit(iso_smooth, signal_smooth)
        deltafof = compute_delta_f_over_f(fitted_iso, signal_smooth)
        perc_dfof = deltafof*100
        zscore = z_score (deltafof)
        output_list = np.vstack((time_seconds, isosbestic, raw_signal, iso_smooth, signal_smooth, fitted_iso, perc_dfof,zscore))
        output_dict [key] = output_list.T
    return output_dict
    
all_data_processed = photometry_preprocessing(all_signal, photobleach, framerate, filter_window, type_smoothing)

def export_excel (input_dict, info, save_path):
    output_dfs = []
    output_dict = {}
    # Loop through the dictionary and create a dataframe for each key-value pair
    for key, value in input_dict.items():
        suffix = f"_{key}"
        column_names = [f"{col}{suffix}" for col in ["# seconds", "Control raw", "Signal raw", "Control smooth", "Signal smooth", "fitted iso", "%AF/F", "Z-score"]]
        df = pd.DataFrame(value, columns=[column_names])
        output_dict [key] = df
        output_dfs.append(df)
        
        
    # Concatenate all dataframes together along the columns axis
    output_df = pd.concat(output_dfs, axis=1)
    output_df.columns = [col[0] for col in output_df.columns]
    sorted_df = output_df.T.sort_index().T
    all_summary = pd.concat ([info, sorted_df], axis=1)
    pd.DataFrame(all_summary).to_csv (save_path)
    
    # Print the resulting dataframe
    return output_df, output_dict, all_summary

data_together, data_dict, summary = export_excel (all_data_processed, all_info, save_path_analysis)

def plot_signal_all(dic, df, save_path):
    pdf_pages = PdfPages(save_path)
    for key, value in dic.items():
        time = value[:, 0]
        iso = value[:, 1]
        raw_signal = value[:, 2]
        signal_smooth = value[:, 4]
        fitted_iso = value[:, 5]
        dfof = value[:, 6]
        fig, axs = plot_processed_signal(time, iso, fitted_iso, raw_signal, signal_smooth, dfof)
        if key in df.iloc[:, 0].values:
        # get the corresponding row
            row = df[df.iloc[:, 0] == key].iloc[0]
            # create the variable with mouse # + condition + date
            variable = f"{row['mouse #']}_{row['condition']}_{row['date']}"
        fig.suptitle(str(variable))
        pdf_pages.savefig(fig, bbox_inches='tight', dpi=150, metadata={'Title': str(key)})
    
    pdf_pages.close()
    return fig
        
    
plot_signal_all(all_data_processed, all_info, save_path_pdf)

