# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:25:55 2023

@author: pgomez
"""
from pathlib import Path
import pandas as pd
import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from matplotlib.backends.backend_pdf import PdfPages
from scipy.signal import find_peaks
from scipy.signal import butter
from scipy.signal import filtfilt
from process_allfiles_reportFP import data_dict, summary, all_info
from configuration import dir_path, filename_analysis, filename_behavior, framerate, window_size, bandpass, start_behavior, name_behav,distance, threshold, baseline, time_after_onset, time_before_onset, time_exposure

filename_extracted = dir_path.stem + ("_extracted.csv")
save_path_analysis= Path (dir_path.joinpath(filename_analysis))
save_path_behaviorpdf = Path (dir_path.joinpath(filename_behavior))
save_path_extracted = Path (dir_path.joinpath(filename_extracted))


def open_behavior (file_path):
    """
    Opens the behavior from a boris csv file 
    """
    behavior = pd.read_csv(file_path)
    return behavior #specific to open keys from a h5py

def take_all_files (dir_path, name_behav):
    """
    Function to take all the files that contain files with item specified in config as "name_behav" in the parent folder, and takes BORIS binary data from all behaviors
    Returns a dictionary with the mouse number and the behaviors extracted for each file

    """
    # Create an empty dictionary to store the results
    behavior_dict = {}
    # Iterate over all files in the directory
    for file_path in Path(dir_path).glob(name_behav):
        # Extract the mouse number from the file name
        mouse_num = file_path.name.split("_")[0]
        # Call the extract_behaviors_file function to get a dictionary of variables for this file
        variables_dict = open_behavior(str(file_path))
        # Check if the mouse number is already a key in the results dictionary
        if mouse_num in behavior_dict:
            # If it is, append the variables dictionary to the list of dictionaries for that mouse
            behavior_dict[mouse_num].append(variables_dict)
        else:
            # If it's not, create a new list with the variables dictionary and store it in the results dictionary
            behavior_dict[mouse_num] = [variables_dict]
       
            
    return behavior_dict

behavior_dict = take_all_files (dir_path, name_behav)

#-------------------------------
def baseline_z_score(signal, baseline_start, baseline_end):
    """
    Calculates the baseline z-score of a signal by subtracting the mean and dividing by the standard deviation of
    a designated baseline period.

    Args:
        signal (numpy array): The signal to calculate the baseline z-score of.
        baseline_start (int): The starting index of the baseline period in the signal.
        baseline_end (int): The ending index of the baseline period in the signal.

    Returns:
        numpy array: The baseline z-scored signal.
    """
    baseline_signal = signal[baseline_start:baseline_end]
    mean_baseline = np.mean(baseline_signal)
    std_baseline = np.std(baseline_signal)
    baseline_z_score = (signal - mean_baseline) / std_baseline

    return baseline_z_score
#-------------------------------


def butter_bandpass(lowcut, highcut, fs, order):
    nyq = 0.5 * fs  # Nyquist frequency
    low = lowcut/nyq  # Normalized lower cutoff frequency
    high = highcut/nyq  # Normalized upper cutoff frequency
    b, a = butter(order, [low, high], btype='band')  # Compute filter coefficients
    return b, a

def apply_butterpass (signal, lowcut, highcut, framerate, order):
    b, a = butter_bandpass(lowcut, highcut, framerate, order)
    filtered = filtfilt(b, a, signal)
    return filtered 



#------------------------------
def detect_transients_prominence(signal, window_size, threshold, excel_filename):
    """
    Function that detects calcium peak transients using a dynamic threshold calculated as Threshold * MAD during the specified window size
    Uses the peak prominence, meaning how much a peak stands out from the surrounding baseline of the signal and is defined as the vertical 
    distance between the peak and its lowest contour line, to determine whats a peak. 

    """
    # Calculate the median absolute deviation (MAD) in the moving window
    kernel = np.ones(window_size) / window_size
    
    # Calculate the median absolute deviation (MAD) in the moving window
    mad = np.zeros(signal.shape)
    for i in range(signal.shape[0]):
        start = max(0, i - window_size // 2)
        end = min(signal.shape[0], i + window_size // 2)
        mad[i] = np.median(np.abs(signal[start:end] - np.median(signal[start:end])))
    
    # Convolve the signal with the kernel to obtain the moving average
    moving_avg = np.convolve(signal, kernel, mode='same')
    # Use the MAD as a threshold to detect calcium events in the signal
    thres = threshold * mad
    peaks, properties = find_peaks(signal, prominence=thres)
    
    # # Calculate the time intervals between consecutive peaks
    # time_intervals = np.diff(peaks)
    # # Calculate the frequencies
    # frequencies = 1 / (time_intervals / len(signal))
    
    # # Create a DataFrame to store the results
    # df = pd.DataFrame({'Index': peaks, 'Value': signal[peaks], 'Frequency': frequencies})

    # # Save the DataFrame to an Excel file
    # df.to_excel(excel_filename, index=False)
    
    # Return the indices and values of the detected calcium events
    event_values = signal[peaks]
    
    return peaks, event_values

#-----------------------------
def detect_transients_by_height(data, threshold, distance):
    """
    Function that detects peaks by height from zero, and within a specific distance 

    """
    # Find peaks in the signal
    mad = np.median(np.abs(data - np.median(data)))
    thres = threshold * mad
    peaks, _ = find_peaks(data, height=thres, distance=distance)
    event_values = data [peaks]
    return peaks, event_values

#--------------------

def match_behavior_data (behavior_dict, data_dict):
    extracteddata_dict = {}
    
    for key, value in data_dict.items():
            extracteddata_dict[key] = value.iloc[:, -2:]
    
    matched_dict = {}
    
    for key in behavior_dict:
        if key in extracteddata_dict:
            value1 = behavior_dict[key]
            value2 = extracteddata_dict[key]
            matched_dict[key] = [value1, value2]
    
    matched_dict_corrected = {}
    for key in matched_dict:
        length_signal = len (matched_dict[key][1])
        difference = len (matched_dict[key][0][0])- length_signal -1
        behavior_corrected = matched_dict[key][0][0].iloc [difference:-1, :]
        matched_dict_corrected[key]= [behavior_corrected, matched_dict[key][1]]
        
    return matched_dict_corrected

dict3 = match_behavior_data (behavior_dict, data_dict)

#------------------------

    
def plot_signal_behavior_peaks(timestamp, signal, behavior, peaks, window_size, threshold, distance):
    fig, ax = plt.subplots(figsize=(10, 3), sharex = True)
    ax.plot (peaks, signal[peaks], "r+")
    ax.plot(signal, c = 'black', label='Zscore')
    legends = []
    # Get the column names starting from the second column
    for column in behavior:
        if column == "time":
            columns = behavior.columns[1:]
            break
        else:
            columns = behavior.columns

    for i, column in enumerate(columns):
        events =np.where(behavior[column]>0)[0]
        for j in range(len(events)//2):
            start_idx = events[2*j]
            color = plt.cm.get_cmap('tab10')(i)  # Get a color based on column index
            rect = plt.Rectangle((start_idx, np.min(signal)), 1, np.max(signal)-np.min(signal), color=color, alpha=0.2)
            ax.add_patch(rect)
    
        legends.append((rect, column))
    
    plt.axhline(y=0, color='black', linestyle='-', alpha =0.4)
        
    ax.xaxis.set_major_locator(MultipleLocator(1200))
    ax.set_ylabel("%AF/F")
    ax.set_xlabel("Time(s)")
    ax.legend(*zip(*legends), loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    return fig, ax

#----------------------------


def all_plots_behavior(merged_dict, save_path_behaviorpdf, window_size, threshold, distance, framerate, bandpass, baseline):
    pdf_pages = PdfPages(save_path_behaviorpdf)
    for key in merged_dict:
        behavior = merged_dict [key][0]
        timestamp = behavior.iloc [:, 0]
        signal = merged_dict [key][1]
        deltafof = signal.iloc [:, 0]
        
        if baseline is None:
            zscore =  signal.iloc [:, 1]
        else:
            zscore = baseline_z_score(deltafof, 0, baseline)

        filtered_zscore = apply_butterpass (np.array(zscore).T, bandpass['lowcut'], bandpass['highcut'], framerate, bandpass['order'])
        peaks_idx, event_values = detect_transients_prominence(filtered_zscore, window_size, threshold,save_path_behaviorpdf.parent.joinpath(f"{key}_peaks.csv") )
        fig, ax = plot_signal_behavior_peaks(timestamp, deltafof, behavior, peaks_idx, window_size, threshold, distance)
        ax.set_xticklabels(ax.get_xticks() // framerate)
        fig.suptitle(str(key))
        pdf_pages.savefig(fig, bbox_inches='tight', dpi=150, metadata={'Title': str(key)})
    
    pdf_pages.close()
    plt.show()
        
        
all_plots_behavior(dict3, save_path_behaviorpdf, window_size, threshold, distance, framerate, bandpass, baseline)

#------------------
def extract_fiber_photometry_signal(signal, time_before_onset, time_after_onset, behavior1, behavior2, framerate, save_path_extracted):
    # Find the onset and offset of the first and last behavioral events

    onset1 = np.where(behavior1 == 1)[0][0]

    onset2 = np.where(behavior2 == 1)[0][0]



    # Extract the signal between 15s before the onset and 60s after the offset
    start1 = max(0, onset1 - time_before_onset*framerate)
    end1 = min(len(signal), onset1 +  time_after_onset*framerate)
    extracted_signal1 = signal[start1:end1]
    
    start2 = max(0, onset2 - time_before_onset*framerate)
    end2 = min(len(signal), onset2 +  time_after_onset*framerate)
    extracted_signal2 = signal[start2:end2]
    
    
    
    startbaseline = time_before_onset*framerate
   # Calculate the baseline zscore of the extracted signals 
    baseline_signal1 = signal[onset1-startbaseline:onset1]
    mean_baseline1 = np.mean(baseline_signal1)
    std_baseline1 = np.std(baseline_signal1)
    extracted_zscored1 = (extracted_signal1 - mean_baseline1)
    
    
    baseline_signal2 = signal[onset2-startbaseline:onset2]
    mean_baseline2 = np.mean(baseline_signal2)
    std_baseline2 = np.std(baseline_signal2)
    extracted_zscored2 = (extracted_signal2 - mean_baseline2)
    
    
    return extracted_zscored1, extracted_zscored2

#-----------------
def extract_fiber_photometry_signals(signal, time_before_onset, time_after_onset, behavior, framerate, df):
    total_event = 0
    extracted_signals = []

    for i in range(len(df)):
        if i == 0:
            continue
        else:
            # Check if the 'behavior' column matches the specified value
            if df.iloc[i, 0] == behavior:
                if df.iloc[i - 1, 0] != behavior:
                    start = max(0, i - time_before_onset * framerate)
                    end = min(len(signal), i + time_after_onset * framerate)
                    extracted_signal = signal[start:end]

                    # Calculate the baseline z-score of the extracted signal
                    start_baseline = max(0, i - time_before_onset * framerate)
                    baseline_signal = signal[start_baseline:i]
                    mean_baseline = np.mean(baseline_signal)
                    std_baseline = np.std(baseline_signal)
                    extracted_zscored = (extracted_signal - mean_baseline) / std_baseline

                    # Append the extracted z-scored signal to the list
                    extracted_signals.append(extracted_zscored)

    # Create a dataframe with columns for each extracted signal
    df_result = pd.DataFrame({f'Signal_{i+1}': signal for i, signal in enumerate(extracted_signals)})

    return df_result


def save_signals_to_excel(file_path, signal_behavior_dict, time_before_onset, time_after_onset, framerate):
    with pd.ExcelWriter(file_path, engine='xlsxwriter') as writer:
        for key, value in signal_behavior_dict.items():
            signal_df = value['signal']
            behavior_df = value['behavior']

            # Extract signals using the provided function
            result_df = extract_fiber_photometry_signals(signal_df, time_before_onset, time_after_onset, behavior_df, framerate)

            # Save the result_df to the Excel file with a sheet named after the key
            result_df.to_excel(writer, sheet_name=key, index=False)


def AUC_extracted_signal (extracted_signal, time_exposure, time_before_onset):
    start_baseline = 0
    stop_baseline = time_before_onset*framerate

    
    start_exposure = time_before_onset*framerate
    stop_exposure = start_exposure + (time_exposure*framerate)

    
    # Calculate area under the curve for baseline and exposure periods
    baseline_area = np.trapz(extracted_signal[start_baseline:stop_baseline])
    exposure_area = np.trapz(extracted_signal[start_exposure:stop_exposure])
    
    return baseline_area, exposure_area
#----------------
def plot_extracted_signals (extracted_signals, time_exposure, time_before_onset, framerate):
    fig, axs = plt.subplots(nrows=len(extracted_signals.columns), figsize=(10, 8), sharex=True)

    # Loop over each column in the extracted_signals DataFrame
    for i, column_name in enumerate(extracted_signals.columns):
        ax = axs[i]
        ax.plot(extracted_signals[column_name], c='black', label='Zscore')
    
        start_baseline = 0
        stop_baseline = time_before_onset * framerate
        rect_baseline = plt.Rectangle((start_baseline, np.min(extracted_signals[column_name])),
                                      stop_baseline - start_baseline,
                                      np.max(extracted_signals[column_name]) - np.min(extracted_signals[column_name]),
                                      color="grey", alpha=0.1)
        ax.add_patch(rect_baseline)
    
        start_exposure = time_before_onset * framerate
        stop_exposure = start_exposure + (time_exposure * framerate)
        rect_exposure = plt.Rectangle((start_exposure, np.min(extracted_signals[column_name])),
                                       stop_exposure - start_exposure,
                                       np.max(extracted_signals[column_name]) - np.min(extracted_signals[column_name]),
                                       color="green", alpha=0.1)
        ax.add_patch(rect_exposure)
    
        ax.set_ylabel("baseline Z-Score")
        ax.set_xlabel("frames")
        ax.set_title("Extracted Signal {}".format(column_name))
        ax.legend(loc='best')
    
    # Adjust the subplot layout
    plt.tight_layout()
    plt.show()

def extract_all(merged_dict, time_before_onset,time_exposure, time_after_onset, framerate, save_path_extracted):
    for key in merged_dict:
        behavior = merged_dict [key][0]
        timestamp = behavior.iloc [:, 0]
        signal = merged_dict [key][1]
        deltafof = signal.iloc [:, 0]
        
        if baseline is None:
            zscore =  signal.iloc [:, 1]
        else:
            zscore = baseline_z_score(deltafof, 0, baseline)

        filtered_zscore = apply_butterpass (np.array(zscore).T, bandpass['lowcut'], bandpass['highcut'], framerate, bandpass['order'])
        
        stress = pd.DataFrame()
        saline = pd.DataFrame()
        
        # Iterate through the columns
        for column in behavior.columns:
            if 'stress' in column:
                stress[column] = behavior[column]
            elif 'saline' in column:
                saline[column] = behavior[column]

        extracted_zscore_saline, extracted_zscore_stress = extract_fiber_photometry_signal(zscore, time_before_onset, time_after_onset, saline, stress, framerate, save_path_extracted)
        extracted_zscore_saline.reset_index(drop=True, inplace=True)
        extracted_zscore_stress.reset_index(drop=True, inplace=True)
        extracted_signal = pd.concat ([pd.DataFrame (extracted_zscore_saline), pd.DataFrame (extracted_zscore_stress)], axis = 1)
        extracted_signal.columns = ['saline', 'SBT']
        
        save_path_extracted = Path (dir_path.joinpath(f"{key}_extracted.csv"))
        extracted_signal.to_csv(save_path_extracted, index = False)
        
        baseline_areasaline, exposure_areasaline = AUC_extracted_signal (extracted_zscore_saline, time_exposure, time_before_onset)

        baseline_areastress, exposure_areastress = AUC_extracted_signal (extracted_zscore_stress, time_exposure, time_before_onset)

        print (f"{key}_saline :{baseline_areasaline,  exposure_areasaline}")
        print (f"{key}_stress :{baseline_areastress,exposure_areastress}")
        
        plot_extracted_signals (extracted_signal, time_exposure, time_before_onset, framerate)
        
        
extract_all(dict3, time_before_onset,time_exposure, time_after_onset, framerate, save_path_extracted)
