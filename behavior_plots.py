# -*- coding: utf-8 -*-
"""
Created on Sun Aug  6 18:42:56 2023

@author: pgomez
"""


from pathlib import Path
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.stats import sem
from matplotlib.ticker import MultipleLocator
import csv
import math

root_directory = Path(r"E:\Experiments\2023\20231201_Reanalyzed partner social STS\OBS")
summary_all = r"E:\Experiments\2023\20231201_Reanalyzed partner social STS\Summary all.csv"

def open_behavior(file_path):
    open_test = pd.read_csv(file_path)
    first_row = open_test.columns[0]
    if first_row == "Behaviour":
        annotations = pd.read_csv(file_path).iloc[:, 0:2]
        print ("Correct delimiter")
    else:
        print ("Bad delimiter")
        sniffer = csv.Sniffer()
        detected_delimiter = sniffer.sniff(first_row).delimiter
        annotations_temporary = pd.read_csv(file_path, delimiter=detected_delimiter)
        print (annotations_temporary)

        if len(open_test.columns) > 1:
            decimal_separator = sniffer.sniff(annotations_temporary.iloc[2, 1]).decimal
            annotations = pd.read_csv(file_path, delimiter=detected_delimiter, decimal=decimal_separator).iloc[:, 0:2]
        else:
            annotations = pd.read_csv(file_path, delimiter=";", decimal=".").iloc[:, 0:2]
            print("Number of columns is not sufficient to infer the decimal separator.")
          

    return annotations



#EDIT ANNOTATIONS 
def replace_value_in_dataframe(df, column_name, old_value, new_value):
    if old_value in df[column_name].values:
        df[column_name] = df[column_name].replace(old_value, new_value)
    return df

def cut_dataframe_to_activity_window(df):
    """
    Cut the annotations from no_activity to no_activity so it skips the END_video mark

    """
    no_activity_indices = df[df['Behaviour'] == 'no_activity'].index
    end_video = df[df['Behaviour'] == 'VIDEO_END'].index[0]
    # Find the first occurrence of a behavior that isn't 'no_activity'
    if len(no_activity_indices) == 1:
        print ("only one")
        first_non_no_activity = df[df['Behaviour'] != 'no_activity'].iloc[0][' start_time(ms)']
        if df.loc[end_video-1, "Behaviour"] == "no_activity":
            last_non_no_activity = df.loc[end_video-1][ " start_time(ms)"]
            new_df = df[(df[' start_time(ms)'] >= first_non_no_activity) & (df[' start_time(ms)'] <= last_non_no_activity)]
        else:
            end_video = df[df['Behaviour'] == 'VIDEO_END'].iloc[0][' start_time(ms)']
            new_df = df[(df[' start_time(ms)'] >= first_non_no_activity) & (df[' start_time(ms)'] <= end_video)]
        
    elif len(no_activity_indices) > 1:
        print ("several")
        first_non_no_activity =df[df['Behaviour'] != 'no_activity'].iloc[0][' start_time(ms)']
        if df.loc[end_video-1, "Behaviour"] == "no_activity":
            last_non_no_activity = df.loc[end_video-1][ " start_time(ms)"]
        else:
            last_non_no_activity = df[df['Behaviour'] == 'VIDEO_END'].iloc[0][' start_time(ms)']
            
        new_df = df[(df[' start_time(ms)'] >= first_non_no_activity) & (df[' start_time(ms)'] <= last_non_no_activity)]
    else: 
        new_df = df.loc[:end_video, :]
    

    return new_df

# Function to round up numbers to one decimal place
def round_up_to_full_seconds(df):
    """
    We round up the values to full seconds by first rounding up to decimals and then to full numbers

    """
    df = df.copy()

    df[' start_time(ms)'] = df[' start_time(ms)'].apply(lambda x: int(x) + 1 if (x * 10) % 10 >= 5 else int(x))

    behavior_column = df['Behaviour']
    time_column = df.iloc[:, 1]
    
    
    first_time = time_column.iloc[0]


    first_time = time_column.iloc[0]# Find the last occurrence of 'no_activity'
    adjust_time = time_column -first_time
        
    result_df = pd.concat([behavior_column, adjust_time], axis=1)
    result_df.reset_index(drop=True, inplace=True)
    return result_df


def behavior_for_each_second(df):
    """
    We establish a cell per second, so we can have a timeline per second of the behavior

    """
        # Initialize an empty DataFrame to store the result
    if 'no_activity' in df['Behaviour'].values:
    # Replace 'no_activity' with 'sit'
        df['Behaviour'] = df['Behaviour'].replace('no_activity', 'sit')
        
    max_time = df[' start_time(ms)'].max()
    new_df = pd.DataFrame({'Time': range(max_time), 'event ID': [None] * max_time})
    
    for index, row in new_df.iterrows():
        time_new = new_df.loc[index, "Time"]
        time_og_annotations = df[df[' start_time(ms)'] == time_new]
        if len(time_og_annotations) > 1:
            if (time_og_annotations[" start_time(ms)"] == "sniff ano-gen").any():
                new_behavior = "sniff ano-gen"
            else:
                new_behavior = time_og_annotations.iloc[-1]["Behaviour"]
        elif len(time_og_annotations) == 0:
            continue
        else:
            new_behavior = time_og_annotations["Behaviour"].values[0]  # Extract the first value as a scalar
        new_df.at[index, "event ID"] = new_behavior
        
    new_df['event ID'] = new_df['event ID'].ffill()
    
    last_time = int(new_df['Time'].max())
    if last_time < 300:
        # Create a DataFrame to fill the missing values with 'sit'
        missing_df = pd.DataFrame({'Time': range(last_time + 1, 301), 'event ID': 'sit'})

        # Concatenate the original DataFrame and the missing values DataFrame
        result_df = pd.concat([new_df, missing_df], ignore_index=True)
        result_df = result_df[result_df['Time'] <= 300]
    elif last_time >= 300:
        result_df = new_df[new_df['Time'] <=300]

    return result_df

#DETECT BEHAVIORS
def detect_behaviors(dataframe, behaviors):
    result_dict = {}
    for column in dataframe.columns:
        column_dict = {}
        for behavior in behaviors:
            behavior_series = (dataframe[column] == behavior)
            column_dict[behavior] = behavior_series
        result_dict[column] = column_dict
    return result_dict

#CALCULATE BY TIME
def count_behavioral_cum_events(df, behavior):
    events = []
    total_event = 0
    for i in range(len(df)):
        if i == 0:
            if df.iloc[i, 0] == behavior:
                total_event += 1
            else:
                total_event = 0
        else:
            # Check if the 'behavior' column matches the specified value
            if df.iloc[i, 0] == behavior:
                if df.iloc[i - 1, 0] != behavior:
                    total_event += 1
                    
        events.append(total_event)
    cum_events = pd.DataFrame (events, columns =["cumulative events"])
    
    median = float(cum_events.max()/2)
    if median % 2 == 0:  # If the number is even
    # Find all values in the DataFrame that correspond to the even number
        even_number_values = cum_events[cum_events["cumulative events"] == int(median)]
    # Find the index of the last row that has the even number
        median_time = even_number_values.index[-1]
    else:  # If the number is odd
        # Find all values in the DataFrame that correspond to the odd number + 1
        odd_number_values = cum_events[cum_events["cumulative events"] == int(median + 1)]
        median_time = odd_number_values.index[-1]

    return cum_events, median_time



def calculate_freq_and_events(dataframe):
    behaviors = ["sniff ano-gen", "sniff head-torso", "allogroom", "groom", "rear", "dig", "walk", "sit", "fight"]
    data = detect_behaviors(dataframe, behaviors)
    freq_individual_result = {}
    event_cum_norm_result = {}
    events_result ={}
    events_individual_result = {}
    for mouse, values in data.items():
        for event, values_dict in values.items():
            events_individual_dict ={}
            freq_individual_dict ={}
            for row in range(len(values_dict)):
                events_individual = [int(b) for b in [data[mouse][event][row] for mouse in data]]
                events_individual_dict[row] = events_individual
                freq_individual_dict [row] =[i / len(dataframe.columns) for i in events_individual]
             
            events_individual_result[event] = pd.DataFrame.from_dict(events_individual_dict, orient='index', columns=data.keys())
            freq_individual_result[event] = pd.DataFrame.from_dict(freq_individual_dict, orient='index', columns=data.keys())
  
    for key, data_all in events_individual_result.items():
        cumulative_df = data_all.cumsum()
        event_cum_norm_result[key] = cumulative_df
            
    return events_individual_result, freq_individual_result, event_cum_norm_result

def calculate_bout_duration_behavior (input_df, behavior):
    contact_number = 0
    duration = 0
    bout_durations = []

    for index, row in input_df.iterrows():
        if row['Behaviour'] == behavior:
            if index + 1 < len(input_df):
                next_row = input_df.loc[index + 1]
                duration = next_row [1]- row[1]
                contact_number += 1
                bout_durations.append((contact_number, duration))
                
    bout_df = pd.DataFrame(bout_durations, columns=["bout number", "duration"])
    return bout_df

def calculate_latency_to_bouts (df, behavior):
    contact_number = 0
    duration = 0
    bout_latency = []
    first_row = df.loc[1]
    for index, row in df.iterrows():
        if row['Behaviour'] == behavior:
            if index + 1 < len(df):
                latency = row[1]-first_row[1]
                contact_number += 1
                bout_latency.append((contact_number, latency))
    bout_latency_df = pd.DataFrame(bout_latency, columns=["bout number", "latency to start"])   
    return bout_latency_df


def calculate_latency_between_bouts(df, behavior):
    bout_latency = []
    first_bout_latency = None
    prev_row = None
    found_behavior = False
    first_row = df.loc[1]
    for index, row in df.iterrows():
        if row['Behaviour'] == behavior:
            if not found_behavior:
                first_bout_latency = row[1]-first_row[1]
                found_behavior = True
            if prev_row is not None:
                latency = row[' start_time(ms)'] - prev_row[' start_time(ms)']
                bout_latency.append(latency)
            prev_row = row
        
    bout_latency_all = [first_bout_latency] + bout_latency 
    bout_latency_df = pd.DataFrame(bout_latency_all, columns=["interval between bouts"])
    return bout_latency_df

    

def calculate_half_time(df, behavior):
    latency_to_bouts = calculate_latency_to_bouts (df, behavior)
    duration_bouts= calculate_bout_duration_behavior(df, behavior)
    
    both = pd.concat ([duration_bouts, latency_to_bouts['latency to start']], axis=1)
    
    total_time = duration_bouts['duration'].sum()
    
    half_time = total_time / 2
    
    moment = both[both["duration"].cumsum() >= half_time]['latency to start'].iloc[0]
    
    index = both.loc[both['latency to start'] == moment].index[0]
    
    sum_until_half = both["duration"].iloc[:index+1].sum()
    
    real_moment = moment - (sum_until_half-half_time)
    return real_moment


def calculate_number_duration_totalevents (df, behavior):
    duration_bouts= calculate_bout_duration_behavior(df, behavior)
    
    number_events = duration_bouts.shape[0]
    total_time = duration_bouts["duration"].sum()
    
    return number_events, total_time


    
#PLOT ALL
def plot_frequency(dataframe, title_suffix):
    freq_individual = calculate_freq_and_events(dataframe)[1]
    freq_dict ={}
    num_time_points = len(freq_individual[list(freq_individual.keys())[0]])
    for key, value in freq_individual.items():
        sum_df = value.sum(axis=1)
        freq_dict[key] = sum_df

    num_plots = len(freq_dict)
    fig, axs = plt.subplots(num_plots, 1, figsize=(10, num_plots * 2), sharex=True, sharey=True)
    fig.tight_layout(pad=2.0)
    plt.subplots_adjust(hspace=0.4)

    for i, (behavior, data) in enumerate(freq_dict.items()):
        axs[i].imshow(data.values.reshape(1, -1), cmap='inferno', vmin=0, vmax=1, aspect='auto', extent=[0, num_time_points, 0, 1])
        axs[i].set_yticks([0.5])
        axs[i].set_yticklabels([behavior])
        axs[i].set_xlabel('Time (s)')
        axs[i].set_title(f'{behavior} Timeline')

    plt.savefig(root_directory.joinpath(f'{title_suffix}.svg'), format='svg')
    plt.show()
    plt.close()  # Close the plot



def plot_behaviors(dataframe, title_suffix, behavior_colors):
    behavior_dict = detect_behaviors(dataframe, behavior_colors.keys())
    
    fig, ax = plt.subplots(figsize=(10, 3))
    y_ticks = []
    y_tick_labels = []
    
    for i, (key, value) in enumerate(behavior_dict.items()):
        y_ticks.append(i + 0.5)
        y_tick_labels.append(key)
        
        for behavior, values in value.items():
            events = np.nonzero(np.array(values))[0]
            
            for event in events:
                color = behavior_colors.get(behavior, "white")
                rect = plt.Rectangle((event, i), 1, 1, color=color)
                ax.add_patch(rect)
    
    ax.set_yticks(y_ticks)
    ax.set_yticklabels(y_tick_labels)
    ax.set_xlabel("Time (seconds)")
    ax.set_ylabel("Mouse")
    ax.set_ylim(0, len(behavior_dict))
    ax.set_xticks(np.arange(0, 301, 30))
    fig.suptitle(title_suffix, fontsize=14)
    plt.gca().invert_yaxis()
    plt.savefig(root_directory.joinpath(f'{title_suffix}.svg'), format='svg')
    plt.show()
    plt.close()  # Close the plot

def plot_frequency_graph(dataframe, behaviors, title_suffix):
    freq_individual = calculate_freq_and_events(dataframe)[1]
    freq_dict ={}
    for key,value in freq_individual.items():
        freq_total= value.sum(axis=1)
        freq_dict[key]=freq_total
    num_plots = len(freq_dict.keys())
    fig, axs = plt.subplots(num_plots, 1, figsize=(10, num_plots*2), sharex=True)
    fig.tight_layout(pad=2.0)
    plt.subplots_adjust(hspace=0.4)

    for i, behavior in enumerate(behaviors):
        if behavior in freq_dict:
            data = freq_dict[behavior]
            time_frames = np.arange(0, len(data), 1)
            axs[i].fill_between(time_frames, 0, data, color='lightblue')
            axs[i].plot(time_frames, data, color='lightblue', linewidth=1)
            axs[i].set_xlabel('Time (frames)')
            axs[i].set_ylabel('Frequency')
            axs[i].set_title(f'{behavior} Frequency')
            axs[i].set_xlim(0, len(data) - 1)
            axs[i].set_ylim(0, 1)  # Set y-axis limits to [0, 1]
            axs[i].set_yticks(np.arange(0, 0.61, 0.2)) 
        else:
            # Plot an empty graph with no data
            axs[i].set_xlabel('Time (frames)')
            axs[i].set_ylabel('Frequency')
            axs[i].set_title(f'{behavior} Frequency (Empty)')
            axs[i].set_xlim(0, 1)  # Set minimal x-axis range
            axs[i].set_ylim(0, 0.6)  # Set y-axis limits to [0, 0.6] with ticks of 0.2
            axs[i].set_yticks(np.arange(0, 0.61, 0.2))
    
    plt.savefig(root_directory.joinpath(f'{title_suffix}.svg'), format='svg')
    plt.show()
    plt.close()  # Close the plot


#CREATE NEW DICTIONARIES
def process_data(root_directory, summary_all):
    all_data = {}
    batch_dict = {}
    
    open_test = pd.read_csv(summary_all)
    first_row = open_test.columns[0]
    if first_row == "Batch":
        summary_df = pd.read_csv(summary_all)
        print ("Correct delimiter")
    else:
        sniffer = csv.Sniffer()
        detected_delimiter = sniffer.sniff(first_row).delimiter
        summary_df = pd.read_csv(summary_all, delimiter=detected_delimiter, decimal=",")
    
    for batch_folder in os.listdir(root_directory):
        if batch_folder.startswith("Batch") and os.path.isdir(os.path.join(root_directory, batch_folder)):
        # Check if it starts with "Batch" and is a directory
            batch_dict[batch_folder] = os.path.join(root_directory, batch_folder)

        filtered_data = {}    
        # Iterate through the keys in the dictionary
        for batch, directory in batch_dict.items():
            # Extract the batch number from the key (e.g., 'Batch 1' => 1)
            batch_number = int(batch.split()[-1])
    
            # Filter the dataframe based on the 'Batch' column
            filtered_df = summary_df[summary_df['Batch'] == batch_number]
    
            # Add the filtered dataframe to the results dictionary
            filtered_data[batch] = {
                'Directory': directory,
                'DataFrame': filtered_df}
    return filtered_data

    
def find_and_store_file(dictionary_all):
    all_results = {}
    for batch_name, batch_info in dictionary_all.items():
        path_dict ={}
        directory_list_by_batch = []
        directory_path = batch_info.get('Directory')
        summary_animals = batch_info.get('DataFrame')
        for item in os.listdir(directory_path):
            item_path = os.path.join(directory_path, item)
            if os.path.isdir(item_path):  # Check if the item is a directory
                directory_list_by_batch.append(item_path)
        for path in directory_list_by_batch:
            path_elements = os.path.splitext(os.path.basename(path))[0].split("_")

            for index, row in summary_animals.iterrows():
                mouse_pair = row['mouse_pair']
                condition = row['condition']
                genotype = row['genotype']
        
                if mouse_pair in path_elements:
                    sniff_file_path = os.path.join(path, "annotations.csv")
                    if os.path.exists(sniff_file_path):  # Check if the file exists
                        if condition not in path_dict:
                            path_dict[condition] = {}
                        if genotype not in path_dict[condition]:
                            path_dict[condition][genotype] = {}
                        if mouse_pair not in path_dict[condition][genotype]:
                            path_dict[condition][genotype][f"{batch_name}_{mouse_pair}"] = {}

                            path_dict[condition][genotype][f"{batch_name}_{mouse_pair}"] = sniff_file_path
                        else: 
                            print ("ERROR")
                        
        for condition, genotypes in path_dict.items():
            if condition not in all_results:
                all_results[condition] = {}
            for genotype, mouse_pair_paths in genotypes.items():
                if genotype not in all_results[condition]:
                    all_results[condition][genotype] = {}
                all_results[condition][genotype].update(mouse_pair_paths)
        
    return all_results

def concat_all_processed_annotations_events (all_results): 
    all_batch = {}
    for condition, genotypes in all_results.items():
        concatenated_data = {}  # Dictionary to store data for each genotype in the condition
        for genotype, mouse_pairs in genotypes.items():
            concatenated_df = None  # Initialize the DataFrame to accumulate data for this genotype
            for mouse_name, path in mouse_pairs.items():
                annotations = open_behavior(path)
                print (path)
                print (annotations)
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "anogenital sniff", "sniff ano-gen")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "sit/still", "sit")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "rearing", "rear")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "self-groom", "groom")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "escape", "walk")
                annotations_active = cut_dataframe_to_activity_window(annotations)
                annotations_rounded = round_up_to_full_seconds(annotations_active)
                annotations_5min = behavior_for_each_second(annotations_rounded)[["event ID"]]
                batch_name = path.split('\\')[-3]
                new_df_name = f"event ID_{mouse_name}"
                annotations_renamed = annotations_5min.rename(columns={'event ID': new_df_name})
                if concatenated_df is None:
                    concatenated_df = annotations_renamed
                else:
                    concatenated_df = pd.concat([concatenated_df, annotations_renamed], axis=1)  # Concatenate along columns
                print (f"{new_df_name} has {len(concatenated_df)} rows")
            concatenated_data[genotype] = concatenated_df.T.sort_index().T  # Store the accumulated data for this genotype
        all_batch[condition] = concatenated_data  
    return all_batch


def cumulative_timeline_events (all_results, behavior): 
    all_batch = {}
    median_dict ={}
    for condition, genotypes in all_results.items():
        concatenated_data = {} 
        concatenated_median ={}# Dictionary to store data for each genotype in the condition
        for genotype, mouse_pairs in genotypes.items():
            concatenated_df = None  # Initialize the DataFrame to accumulate data for this genotype
            median_list= []
            for mouse_name, path in mouse_pairs.items():
                annotations = open_behavior(path)
           
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "anogenital sniff", "sniff ano-gen")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "sit/still", "sit")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "rearing", "rear")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "self-groom", "groom")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "escape", "walk")
                annotations_active = cut_dataframe_to_activity_window(annotations)
                annotations_rounded = round_up_to_full_seconds(annotations_active)
                annotations_5min = behavior_for_each_second(annotations_rounded)[["event ID"]]
                new_annotations, median = count_behavioral_cum_events(annotations_5min, behavior)
                
                batch_name = path.split('\\')[-3]
                new_df_name = f"cumulative events_{batch_name}_{mouse_name}"
                annotations_renamed = new_annotations.rename(columns={'cumulative events': new_df_name})
                if concatenated_df is None:
                    concatenated_df = annotations_renamed
                else:
                    concatenated_df = pd.concat([concatenated_df, annotations_renamed], axis=1)  # Concatenate along columns
                median_list.append(median)
            concatenated_data[genotype] = concatenated_df.T.sort_index().T  # Store the accumulated data for this genotype
            concatenated_median [genotype] =pd.DataFrame(median_list, columns = ["median time"])
        all_batch[condition] = concatenated_data  
        median_dict[condition] = concatenated_median
    return all_batch, median_dict





def concat_all_processed_annotations_time (all_results):
    all_batch = {}
    for condition, genotypes in all_results.items():
        concatenated_data = {}  # Dictionary to store data for each genotype in the condition
        for genotype, mouse_pairs in genotypes.items():
            concatenated_df = None  # Initialize the DataFrame to accumulate data for this genotype
            annotations_active_dict = {}  # Dictionary to store annotations_active for each mouse_name

            for mouse_name, path in mouse_pairs.items():
                annotations = open_behavior(path)
                print(path)
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "anogenital sniff", "sniff ano-gen")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "sit/still", "sit")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "rearing", "rear")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "self-groom", "groom")
                annotations = replace_value_in_dataframe(annotations, "Behaviour", "escape", "walk")
                annotations_active = cut_dataframe_to_activity_window(annotations)

                annotations_active_dict[mouse_name] = annotations_active

            concatenated_data[genotype] = annotations_active_dict
            
        all_batch[condition]=concatenated_data
    return all_batch



def calculate_all_parameters (dictionary_time, behavior):
    latency_to_bouts_dict = {}
    duration_bouts_dict = {}
    interval_bouts_dict = {}
    half_time_dict ={}
    number_events_dict = {}
    total_time_dict = {}   
    for condition, genotypes in dictionary_time.items():
        latency_to_bouts_dict_batch = {}
        duration_bouts_dict_batch = {}
        interval_bouts_dict_batch = {}
        half_time_dict_batch ={}
        number_events_dict_batch = {}
        total_time_dict_batch = {}   
        for genotype, dict_animals in genotypes.items():
            latency_to_bouts_df = pd.DataFrame()
            duration_bouts_df = pd.DataFrame()
            interval_bouts_df = pd.DataFrame()
            half_time_list =[]
            number_events_list = []
            total_time_list = []  
            for all_animals, annotations in dict_animals.items():
                if behavior in annotations['Behaviour'].values:
                    latency_to_bouts = calculate_latency_to_bouts (annotations, behavior)['latency to start'].rename(f"latency to start_{all_animals}")
                    duration_bouts= calculate_bout_duration_behavior(annotations, behavior)['duration'].rename(f"duration_{all_animals}")
                    interval = calculate_latency_between_bouts(annotations, behavior)
                    half_time = calculate_half_time(annotations, behavior)
                    number_events = calculate_number_duration_totalevents (annotations, behavior)[0]
                    total_time = calculate_number_duration_totalevents (annotations, behavior)[1]
                    print ("behavior exists")
                else:
                    latency_to_bouts =  pd.DataFrame({f"latency to start_{all_animals}": ['no_behavior']})
                    duration_bouts= pd.DataFrame({f"duration_{all_animals}": ['no_behavior']})
                    interval = pd.DataFrame({f"interval_{all_animals}": ['no_behavior']})
                    half_time = 0
                    number_events = 0
                    total_time = 0
                    print ("no behavior")
                    # Add the values to the lists even when behavior is not found


                if latency_to_bouts_df is None:
                    latency_to_bouts_df = latency_to_bouts
                else:
                    latency_to_bouts_df = pd.concat([latency_to_bouts_df, latency_to_bouts], axis=1)
                    
                if duration_bouts_df is None:
                    duration_bouts_df = duration_bouts
                else:
                    duration_bouts_df = pd.concat([duration_bouts_df, duration_bouts], axis=1)
                    
                if interval_bouts_df is None:
                    interval_bouts_df = interval.rename(columns={'interval between bouts': f"interval_{all_animals}"})
                else:
                    interval_bouts_df = pd.concat([interval_bouts_df, interval.rename(columns={'interval between bouts': f"interval_{all_animals}"})], axis=1)
                
                half_time_list.append (half_time)
                half_time_df = pd.DataFrame({'half_time': half_time_list})
                
                number_events_list.append (number_events)
                number_events_df = pd.DataFrame({'number_events': number_events_list})
                
                total_time_list.append (total_time)
                total_time_df = pd.DataFrame({'total_time': total_time_list})
            
            latency_to_bouts_dict_batch[genotype] = latency_to_bouts_df.T.sort_index().T 
            duration_bouts_dict_batch[genotype] = duration_bouts_df.T.sort_index().T
            interval_bouts_dict_batch[genotype] = interval_bouts_df.T.sort_index().T
            half_time_dict_batch [genotype]=half_time_df
            number_events_dict_batch [genotype]=number_events_df
            total_time_dict_batch [genotype]=total_time_df
        
        latency_to_bouts_dict[condition]=latency_to_bouts_dict_batch
        duration_bouts_dict[condition] = duration_bouts_dict_batch
        interval_bouts_dict[condition] = interval_bouts_dict_batch
        half_time_dict[condition] =half_time_dict_batch
        number_events_dict[condition] =number_events_dict_batch
        total_time_dict[condition] =total_time_dict_batch
    
    return latency_to_bouts_dict, duration_bouts_dict, interval_bouts_dict, half_time_dict, number_events_dict, total_time_dict

    
behavior_colors = {
    "sniff ano-gen": "cyan",
    "sniff head-torso": "limegreen",
    "allogroom": "navy",
    "groom": "orange",
    "rear": "thistle",
    "dig": "saddlebrown",
    "walk": "darkgreen",
    "sit": "khaki",
    "fight": "lightcoral"
}

anogen_colors = {
    "sniff ano-gen": "cyan",
    "sniff head-torso": "white",
    "allogroom": "white",
    "groom": "white",
    "rear": "white",
    "dig": "white",
    "walk": "white",
    "sit": "white",
    "fight": "white"
}

headtorso_colors = {
    "sniff ano-gen": "white",
    "sniff head-torso": "limegreen",
    "allogroom": "white",
    "groom": "white",
    "rear": "white",
    "dig": "white",
    "walk": "white",
    "sit": "white",
    "fight": "white"
}


def graph_plotting (dictionary):
    behaviors = ["sniff ano-gen", "sniff head-torso", "allogroom", "groom", "rear", "dig", "walk", "sit", "fight"]
    for condition, genotypes in dictionary.items():
        for genotype, dataframe in genotypes.items():
            plot_frequency_graph(dataframe, behaviors, f"freq_{condition}_{genotype}")
            plot_behaviors(dataframe, f"rastaplot ALL_{condition}_{genotype}", behavior_colors)
            plot_behaviors(dataframe, f"rastaplot ANOGENITAL SNIFFING_{condition}_{genotype}", anogen_colors)
            plot_behaviors(dataframe, f"rastaplot BODY EXPLORATION_{condition}_{genotype}", headtorso_colors)
            plot_frequency(dataframe, f"freq heatmap_{condition}_{genotype}")
            
def generate_excels_by_time (dictionary, root_directory):
    behaviors = ["sniff ano-gen", "sniff head-torso", "allogroom", "groom", "rear", "dig", "walk", "sit", "fight"]
    for condition, genotypes in dictionary.items():
        for genotype, dataframe in genotypes.items():
            detect_behaviors(dataframe, behaviors)
            events, freq, cum_events = calculate_freq_and_events(dataframe)
            
            events_excel = pd.ExcelWriter(root_directory.joinpath(f'events ALL_{condition}_{genotype}.xlsx'))
            for key_events, df_events in events.items():
                df_events.to_excel(events_excel, sheet_name=key_events, index=False)
            events_excel.save()
            events_excel.close()
            
            freq_excel = pd.ExcelWriter(root_directory.joinpath(f'freq ALL_{condition}_{genotype}.xlsx'))
            for key_freq, df_freq in freq.items():
                 df_freq.to_excel(freq_excel, sheet_name=key_freq, index=False)
            freq_excel.save()
            freq_excel.close()
            
            norm_cum_freq_excel = pd.ExcelWriter(root_directory.joinpath(f'normalized_cum_EVENT ALL_{condition}_{genotype}.xlsx'))
            for key_cum_freq, df_cum_freq in cum_events.items():
                 df_cum_freq.to_excel(norm_cum_freq_excel, sheet_name=key_cum_freq, index=False)
            norm_cum_freq_excel.save()
            norm_cum_freq_excel.close()
            
def generate_excel_cum_event (dictionary, root_directory, behavior):
    dict_events, dict_median = cumulative_timeline_events (dictionary, behavior)

    for condition, genotypes in dict_events.items():
        for genotype, dataframe in genotypes.items():
            events_excel = pd.ExcelWriter(root_directory.joinpath(f'CUM events ALL__{behavior}_{condition}_{genotype}.xlsx'))
            dataframe.to_excel(events_excel, sheet_name=genotype, index=False)
            events_excel.save()
            events_excel.close()
    
    excel_median = pd.ExcelWriter(root_directory.joinpath(f'median time_{behavior}.xlsx'))        
    for condition1, genotypes1 in dict_median.items():
        for genotype1, dataframe1 in genotypes1.items():
            sheet_name1 = f"{condition1}_{genotype1}"
            dataframe1.to_excel(excel_median, sheet_name=sheet_name1, index=False)
    excel_median.save()
    
    
            
def generate_excels_parameters (dictionary, root_directory, behavior):
    dict_all = concat_all_processed_annotations_time(dictionary)
    dict_anogen = calculate_all_parameters (dict_all, behavior)
    latency_to_bouts_anogen=dict_anogen[0]
    duration_bouts_anogen = dict_anogen[1]
    interval_bouts_anogen = dict_anogen[2]
    half_time_anogen =dict_anogen[3]
    number_events_anogen =dict_anogen[4]
    total_time_anogen =dict_anogen[5]
    
    excel_latency = pd.ExcelWriter(root_directory.joinpath(f'latency_bouts_{behavior}.xlsx'))
    for condition_genotype, data in latency_to_bouts_anogen.items():
        for genotype, dataframe in data.items():
            sheet_name = f"{condition_genotype}_{genotype}"
            dataframe.to_excel(excel_latency, sheet_name=sheet_name, index=False)
    excel_latency.save()
    
    excel_duration = pd.ExcelWriter(root_directory.joinpath(f'duration_bouts_{behavior}.xlsx'))
    for condition_genotype_1, data1 in duration_bouts_anogen.items():
        for genotype1, dataframe1 in data1.items():
            sheet_name1 = f"{condition_genotype_1}_{genotype1}"
            dataframe1.to_excel(excel_duration, sheet_name=sheet_name1, index=False)
    excel_duration.save()
    
    excel_interval = pd.ExcelWriter(root_directory.joinpath(f'interval_{behavior}.xlsx'))
    for condition_genotype_2, data2 in interval_bouts_anogen.items():
        for genotype2, dataframe2 in data2.items():
            sheet_name2 = f"{condition_genotype_2}_{genotype2}"
            dataframe2.to_excel(excel_interval, sheet_name=sheet_name2, index=False)
    excel_interval.save()
    
    excel_half_time = pd.ExcelWriter(root_directory.joinpath(f'half_time_{behavior}.xlsx'))
    for condition_genotype_4, data4 in half_time_anogen.items():
        new_df4 = None
        for genotype4, df in data4.items():
            df.rename(columns={'half_time': genotype4}, inplace=True)
            if new_df4 is None:
                new_df4 = df
            else:
                new_df4 = pd.concat ([new_df4, df], axis=1)
        sheet_name4 = f"{condition_genotype_4}"
        new_df4.to_excel(excel_half_time, sheet_name=sheet_name4, index=False)
    excel_half_time.save()
    
    excel_number_events = pd.ExcelWriter(root_directory.joinpath(f'number_events_{behavior}.xlsx'))
    for condition_genotype_5, data5 in number_events_anogen.items():
        new_df5 = None
        for genotype5, df5 in data5.items():
            df5.rename(columns={'number_events': genotype5}, inplace=True)
            if new_df5 is None:
                new_df5 = df5
            else:
                new_df5 = pd.concat ([new_df5, df5], axis=1)
        sheet_name5 = f"{condition_genotype_5}"
        new_df5.to_excel(excel_number_events, sheet_name=sheet_name5, index=False)
    excel_number_events.save()
    
    excel_total_time = pd.ExcelWriter(root_directory.joinpath(f'total_time_{behavior}.xlsx'))
    for condition_genotype_6, data6 in total_time_anogen.items():
        new_df6 = None
        for genotype6, df6 in data6.items():
            df6.rename(columns={'total_time': genotype6}, inplace=True)
            if new_df6 is None:
                new_df6 = df6
            else:
                new_df6 = pd.concat ([new_df6, df6], axis=1)
        sheet_name6 = f"{condition_genotype_6}"
        new_df6.to_excel(excel_total_time, sheet_name=sheet_name6, index=False)
    excel_total_time.save()
    



#RUN CODE
dictionary_all = process_data(root_directory, summary_all)
dict_annotations = find_and_store_file(dictionary_all)

dictionary_processed_events = concat_all_processed_annotations_events (dict_annotations)
generate_excels_by_time (dictionary_processed_events,root_directory)

generate_excels_parameters (dict_annotations, root_directory, "sniff head-torso")
generate_excels_parameters (dict_annotations, root_directory, "sniff ano-gen")
generate_excels_parameters (dict_annotations, root_directory, "groom")
generate_excels_parameters (dict_annotations, root_directory, "dig")
generate_excels_parameters (dict_annotations, root_directory, "rear")
generate_excels_parameters (dict_annotations, root_directory, "sit")
generate_excels_parameters (dict_annotations, root_directory, "allogroom")
generate_excels_parameters (dict_annotations, root_directory, "walk")
generate_excel_cum_event (dict_annotations, root_directory, "sniff head-torso")
generate_excel_cum_event (dict_annotations, root_directory, "sniff ano-gen")
generate_excel_cum_event (dict_annotations, root_directory, "sniff head-torso")
generate_excel_cum_event (dict_annotations, root_directory, "groom")
generate_excel_cum_event (dict_annotations, root_directory, "dig")
generate_excel_cum_event (dict_annotations, root_directory, "rear")


graph_plotting (dictionary_processed_events)


