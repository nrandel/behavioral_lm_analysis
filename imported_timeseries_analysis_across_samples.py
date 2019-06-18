#%% [markdown]
# IMPORTANT: behavior annotation with Chris's fran code are zero-based!!!
# 
# addition to the code: 
# 
# 1) read in a single/ few data sets and read all datasets (for cellConfig and behavior_transition)
# 
# 2) check if I really have solved the timestamp issue
# 
# 3) check if all csv files exists
# 
# 4)!!! ratio to entcounter for less activity in the end
# 
# ------------------------------------------------------------------------------
#%% [markdown]
# This script was build for analysis of neuronal avtivity data in combination with behavioral annotations. 
# 
# Multiple csv-files (behavior, dff, time_stamp, neuronal activity (lm-data)) are combined in a single dataframe per experiment (sample_df). The identity of the sample is kept by introducing a sample-ID (date, number of sample) and an experiment-ID (imaging acquisition type <close and open loop>). Kind of behavior is predefined and extended with 'quiet'. Data of behavior, as well time-data are combined with all lm data. To define the start end end, as well accoiunt for overlapping events of the same behavior, a start, end, and overlap column is added for each behavior.

#%%
# Import libraries and implement error messages (Default should be True)
# To make the plot in the notebook and not in an extra window
get_ipython().run_line_magic('matplotlib', 'notebook')

import ast
import csv
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
import pandas as pd
import glob, os
import re

error_on_missing_timestamps = False
error_on_time_light_mismatch = False
error_on_time_behavior_mismatch = False
error_on_missing_behaviors = False
error_on_invalid_behavior_range = False


#%%
#Open multiple .csv from single directory. Define existing behaviors. Define sample_ID and experiment_ID.

#directory for behavior data
behavior_directory = r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/behavior_csv/' # directory for behavioral annotations
behavior_files = glob.glob(os.path.join(behavior_directory, "*.csv")) #join pathname with filename

# Behavior columns available in CSV files
available_behaviors = ('fw', 'bw', 'stim', 'hunch', 'turn', 'other', 'HP', 'left turn', 'right turn')

#Regular expression (define the expression filenames are searched for)
#'.' single character, matched everything, '*' 0>> occurences, '/' path delimiter, '\d' 0-9 digit,
#'+' 1>> occurences, 'L' here character from filename
#() outcome here: 2 groups, useful for extraction
#[] optional list, eg 1 --> 1
#? character non or once 
#Behavior reg-ex (regular expression)
behavior_sample_re = re.compile('.*/(\d\d-\d\d-\d\dL\d+(-\d+)?)-behavior-(.+).csv')

# Function: readall_behavior iterates through all csv (sorted) 
# and appends the files into the list (ls) and returns dictionary
def readall_behavior(all_files, printit=False):
    data = {}
    for filename in sorted(all_files):
        # Find sample ID, file name pattern: YY-MM-DDLXDETAIL.csv,
        # exp_id = DETAIL: several measurements of same sample 
        # (cl (closeloop, RGECO/ Chronos), ol (openloop, RGECO/ Chronos), 
        # blocks (Raghav: GCaMP/Chrimson))
        # Larva ID: YY-MM-DDLX
        # Look for filename_components, which are true for pattern
        match = behavior_sample_re.match(filename)
        if not match:
            raise ValueError('Unexpected filename format: {}'.format(filename))
        filename_components = match.groups()
        #define filename_components sample_id (first group), and exp_id (sec group)
        part_sample_id, _, exp_id = filename_components         
        sample_id = "{}-{}".format(part_sample_id, exp_id)
        
        df = pd.read_csv(filename, index_col=None, header=0, delimiter = ';')
        df.fillna(0, inplace=True) #replace NaN with zero
        df['sample_id'] = sample_id  #add sample_id column
        df['exp_id'] = exp_id #add exp_id column
        data[sample_id] = df
        #Count 'True' for each column ('behavior') in each single behavior.csv)
        #print(filename, df[df == 1].count()) 
        #print(df)
    return data

from functions import readall_behavior
behavior_data = readall_behavior(behavior_files)
#print(behavior_data['17-08-26L3-cl'])
print(list(behavior_data.keys()))


#%%
# Frequency of each behavior in all imported behavior.csv by using the returned 'ls' from 
# the function readAll: concatenate the 'behavior_files' (global variable). 'True' for each 
# column ('behavior_type') in the concatenated file (df_behavior).
# Sorting has to be = False (warning message without 'sort')
df_behavior = pd.concat(behavior_data.values(), axis = 0, ignore_index = True, sort = False) #add sorting
print(df_behavior[df_behavior == 1].count()) 


#%%
# Import and merge fluorescence data: Several LM files for the same sample_id exists, but differ in cell_id).
# List of LM data with two extra columns: sample_id and cell_id

    #Mapping of sampleID vs lists of LM dataframes
    #Convert list to a single dataframe
    #Map with df_behavior (later done)

# Open LM files from different directories
lightmicroscope_directories = [r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/Basin_traces/', 
                               r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/Handle-like_Traces',
                               r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/a00c_traces'
                              ] 

# Iterate through LM data and extend files in a list from within and between directory and 
# build a list of files from all directories
# (Note: append would 'extend' lists and not single files)
lightmicroscope_files = []
for d in lightmicroscope_directories:
    lightmicroscope_files.extend(
        glob.glob(os.path.join(d, "*.csv"))) #join pathname with filename, 
                                                                       
# Regular expression (define the expression filenames are searched for)
# '.' single character, matched everything, '*' 0>> occurences, '/' path delimiter, '\d' 0-9 digit,
# '+' 1>> occurences, 'L' here character from filename
# () outcome here: 2 groups, useful for extraction

# Lightmicroscopic data reg-ex (regular expression)
lightmicroscope_sample_re = re.compile('.*/(\d\d-\d\d-\d\dL\d+(-\d+)?)-(.*)-(.*).csv')

# Function: readall_lm iterates through all LM_csv (sorted) 
# and returns a dictionary{key:value} 
# samples = {sample_id:cell-id}
def readall_lm(all_files):
    samples = {}
    for filename in sorted(all_files):
        # Find sample ID, file name pattern: YY-MM-DDLXDETAIL.csv,
        # Larva ID: YY-MM-DDLX, DETAIL = cell_id
        #look for filename_components, which are true for pattern
        match = lightmicroscope_sample_re.match(filename)
        if not match:
            raise ValueError('Unexpected filename format: {}'.format(filename))
        filename_components = match.groups()
        part_sample_id, _, cell_id, exp_id = filename_components
        
        sample_id = "{}-{}".format(part_sample_id, exp_id)
        
        # Read LM.files 
        df = pd.read_csv(filename, index_col=None, header=0, delimiter = ',')
        # Replace NaN with zero
        df.fillna(0, inplace=True)
        
        # Add cellname to each column as prefix
        # lambda is a non defined function (longer version: def lambda(x):)
        # Rename of columns after the format cell_id, name) eg: Basin A9
        # inplace = True: column names are overwritten (if False: new dataframe)
        df.rename(lambda x: '{}_{}'.format(cell_id, x), axis = 'columns', inplace = True)
        
        # Get the sample_id (key) from the dictionary? to make a list [sample_cells] and 
        # if sample_id exists, append the list
        # if sample_id does not exists, start a new list
        # reminder: there can be several cell_id per sample_id
        sample_cells = samples.get(sample_id)
        if not sample_cells:
            samples[sample_id] = sample_cells = {
                'data': [],
                'exp_id': exp_id,
            }
        sample_cells['data'].append(df)
        
    return samples

lm_samples = readall_lm(lightmicroscope_files)

# New dictionary: lm_data{} to build a single dataframe with all cell_ids combined 
# for a single sample. Iterate over dict (samples)? and data from same sample in 
# one dataframe. 
# df.items iterate over pairs and build a list

lm_data = {}

# Iterate over all light samples and merge all found files
# for each sample into a single data frame (per sample)
for sample_id, sample_info in lm_samples.items():
    cells_dataframes = sample_info['data']
    #check if number of cells >= 1
    if not cells_dataframes:
        raise ValueError('No cells found for sample {}'.format(sample_id))
    #first element in the list
    lm_df = None

    #iteration through other df
    for cdf in cells_dataframes:
        if lm_df is None:
            lm_df = cdf
        else:
            if len(lm_df.index) != len(cdf.index):
                raise ValueError('Data frame frame to merge has not same row count as target', sample_id)
            lm_df = pd.merge(lm_df, cdf, left_index = True, right_index = True)
            
    lm_df['sample_id'] = sample_id  #add sample_id column
    lm_df['exp_id'] = sample_info['exp_id']
    lm_data[sample_id] = lm_df
#print(list(lm_data.keys()))


#%%
# Import txt-files from of the absolute time/frame from the Ca-imaging (lm-data). 
# All txt-files have to be transposed, which is a memory intensive step. After the 
# data are complete, the transposed files should be exported (ToDo). Time-data are 
# combined with sample-ID and experiment-ID.

timelapse_directory =(r'/Users/randeln/Documents/Zlatic_lab/close-loop/Notes/timelapse/') 
timelapse_files = glob.glob(os.path.join(timelapse_directory, "*.txt")) #join pathname with filename


# Regular expression (define the expression filenames are searched for)
# '.' single character, matched everything, '*' 0>> occurences, '/' path delimiter, '\d' 0-9 digit,
# '+' 1>> occurences, 'L' here character from filename
# () outcome here: 2 groups, useful for extraction
# [] optional list, eg 1 --> 1
# ? character non or once 

# Behavior reg-ex (regular expression)
time_sample_re = re.compile('.*/(\d\d-\d\d-\d\dL\d+(-\d+)?)-time-(.+).txt')

# Function: readall_timelapse iterates through all txt (sorted) and appends the 
# files into the dict (data) and returns ls
def readall_time(all_files, printit=False):
    data = {}
    for filename in sorted(all_files):
        # Find sample ID, file name pattern: YY-MM-DDLXDETAIL.csv,
        # exp_id = DETAIL: several measurements of same sample (cl (closeloop), ol (openloop), blocks (Raghav))
        # Larva ID: YY-MM-DDLX
        #look for filename_components, which are true for pattern
        match = time_sample_re.match(filename)
        if not match:
            raise ValueError('Unexpected filename format: {}'.format(filename))
        filename_components = match.groups()
        part_sample_id, _, exp_id = filename_components #define filename_components sample_id (first group), and exp_id (sec group)  
        sample_id = "{}-{}".format(part_sample_id, exp_id)
        
        df = pd.read_csv(filename, header=1, index_col=None, delim_whitespace = True)
        df = df.T #transposing because read_csv imports as row
        df = df.reset_index() #transpose function sets data as index
        df.rename(columns={'index':'time'}, inplace=True) #rename reset index column to time
        df['time'] = df.time.astype(float)
        data[sample_id] = df
        
    return data


#%%
timelapse_cache = 'timelapse.cache'

try:
    with open(timelapse_cache, 'r') as timelapse_cache_file:
        # TODO
        cache_data = timelapse_cache_file.read()
        time_data = ast.literal_eval(cache_data)
except FileNotFoundError as e:
    print('No cache file found, recomputing')
    # No cache file found, recompute
    time_data = readall_time(timelapse_files)
    # Write cache
    


#%%
sample_data = {}

# Time data are merged into light data and checked if number length of lm = timestamp.  
# Due to technical conditions, some time.txt-file have too many or not enough time data compared
# to the corresponding LM data. The discrepancy is fixed by either dropping the extra timepoints or 
# by taking the average of the difference between each timepoint and extend the dataframe. 
# The first 10 timepoints are not included to account for instability of the microscope in 
# the beginning due to the moving parts. 
# Maximal difference between timepoints fyi.

for sample_id, sample_df in lm_data.items():
    # Add time stamps to data frame of current sample by merging
    # The time data frame for the current sample, which is expected
    # to match the light data (based on index).
    timestamp_df = time_data.get(sample_id)
    if timestamp_df is None:
        msg = '{}: could not find timestamp data for sample'.format(sample_id)
        if error_on_missing_timestamps:
            raise ValueError(msg)
        # Ignore, if missing data shouldn't cancel the whole process.
        print(msg)
        continue
        
    n_timestamps = len(timestamp_df)
    n_lightdata = len(sample_df)
    
    # The timestamp and light recordings are done by different systems.
    # This can cause the existence of additional time points/ or missing time points in a
    # dataset, which will be filtered out in the merge operation below.
    if n_lightdata != n_timestamps:
        msg = '{}: time data ({} entries) doesn\'t match light data ({} entries)'.format(
                sample_id, n_timestamps, n_lightdata)
        if error_on_time_light_mismatch:
            raise ValueError(msg)
        print(msg)
        diffs = np.diff(timestamp_df['time'])[10:] #from 10th row onwards
        diffs_avg = diffs.mean(axis=0)
        #diff between timedata and lightdata
        missing_data = len(sample_df) - len(timestamp_df)
        
        #add 'diffs_avg' to fill in missing_timedata
        if missing_data > 0:
            last_valid_index = len(timestamp_df) - 1
            last_timestamp = timestamp_df.iloc[last_valid_index]['time']
            if pd.isna(last_timestamp):
                raise ValueError('Unexpected last valid timestamp for sample {} at index {}'.format(
                        sample_id, last_valid_index))
            for i in range(0, missing_data):
                last_valid_index += 1
                timestamp_df.loc[last_valid_index] = timestamp_df.iloc[last_valid_index - 1]['time'] + diffs_avg
        elif missing_data < 0:
            drop_start = len(timestamp_df) + missing_data
            drop_end = len(timestamp_df)
            timestamp_df.drop(list(range(drop_start, drop_end)))

    # Merge timedata into light data
    # Use an 'inner' join/merge to exclude time points that don't have matching light data.
    new_sample_df = pd.merge(sample_df, timestamp_df, left_index = True, right_index = True, how='inner')
    
    # Store newly created data frame for sample (dictionary)
    sample_data[sample_id] = new_sample_df
    
print('Matched {} light data sets with their respective time points'.format(len(sample_data)))

# Max.diffs for timestamps
diffs = np.diff(timestamp_df['time'])[10:] #from 10th row onwards
mx = diffs.max()
#print(mx)


#%%
# Combine behavior data with light data into a single data frame
# per sample ID. To do so, add behavior data to light data frames,
# because the light data is already organizes by frame. To accomodate
# frame ranges without an behavior data, a column named "quiet" is
# added which is True in these cases and False otherwise. Additionally,
# for each behavior column, a behavior start and end column as well as
# an overlap column is added so that parallel and successive behaviors
# of the same type can be differentiated.

for sample_id, sample_df in sample_data.items():
    sample_behavior = behavior_data.get(sample_id)
    if sample_behavior is None:
        msg = 'Could not find behavior data for sample "{}"'.format(sample_id)
        if error_on_missing_behaviors:
            raise ValueError(msg)
        print(msg)
        continue

    # Add extra columns for behavior
    for behavior in available_behaviors:
        sample_df[behavior] = False
        sample_df['{}_start'.format(behavior)] = False
        sample_df['{}_end'.format(behavior)] = False
        sample_df['{}_overlap'.format(behavior)] = False
    
    # Add 'quiet' column. Set it initially to True and mark frames
    # with actual behavior as quiet = False.
    sample_df['quiet'] = True
    
    n_light_entries = len(sample_df)

    # Iterate over behavior data and add data to target data frame
    for i, row in sample_behavior.iterrows():
        # Start ane end are 1-based, make them 0-based
        start = int(row['START'])
        end = int(row['END'])
        
        if type(row['START']) == str:
            print(sample_id)
            print(start, end)
        
        if start >= end:
            msg = "{}: start ({}) needs to be strictly smaller than end ({})".format(sample_id, start, end)
            if error_on_invalid_behavior_range:
                raise ValueError(msg)
            print(msg)
            continue
        
        # Make sure we capture start/end times that are a fractional number.
        if row['START'] - start > 0 or row['END'] - end > 0:
            raise ValueError('{}: start and end frame number can\'t contain fractions'.format(sample_id))
            
        # Ignore behavior entries with an end frame higher than available light data.
        # The behavior data is one-based, which is why a strict larger than test should
        # be correct.
        if end > n_light_entries:
            msg = 'Sample: {} - Behavior row with range {}-{} exceeds light time points ({}): {}'.format(
                sample_id, start, end, n_light_entries, row)
            if error_on_time_behavior_mismatch:
                raise ValueError(msg)
            print(msg)
            continue
            
        # Find behavior observed in row
        observed_behaviors = []
        for behavior in available_behaviors:
            if row[behavior]:
                observed_behaviors.append(behavior)
        
        # We assume that not more than two behaviors are observed at the same time
        if len(observed_behaviors) > 2:
            raise ValueError('Found multiple behaviors in row {} of sample {}'.format(i, sample_id))
        
        # Add observed behavior information to target data frames in all
        # rows in behavior range.
        for b in observed_behaviors:
            # Iterate over frames valid for current behavior. Every valid
            # frame is mapped into the canonical (light/cell) data frame,
            # which is 0-indexed.
            for j in range(start, end + 1):
                # Behavior ranges are 1-indexed
                current_frame = j - 1
                # If the current behavior has already been observed at this frame,
                # set overlap to True, because we are about to mark this behavior
                # again as observed for this frame.
                if sample_df.at[current_frame, b]:
                    sample_df.at[current_frame, '{}_overlap'.format(b)] = True
                else:
                    sample_df.at[current_frame, b] = True
                
                # Mark this row as not quiet, because we observed
                # a behavior in the current frame.
                sample_df.at[current_frame, 'quiet'] = False

            sample_df.at[start - 1, '{}_start'.format(b)] = True
            sample_df.at[end - 1, '{}_end'.format(b)] = True
            
    # Mark quiet ranges with _start, _end and _overlap. By definion,
    # quiet_overlap is always False.
    sample_df['quiet_start'] = False
    sample_df['quiet_end'] = False
    sample_df['quiet_overlap'] = False
    last_sample_idx = n_light_entries - 1
    for i, row in sample_df.iterrows():
        sample_df.at[i, 'quiet_start'] = row['quiet'] and (i == 0 or not sample_df.at[i - 1, 'quiet'])
        sample_df.at[i, 'quiet_end'] = row['quiet'] and (i == last_sample_idx or not sample_df.at[i + 1, 'quiet'])
        print(sample_df)

#%% [markdown]
# The above generated dataframe per sample (sample_df) including lm data, behavior, time and sample_id/exp_id, and the extended cell-names and extra behavioral columns can be analysed in the following section.
# 
# Class will be defined, including sample_id, cell_type, event (type), and filter pattern. This 
# allows to extract information throughout all samples about activity pattern (lm-data) 
# of a specific celltype (including sub-type) and event_start (=behavior). 
# 
# For single sample: cell_subset_df 
#         allows visualisation of whole experimantal time, and extract/ visualized information around a specific
#         event, which are aligned and normalized (event_start = zero). An adjustable time-window around the 
#         event_start of interested is included. 
#     - avg, max, min, stdev and sem can be direct calculated from the data (for whole timeseries)
#     - But if events are aligned, same problem as with multiple samples
#         
# For multiple combined samples (no additional processing): all_events
#         extract/ visualized information around a specific
#         event, which are aligned and normalized (event_start = zero). An adjustable time-window around the
#         event_start of interested is included. 
# 
#     The samples are imaged (Ca-imaging) with different imaging speed, and therefore a direct comparison between the 
#     samples is not possible. Following things have to be considered:
#     
#                     ToDo!!!
#     
#         - If the data are analysed as raw data (all_events) there are many NaN in the data set, and also the same 
#         timestamp is duplicated for each individual sample
#         - avg, max, min, stdev and sem can not be direct calculated from the data, because(?) vergessen:(
#         
#         - For some analysis (), the samples should be aligned: For now we use interpolation after index. 
#         THAT HAS TO BE CHECKED!!!  
#         
#         - alternatives: underlying fitting curve

#%%
#Define a class with sample_id, cell_type, event_name and filter_pattern
class CellTraceConfig:
    
    def __init__(self, sample_id, cell_type, event_name, filter_pattern=None):
        self.sample_id = sample_id
        self.cell_type = cell_type
        self.event_name = event_name
        self.filter_pattern = filter_pattern
        
    def get_filter_regex(self):
        filter_regex = '^{}_'.format(self.cell_type)
        if self.filter_pattern:
            filter_regex += '.*{}.*'.format(self.filter_pattern)
        return filter_regex
    
    def get_event_start_col(self):
        return '{}_start'.format(self.event_name)

    def add_event_time_points_to_plot(self, source_df, plot):
        for idx, row in source_df.iterrows():
            plot.annotate(self.event_name, xy=(row['time'], 1))
            plt.axvline(row['time'], color='k', linestyle='-')  
            
#Define a class with sample_id, cell_type, event_time and filter_pattern
class CellTransConfig:
    
    def __init__(self, sample_id, cell_type, event_time, filter_pattern=None):
        self.sample_id = sample_id
        self.cell_type = cell_type
        self.event_time = event_time
        self.filter_pattern = filter_pattern
        
    def get_filter_regex(self):
        filter_regex = '^{}_'.format(self.cell_type)
        if self.filter_pattern:
            filter_regex += '.*{}.*'.format(self.filter_pattern)
        return filter_regex


#%%
# Allows to load specific samples (single samples) with specific filter pattern

cell_trace_configs = [
    #CellTraceConfig('17-08-26L6-cl', 'basin', 'stim'),
    CellTraceConfig('17-08-26L5-cl', 'A00c', 'stim', 'mid'),
    CellTraceConfig('17-08-26L2-cl', 'A00c', 'stim', 'mid'),
    #CellTraceConfig('17-08-23L2-cl', 'A00c', 'stim', 'mid'),
    CellTraceConfig('17-08-26L6-cl', 'A00c', 'stim', 'mid'),
    #CellTraceConfig('17-08-24L2-1-cl', 'A00c', 'stim'),
    #CellTraceConfig('17-08-24L2-2-cl', 'A00c', 'stim', 'mid'),
    #CellTraceConfig('17-08-24L5-cl', 'A00c', 'stim', 'mid')
]
'''
# Allows to load all samples with specific filter pattern
cell_trace_configs = [
    CellTraceConfig(name,'A00c', 'fw') for name in lm_data]
'''
#put '' [empty string] if you dont want any cell type

# Allow to load all samples with specific filter pattern
# ToDo

all_events = [] #List of events, with raw dff data (no interpolation or other 
                #processing done at this point). Sample_id is added to the cell name. 

for ctc in cell_trace_configs:
    sample_df = sample_data.get(ctc.sample_id)
    if sample_df is None:
        raise ValueError('{}: could not find sample data'.format(ctc.sample_id))
        continue    
    # Extract columns matching our cell type and the optional filter pattern.
    # Pandas' filter() operations works on columns for DataFrames by default.
    cell_subset_df = sample_df.filter(regex=ctc.get_filter_regex()) #Get subset of cells 
    cell_subset_df.set_index(sample_df.time, inplace=True) #Set time to index (essential for min/max...)
    cell_subset_df.reset_index(inplace = True) # Add index and time = column
    #print(ctc.sample_id, cell_subset_df)   
    # Get rows where current event starts.
    event_df = sample_df[sample_df.loc[:,ctc.get_event_start_col()]]
    # Gives the timestamp for the event_df (start)
    for idx, row in event_df.iterrows():
        print('TP of {} ='.format(ctc.event_name), row['time'])
        
    # Extract for specific time window and align several events. 
    # Define timepoints pre and post an event (event_df). 
    # This works for single sample or multiple samples aligned 
    # Note: In cell_subset_df, time was set to index, because for the previous avg calculation 
    # Add index and time = column

    # Set the window range left and right from the event
    left_half_window_size = 10.0 #in seconds
    right_half_window_size = 50.0

    # Event_df defined in pargraph before 
    windows = []
    n_behavior = 0
    for i,row in event_df.iterrows():
        n_behavior += 1
        window_start = row['time'] - left_half_window_size
        window_end = row['time'] + right_half_window_size
        
        # Get subset of rows between window_start and window_end       
        event = cell_subset_df[(cell_subset_df.time >= window_start) & (cell_subset_df.time <= window_end)]
        # Normalizing the data to align on beginning of selected
        # behavior (event_df = Zero) by substracting events in window
        # around start of event of interest from start of event interest.
        # Note: using ":" in event.loc[] will select "all rows" in our window.
        event.loc[:, 'time'] = event['time'] - row['time']

        # Add sample_id to each column as prefix and n_behavior as suffix to distinguish events within a sample
        event.rename(lambda x: '{}_{}_{}'.format(ctc.sample_id, x, n_behavior), 
                     axis = 'columns', inplace = True) 

        # Rename time collum to time
        event.rename(columns={ event.columns[0]: 'time' }, inplace = True) 
        all_events.append(event) # Append a list with all event
        
        #Round (NR)
        #decimals = 1    
        #event['time'] = event['time'].apply(lambda x: round(x, decimals))
        
        
# Removes first event and takes it as left_window in pd.merge_ordered and iterates than through all_events
all_df = all_events.pop(0)
for right_df in all_events:
    all_df = pd.merge_ordered(all_df, right_df, on="time", how="outer")

# Resets the index as time and drops time column (sollte spaeter kommen)
all_df.index = all_df["time"]
del all_df["time"]        
#print(all_df)

# Index intepolation (linear interpolatione not on all_df, because index [=time] is not eaqually distributed)
int_all_df = all_df.interpolate(method='index', axis=0, limit=None, inplace=False, limit_direction='both')
#print(int_all_df)

#%% [markdown]
# single sample, over whole experimental time
# - does not have much meaningfulness
# 
# single sample, aligned for events
# ...

#%%



#%%
# Single sample -analysis
# For single sample over the whole experimental time
# Note: multiple sample-comparison need more pre-processing(see below)
# Calculate min, max, avg, stddev, sem from cell_subset_df (defined earlier)
cell_subset_df.set_index(sample_df.time, inplace=True) #Set time to index (essential for min/max...)
del cell_subset_df['time'] # delete time_column
cell_subset_df.index.name = None # delete index name
cell_avg_df = cell_subset_df.mean(axis=1)
cell_min_df = cell_subset_df.min(axis=1)
cell_max_df = cell_subset_df.max(axis=1)
# Standard deviation (distribution)
cell_std_df = cell_subset_df.std(axis = 1)
#standard error of mean
cell_sem_df = cell_subset_df.sem(axis = 1)
#print(ctc.sample_id, cell_avg_df) #OK


# For single or multiple sample, aligned for certain event
#Average is wrongly applied, because it avg all events and all cells pro tp
# Good! NaN are ignored and the correct avg is calculated
all_cell_avg_df = int_all_df.mean(axis=1) # Interpolated data used
all_cell_min_df = int_all_df.min(axis=1)
all_cell_max_df = int_all_df.max(axis=1)
# Standard deviation (distribution)
all_cell_std_df = int_all_df.std(axis = 1)
#standard error of mean
all_cell_sem_df = int_all_df.sem(axis = 1)
#print(all_cell_avg_df) #wrong zur haelfte: Want to have avg per celltyp over time point, 
                    #and not avg over all cells per timepoint


#%%
get_ipython().run_line_magic('matplotlib', 'notebook')
from matplotlib import pyplot as plt
plt.rcParams["axes.grid"] = False

# Plotting - single sample 
def layout_plot(plot, tick_spacing=100, fov=(0, 2500, 0, 1.2), legend=False): 
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False) 
    plot.axis(fov)
    plot.legend().set_visible(legend)

# Get rows where current event is active and draw a vertical 
# line to indicate the event in the plot
event_df = sample_df[sample_df.loc[:,ctc.get_event_start_col()] == 1]
fig = plt.figure()
fig.set_facecolor("white")

# Plot all cells from cell_subset_df over entire time (specified in Cell_Trace_Config).
sub1 = fig.add_subplot(211) #211
cell_subset_df.plot(ax=sub1)
ctc.add_event_time_points_to_plot(event_df, sub1)
layout_plot(sub1)

# Avg, min, max, std-dev for multiple cells in single sample over whole time
sub2 = fig.add_subplot(212)#212
ctc.add_event_time_points_to_plot(event_df, sub2)
cell_avg_df.plot(ax=sub2, color = 'g', label = ctc.cell_type, linewidth=1)
cell_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
cell_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
#cell_avg_df.plot.line(yerr=cell_std_df, ax=sub2, color = 'r', alpha = 0.1)
#cell_avg_df.plot.line(yerr=cell_sem_df, ax=sub2, color = 'c', alpha = 0.1)
layout_plot(sub2)


#%%
# Note: HERE FOR PLOTTING THE ALIGNED EVENT, INDEPENDENT OF PRO AND / OR POST_EVENT
# (should be after transition_event)

# Plotting for multi-events (all_df) (raw_dff_data)
# If a dataframe with NANs is plotted, use 
# marker = '+', or 'o', since the line in the lineplot only connects 
# consecutive data points
def aligned_layout_plot(plot, tick_spacing=0.5, fov=(-20, 50, -0.05, 1.9), legend=False): 
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False) 
    plot.axis(fov)
    plot.legend().set_visible(legend)

fig = plt.figure()

# Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
sub1 = fig.add_subplot(211)
all_df.plot(ax=sub1, marker = '*', label = ctc.cell_type)
aligned_layout_plot(sub1)

sub2 = fig.add_subplot(212)
all_cell_avg_df.plot(ax=sub2, color = 'k', label = ctc.cell_type) #use interpolated df to calculate average...
#all_cell_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
#all_cell_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
all_cell_avg_df.plot.line(yerr=all_cell_std_df, ax=sub2, color = 'r', alpha = 0.1)
#all_cell_avg_df.plot.line(yerr=all_cell_sem_df, ax=sub2, color = 'c', alpha = 0.1)
aligned_layout_plot(sub2)



#%%


#%% [markdown]
# The following part extract the information about behavior transition under certain limitation:
# 1) Find second_behavior, and extract information if the defined first_behavior happens within a max_delay.
# 2) Find first_behavior, and extract information if the defined second_behavior happens within a max_delay.
# 3) Find second_behavior, and extract information if the defined first_behavior and third_behavior 
#    happens within a max_delays.
# 4) Find first_behavior, and extract information if the same second_behavior happens within a max_delay. Note: So far no overlap cases detected. Code for overlap cases could not be verified.

#%%

class PreBehaviorTransition:
    
    def __init__(self, sample_id, pre_event, event, max_delay=0):
        self.sample_id = sample_id
        self.pre_event = pre_event
        self.event = event
        self.max_delay = max_delay

def find_behavior_after(sample_id, sample_df, first_event, second_event, max_delay=0):
    """For the data frame of a single sample <df>, find all behaviors
    of type <second_event> that follow the event <first_event>,
    separated by <max_delay> time. The start of <first_event> is expected
    to happen strictly before the start of <second_event>. The end time
    of <first_event> however can overlap with the start time of <second_event>.
    In this case, the time difference is negative, and still smaller than
    <max_delay>. The end time of <first_event> can be before, at or after the
    end of <second_event>.
    
    TODO: If <first_event> and <second_event> are the same type of behavior,
    overlaps have to be taken into account to match start and end times
    to the correct event.
    """
    results = []
    first_event_start_col = '{}_start'.format(first_event)
    first_event_end_col = '{}_end'.format(first_event)
    first_event_overlap_col = '{}_overlap'.format(first_event)
    second_event_start_col = '{}_start'.format(second_event)
    second_event_end_col = '{}_end'.format(second_event)
    
    first_event_start_time = None
    first_event_end_time = None
    second_event_start_time = None
    second_event_end_time = None
    
    for i, row in sample_df.iterrows():
        # Look for end of first behavior and remember its time.
        if row[first_event_start_col]:
            first_event_start_time = row['time']
        if row[first_event_end_col] and not row[first_event_overlap_col]:
            first_event_end_time = row['time']
        if row[second_event_start_col]:
            second_event_start_time = row['time']
        if row[second_event_end_col]:
            second_event_end_time = row['time']
        
        # As long as we haven't collected all needed time points,
        # keep on searching.
        if None in (first_event_start_time, first_event_end_time,
                    second_event_start_time, second_event_end_time):
            continue
            
        #NR
        # Define rules for event_start_time and event_end_time
        if first_event_start_time > second_event_start_time:
            continue
        if first_event_start_time > first_event_end_time:
            continue
        if second_event_start_time > second_event_end_time:
            continue
        
            
        if abs(first_event_start_time - second_event_start_time) < 0.00001:
            print('{}: start time (first) event {} and start time of (second) event {} are the same: {}'.format(
                sample_id, first_event, second_event, first_event_start_time))

        # Test time between first event end and second event start. If it
        # is smaller than <max_delay>, store start of second event as result.
        # The first event end time being greater than the second event start
        # time, is explicitly allowed.
        if (second_event_start_time - first_event_end_time) <= max_delay:
            results.append({
                'sample_id': sample_id,
                'first_event_start': first_event_start_time,
                'first_event_end': first_event_end_time,
                'second_event_start': second_event_start_time,
                'second_event_end': second_event_end_time
            })
        
        # Reset behavior tracking variables to find new pattern match.
        first_event_start_time = None
        first_event_end_time = None
        second_event_start_time = None
        second_event_end_time = None
            
    return results

behavior_transitions = [
    PreBehaviorTransition('17-08-26L1-cl', 'turn', 'bw', 11),
    #PreBehaviorTransition('17-08-26L6-cl', 'turn', 'bw', 3)
]

'''
# Open all samples (!see CellConfig!) >ToDo
behavior_transitions = [
    PreBehaviorTransition(name,'fw', 'bw', 10) for name in lm_data]

'''

found_transitions = []
for bt in behavior_transitions:
    sample_df = sample_data.get(bt.sample_id)
    if sample_df is None:
        raise ValueError('No data found for sample {}'.format(bt.sample_id))
    transitions = find_behavior_after(bt.sample_id, sample_df, bt.pre_event, bt.event, bt.max_delay)
    found_transitions.append(transitions)

print(len(found_transitions)) #number of data sets not the actual stim
print(len(transitions)) #not what I want ToDo
print(found_transitions)


#%%
class PostBehaviorTransition:
    
    def __init__(self, sample_id, event, post_event, max_delay=0):
        self.sample_id = sample_id
        self.post_event = post_event
        self.event = event
        self.max_delay = max_delay

def find_behavior_before(sample_id, sample_df, first_event, second_event, max_delay=0):
    """For the data frame of a single sample <df>, find all behaviors
    of type <first_event> that is followed by the event <second_event>,
    separated by <max_delay> time. The end of <second_event> is expected
    to happen strictly after the end of <first_event>. The start time
    of <second_event> however can overlap with the end time of <first_event>.
    In this case, the time difference is negative, and still smaller than
    <max_delay>. The start time of <second_event> can be before, at or after the
    end of <first_event>.
    
    TODO: If <first_event> and <second_event> are the same type of behavior,
    overlaps have to be taken into account to match start and end times
    to the correct event.
    """
    results = []
    first_event_start_col = '{}_start'.format(first_event)
    first_event_end_col = '{}_end'.format(first_event)
    second_event_start_col = '{}_start'.format(second_event)
    second_event_end_col = '{}_end'.format(second_event)
    second_event_overlap_col = '{}_overlap'.format(second_event)
    
    first_event_start_time = None
    first_event_end_time = None
    second_event_start_time = None
    second_event_end_time = None
    
    for i, row in sample_df.iterrows():
        # Look for start of second behavior and remember its time.
        if row[second_event_start_col] and not row[second_event_overlap_col]:
            #print("{} starts at {}".format(second_event, row["time"]))
            second_event_start_time = row['time']
        if row[first_event_end_col]:
            #print("{} ends at {}".format(first_event, row["time"]))
            first_event_end_time = row['time']
        if row[first_event_start_col]:
            #print("{} starts at {}".format(first_event, row["time"]))
            first_event_start_time = row['time']
        for column in sample_df.columns:
            if (first_event_start_time is not None and
                column.endswith("_start") and
                column != first_event_start_col and
                column != second_event_start_col and
                first_event not in column and
                second_event not in column):
                if row[column]:
                    #print("{} ended at {}, but then found {} at {}".format(first_event, first_event_end_time, column, row["time"]))
                    first_event_start_time = None
                    first_event_end_time = None
                    second_event_start_time = None
                    second_event_end_time = None
                
        
        # As long as we haven't collected all needed time points,
        # keep on searching.
        if None in (first_event_start_time, first_event_end_time,
                    second_event_start_time):
            continue
        
        #NR
        # Define rules for event_start_time and event_end_time
        if first_event_start_time > second_event_start_time:
            continue
        if first_event_start_time > first_event_end_time:
            continue
        
        # Test if first_event_start_time = second_event_start_time
        if abs(first_event_start_time - second_event_start_time) < 0.00001:
            print('{}: start time (first) event {} and start time of (second) event {} are the same: {}'.format(
                sample_id, first_event, second_event, first_event_start_time))

        # Test time between first event end and second event start. If it
        # is smaller than <max_delay>, store start of second event as result.
        # The first event end time being greater than the second event start
        # time, is explicitly allowed.
        if (second_event_start_time - first_event_end_time) <= max_delay:
            results.append({
                'sample_id': sample_id,
                'first_event_start': first_event_start_time,
                'first_event_end': first_event_end_time,
                'second_event_start': second_event_start_time
            })
        
        # Reset behavior tracking variables to find new pattern match.
        first_event_start_time = None
        first_event_end_time = None
        second_event_start_time = None
        second_event_end_time = None
            
    return results

# Open single samples 
'''
behavior_transitions = [
    #PostBehaviorTransition('17-08-26L1-cl', 'turn', 'bw', 11),
    #PostBehaviorTransition('17-08-26L2-cl', 'stim', 'fw', 3),
    #PostBehaviorTransition('17-08-26L5-cl', 'stim', 'fw', 3),
    PostBehaviorTransition('17-08-26L6-cl', 'turn', 'bw', 3)
]
'''
# Open all samples (!see CellConfig!) >ToDo
behavior_transitions = [
    PostBehaviorTransition(name,'stim', 'fw', 1) for name in lm_data]

found_transitions = []
for bt in behavior_transitions:
    sample_df = sample_data.get(bt.sample_id)
    if sample_df is None:
        raise ValueError('No data found for sample {}'.format(bt.sample_id))
    transitions = find_behavior_before(bt.sample_id, sample_df, bt.event, bt.post_event, bt.max_delay)
    
    if transitions:
        found_transitions.append(transitions)


print(len(found_transitions)) #number of data sets not the actual stim
print(sum([len(sample_transitions) for sample_transitions in found_transitions]))
#print(len(transitions)) 
print(found_transitions)


#%%
'''
# Test
def test_find_behavior_before():
    data_columns = ['time', 'bw_start', 'bw_end', 'bw_overlap', 'fw_start', 'fw_end', 'fw_overlap', 'turn_start', 'turn_end', 'turn_overlap']
    data = [
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [4, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [5, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    ]
    toy_df = pd.DataFrame(data, columns = data_columns)
    
    behavior_transitions = [
        #PostBehaviorTransition('17-08-26L1-cl', 'turn', 'bw', 11),
        #PostBehaviorTransition('17-08-26L2-cl', 'stim', 'fw', 3),
        #PostBehaviorTransition('17-08-26L5-cl', 'stim', 'fw', 3),
        PostBehaviorTransition('na', 'fw', 'bw', 5)
    ]
    
    found_transitions = []
    for bt in behavior_transitions:
        sample_df = toy_df
        if sample_df is None:
            raise ValueError('No data found for sample {}'.format(bt.sample_id))
        transitions = find_behavior_before(bt.sample_id, sample_df, bt.event, bt.post_event, bt.max_delay)

        if transitions:
            found_transitions.append(transitions)


    print(len(found_transitions)) #number of data sets not the actual stim
    print(len(transitions)) 
    print(found_transitions)
    
test_find_behavior_before()
'''


#%%
# For post_behavior_transition events get event_time and use ceLL_trace_config to filter by celltype and pattern.
# The results are merged and interpolated.
cell_trans_configs = []
all_trans_events = []

for sample in found_transitions:
    sample_ls_trans = []
    for found_transition in sample:
        #print(found_transition["sample_id"], found_transition["second_event_start"])
        #print(found_transition["sample_id"])
        # For all behavior except stimulus
        #sample_ls_trans.append(found_transition["second_event_start"]) 
        #cell_trans_configs.append(CellTransConfig(found_transition["sample_id"], "handle", 
        #                                          found_transition["second_event_start"]))
        
        # For stimulus as first_event
        sample_ls_trans.append(found_transition["first_event_start"]) 
        cell_trans_configs.append(CellTransConfig(found_transition["sample_id"], "A00c", 
                                                  found_transition["first_event_start"]))
        
                        
    #Find row of found_transition['second_event_start'] in sample_df
    #sample_df = sample_data.get(found_transition["sample_id"])
    #trans_df = sample_df.loc[sample_df['time'].isin(sample_ls_trans)]
    #trans_data[found_transition["sample_id"]] = trans_df
    #print(trans_df, found_transition["sample_id"])

# Extract for specific time window and align several events. 
# Define timepoints pre and post an event (event_df). 
# This works for single sample or multiple samples aligned 
# Note: In cell_subset_df, time was set to index, because for the previous avg calculation 
# Add index and time = column


# Set the window range left and right from the event
left_half_window_size = 100.0 #in seconds
right_half_window_size = 200.0

# trans_df defined in pargraph before 
windows = []
n_behavior = 0

for ctc in cell_trans_configs:
    sample_df = sample_data.get(ctc.sample_id)
    if sample_df is None:
        raise ValueError('{}: could not find sample data'.format(ctc.sample_id))
        continue    
    # Extract columns matching our cell type and the optional filter pattern.
    # Pandas' filter() operations works on columns for DataFrames by default.
    cell_subset_df = sample_df.filter(regex=ctc.get_filter_regex()) #Get subset of cells 
    cell_subset_df.set_index(sample_df.time, inplace=True) #Set time to index (essential for min/max...)
    cell_subset_df.reset_index(inplace = True) # Add index and time = column
    #print(cell_subset_df)
    
    n_behavior += 1
    window_start = ctc.event_time - left_half_window_size
    window_end = ctc.event_time + right_half_window_size
        
    # Get subset of rows between window_start and window_end       
    trans = cell_subset_df[(cell_subset_df.time >= window_start) & (cell_subset_df.time <= window_end)]
    #print(trans) #ok
    # Normalizing the data to align on beginning of selected
    # behavior (event_df = Zero) by substracting events in window
    # around start of event of interest from start of event interest.
    # Note: using ":" in event.loc[] will select "all rows" in our window.
    #trans.loc[:, 'time'] = trans['time'] - row['time']
    trans.loc[:, 'time'] = trans['time'] - ctc.event_time
    #print(trans) #ok
    # Add sample_id to each column as prefix and n_behavior as suffix to distinguish events within a sample
    #trans.rename(lambda x: '{}_{}_{}'.format(ctc.sample_id, x, n_behavior), axis = 'columns', inplace = True) 

    # Rename time collum to time
    trans.rename(columns={ trans.columns[0]: 'time' }, inplace = True) 
    #print(trans) # ok
    all_trans_events.append(trans) # Append a list with all event
#print(all_trans_events)   
        
        
# Removes first event and takes it as left_window in pd.merge_ordered and iterates than through all_events
all_trans_df = all_trans_events.pop(0)
for right_df in all_trans_events:
    all_trans_df = pd.merge_ordered(all_trans_df, right_df, on="time", how="outer")

# Resets the index as time and drops time column (sollte spaeter kommen)
all_trans_df.index = all_trans_df["time"]
del all_trans_df["time"]        
#print(all_trans_df)

#NR
# Index intepolation (linear interpolatione not on all_df, because index [=time] is not eaqually distributed)
int_all_Ptrans_df = all_trans_df.interpolate(method='index', axis=0, limit=None, inplace=False, limit_direction='both')
#print(int_all_trans_df)


#%%
# Average and stddev, min, max, sem for post_behavior_transition events
all_Ptrans_avg_df = int_all_Ptrans_df.mean(axis=1) # Interpolated data used
all_Ptrans_min_df = int_all_Ptrans_df.min(axis=1)
all_Ptrans_max_df = int_all_Ptrans_df.max(axis=1)
# Standard deviation (distribution)
all_Ptrans_std_df = int_all_Ptrans_df.std(axis = 1)
#standard error of mean
all_Ptrans_sem_df = int_all_Ptrans_df.sem(axis = 1)
#wrong zur haelfte: Want to have avg per celltyp over time point, 
#and not avg over all cells per timepoint

# Plotting for multi-events (same_behavioral_transition)
# If a dataframe with NANs is plotted (raw-data = non interpolated), use 
# marker = '+', or 'o', since the line in the lineplot only connects 
# consecutive data points
def aligned_layout_plot(plot, tick_spacing=1, fov=(-100, 200, 0.0, 1.0), legend=False): 
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False) 
    plot.axis(fov)
    plot.legend().set_visible(legend)

fig = plt.figure()

# Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
#sub1 = fig.add_subplot(211)
#all_trans_df.plot(ax=sub1, marker = '*', label = ctc.cell_type)
#aligned_layout_plot(sub1)

sub2 = fig.add_subplot(111) #212
all_Ptrans_avg_df.plot(ax=sub2, color = 'c', label = ctc.cell_type) #use interpolated df to calculate average...
#all_Ptrans_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
#all_Ptrans_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
#all_Ptrans_avg_df.plot.line(yerr=all_Ptrans_std_df, ax=sub2, color = 'r', alpha = 0.1)
all_Ptrans_avg_df.plot.line(yerr=all_Ptrans_sem_df, ax=sub2, color = 'grey', alpha = 0.1)
aligned_layout_plot(sub2)


#%%
'''
class BehaviorTransition:
    
    def __init__(self, sample_id, pre_event, event, post_event,
                 pre_max_delay=0, post_max_delay=0):
        self.sample_id = sample_id
        self.pre_event = pre_event
        self.event = event
        self.post_event = post_event
        self.pre_max_delay = pre_max_delay
        self.post_max_delay = post_max_delay

def find_behavior_between(sample_id, sample_df, first_event, second_event,
                          third_event, pre_max_delay=0, post_max_delay=0):
    """For the data frame of a single sample <df>, find all behaviors
    of type <second_event> that a) follows the event <first_event>,
    separated by <pre_max_delay> time. The start of <first_event> is expected
    to happen strictly before the start of <second_event>. The end time
    of <first_event> however can overlap with the start time of <second_event>.
    In this case, the time difference is negative, and still smaller than
    <pre_max_delay>. The end time of <first_event> can be before, at or after the
    end of <second_event>. And b) the behavior <second_event> is followed by the
    event <third_event>, separated by <post_max_delay> time. The end of
    <third_event> is expected to happen strictly after the end of <second_event>.
    The start time of <third_event> however can overlap with the end time of
    <second_event>. In this case, the time difference is negative, and still
    smaller than <post_max_delay>. The start time of <third_event> can be before,
    at or after the end of <second_event>.
    
    The start of <first_event> is expected to happen strictly before the start
    of <third_event>. Apart from this, both <first_event> and <third_event> can
    overlap.
    
    TODO: handle overlaps when behavior types are the same.
    """
    results = []
    first_event_start_col = '{}_start'.format(first_event)
    first_event_end_col = '{}_end'.format(first_event)
    first_event_overlap_col = '{}_overlap'.format(first_event)
    second_event_start_col = '{}_start'.format(second_event)
    second_event_end_col = '{}_end'.format(second_event)
    third_event_start_col = '{}_start'.format(third_event)
    third_event_end_col = '{}_end'.format(third_event)
    third_event_overlap_col = '{}_overlap'.format(third_event)
    
    first_event_start_time = None
    first_event_end_time = None
    second_event_start_time = None
    second_event_end_time = None
    third_event_start_time = None
    third_event_end_time = None

    for i, row in sample_df.iterrows():
        # Look for behaviors and remember its time.
        if row[first_event_start_col] and first_event_start_time is None:
            first_event_start_time = row['time']
            first_event_end_time = None
            continue
        if row[first_event_end_col] and not row[first_event_overlap_col]:
            first_event_end_time = row['time']
        if row[second_event_start_col]:
            second_event_start_time = row['time']
            second_event_end_time = None
        if row[second_event_end_col]:
            second_event_end_time = row['time']
        if row[third_event_start_col] and not row[third_event_overlap_col]:
            third_event_start_time = row['time']
   
        # As long as we haven't collected all needed time points,
        # keep on searching.
        if None in (first_event_end_time, second_event_start_time,
                    second_event_end_time, third_event_start_time):
            continue
            
        #NR
        if first_event_start_time > second_event_start_time:
            continue
        if first_event_start_time > first_event_end_time:
            continue
        if second_event_start_time > second_event_end_time:
            continue
        if second_event_start_time > third_event_start_time:
            continue
            
        if abs(first_event_start_time - second_event_start_time) < 0.00001:
            print('{}: start time (first) event {} and start time of (second) event {} are the same: {}'.format(
                sample_id, first_event, second_event, first_event_start_time))
        if abs(second_event_start_time - third_event_start_time) < 0.00001:
            print('{}: start time (second) event {} and start time of (third) event {} are the same: {}'.format(
                sample_id, second_event, third_event, second_event_start_time))

        #print(first_event_start_time, first_event_end_time, second_event_start_time, second_event_end_time, third_event_start_time, third_event_end_time)
        # Test time between first event end and second event start. If it
        # is smaller than <max_delay>, store start of second event as result.
        # The first event end time being greater than the second event start
        # time, is explicitly allowed.
        if (second_event_start_time - first_event_end_time) <= pre_max_delay \
                and (third_event_start_time - second_event_end_time) <= post_max_delay:
            results.append({
                'sample_id': sample_id,
                'first_event_start': first_event_start_time,
                'first_event_end': first_event_end_time,
                'second_event_start': second_event_start_time,
                'second_event_end': second_event_end_time,
                'third_event_start': third_event_start_time
            })
        
        # Reset behavior tracking variables to find new pattern match.
        first_event_start_time = None
        first_event_end_time = None
        second_event_start_time = None
        second_event_end_time = None
        third_event_start_time = None
        third_event_end_time = None

    return results

# Open single sample (!see CellConfig!) >ToDo
behavior_transitions = [
    BehaviorTransition('17-08-26L1-cl', 'bw', 'turn', 'bw', 10, 10),
]


# Open all samples (!see CellConfig!) >ToDo
behavior_transitions = [
    BehaviorTransition(name, 'bw', 'turn', 'bw', 3, 3) for name in lm_data]


found_transitions = []
for bt in behavior_transitions:
    sample_df = sample_data.get(bt.sample_id)
    if sample_df is None:
        raise ValueError('No data found for sample {}'.format(bt.sample_id))
    transitions = find_behavior_between(bt.sample_id, sample_df, bt.pre_event,
                                        bt.event, bt.post_event, bt.pre_max_delay,
                                        bt.post_max_delay)
    if transitions:
        found_transitions.append(transitions)

print(len(found_transitions))
print(found_transitions)
'''


#%%
#Will

###Not working

# Open single samples 
'''
first_transitions = [
    #PostBehaviorTransition('17-08-26L1-cl', 'turn', 'bw', 11),
    #PostBehaviorTransition('17-08-26L2-cl', 'stim', 'fw', 3),
    #PostBehaviorTransition('17-08-26L5-cl', 'stim', 'fw', 3),
    PostBehaviorTransition('17-08-26L6-cl', 'turn', 'bw', 3)
]
second_transitions = [
    #PostBehaviorTransition('17-08-26L1-cl', 'turn', 'bw', 11),
    #PostBehaviorTransition('17-08-26L2-cl', 'stim', 'fw', 3),
    #PostBehaviorTransition('17-08-26L5-cl', 'stim', 'fw', 3),
    PostBehaviorTransition('17-08-26L6-cl', 'bw', 'turn', 3)
]
'''
# Open all samples
first_transitions = [
    PostBehaviorTransition(name,'bw', 'turn', 3) for name in lm_data]
    
second_transitions = [
    PostBehaviorTransition(name,'turn', 'bw', 3) for name in lm_data]    

found_transitions = []
for first_bt, second_bt in zip(first_transitions, second_transitions):
    transitions=[]
    assert first_bt.sample_id == second_bt.sample_id, "{} does not match {}".format(first_bt.sample_id, second_bt.sample_id)
    sample_df = sample_data.get(first_bt.sample_id)
    if sample_df is None:
        raise ValueError('No data found for sample {}'.format(bt.sample_id))
    first_transitions = find_behavior_before(first_bt.sample_id, sample_df, first_bt.event, first_bt.post_event, first_bt.max_delay)
    second_transitions = find_behavior_before(second_bt.sample_id, sample_df, second_bt.event, second_bt.post_event, second_bt.max_delay)
    
    print("{} transitions from {} to {}".format(len(first_transitions), first_bt.event, first_bt.post_event))
    print("{} transitions from {} to {}".format(len(second_transitions), second_bt.event, second_bt.post_event))
    
    for ft in first_transitions:
        for st in second_transitions:
            if abs(ft["second_event_start"] - st["first_event_start"]) < 0.00001:
                transitions.append({
                    "sample_id":ft["sample_id"], "first_event_start":ft["first_event_start"], "first_event_end":ft["first_event_end"],
                    "second_event_start": st["first_event_start"], "second_event_end": st["first_event_end"],
                    "third_event_start": st["second_event_start"]
                })
    if transitions:
        print("{} transition triples found".format(len(transitions)))
        found_transitions.append(transitions)
    


print(len(found_transitions)) #number of data sets not the actual stim
print(sum([len(sample_transitions) for sample_transitions in found_transitions]))
print(found_transitions)


#%%
# All behavior_transitions, which were considered before assume that subsequent behaviors 
# are not the same. For same pairwise (2 behavior), we need to access also the data, when 
# the behaviors are same

class SamePairBehaviorTransition:
    
    def __init__(self, sample_id, pre_event, event, max_delay=0):
        self.sample_id = sample_id
        self.pre_event = pre_event
        self.event = event
        self.max_delay = max_delay

def find_behavior_next(sample_id, sample_df, first_event, second_event, max_delay=0):
    """For the data frame of a single sample <df>, find all behaviors
    of type <first_event> that will be followed by the same event <second_event>,
    separated by <max_delay> time. The start of <first_event> is expected
    to happen strictly before the start of <second_event>. The end time
    of <first_event> however can overlap with the start time of <second_event>.
    In this case, the time difference is negative, and still smaller than
    <max_delay>. The end time of <first_event> can be before, at or after the
    end of <second_event>.
    
    If <first_event> and <second_event> are the same type of behavior,
    overlaps have to be taken into account differently to match start and end times
    to the correct event. During iteration for one loop, we have to exclude the 
    fact that the first_event == second_event.
    """
    
    print("finding same behaviors only")
    
    results = []
    first_event_start_col = '{}_start'.format(first_event)
    first_event_end_col = '{}_end'.format(first_event)
    first_event_overlap_col = '{}_overlap'.format(first_event)
    second_event_start_col = '{}_start'.format(second_event)
    second_event_end_col = '{}_end'.format(second_event)
    second_event_overlap_col = '{}_overlap'.format(second_event) 
    
    first_event_start_time = None
    first_event_end_time = None
    second_event_start_time = None
    second_event_end_time = None
    
    # Check for overlap between the same behaviors (print index, where 'True') and use
    # it as a check that there is not this error in the behavior data
    #print(sample_id, sample_df.index[sample_df['bw_overlap']].tolist())
    
    # Note: The overlap statement was removed. This part has to be 
    # checked if overlapping events are found in the data
    for i, row in sample_df.iterrows():
        # Look for start of first behavior and remember its time.
        if row[first_event_start_col]and first_event_start_time is None:
            first_event_start_time = row['time']
        if row[first_event_end_col] and first_event_end_time is None: 
            first_event_end_time = row['time']
        if row[second_event_start_col] and first_event_start_time is not None:
            second_event_start_time = row['time']
        if row[second_event_end_col] and first_event_end_time is not None:
            second_event_end_time = row['time']
        for column in sample_df.columns:
            if (first_event_start_time is not None and
                column.endswith("_start") and
                column != first_event_start_col and
                column != second_event_start_col and
                first_event not in column and
                second_event not in column):
                if row[column]:
                    #print("{} ended at {}, but then found {} at {}".format(first_event, first_event_end_time, column, row["time"]))
                    first_event_start_time = None
                    first_event_end_time = None
                    second_event_start_time = None
                    second_event_end_time = None
        # As long as we haven't collected all needed time points,
        # keep on searching.
        if None in (first_event_start_time, first_event_end_time,
                    second_event_start_time, second_event_end_time):
            continue
        
        #NR
        if first_event_start_time == second_event_start_time:
            continue
        if first_event_end_time == second_event_end_time:
            continue
        if first_event_start_time > first_event_end_time:
            continue
        if first_event_start_time > second_event_start_time:
            continue

        
        # Test time between first event end and second event start. If it
        # is smaller than <max_delay>, store start of second event as result.
        # The first event end time being greater than the second event start
        # time, is explicitly allowed. During iteration the first_event == second_event. 
        if (second_event_start_time == first_event_end_time): #NR
            continue
        if (second_event_start_time - first_event_end_time) <= max_delay:
            results.append({
                'sample_id': sample_id,
                'first_event_start': first_event_start_time,
                'first_event_end': first_event_end_time,
                'second_event_start': second_event_start_time,
                'second_event_end': second_event_end_time
            })
        
        # Reset behavior tracking variables to find new pattern match.
        first_event_start_time = second_event_start_time
        first_event_end_time = second_event_end_time
        second_event_start_time = None
        second_event_end_time = None
            
    return results


# Open single sample (!see CellConfig!) >ToDo

behavior_transitions = [
    #SamePairBehaviorTransition('17-08-26L1-cl', 'turn', 'turn', 2),
    SamePairBehaviorTransition('17-08-31L1-cl', 'bw', 'bw', 3),
    SamePairBehaviorTransition('17-11-06L3-cl', 'bw', 'bw', 3),
    SamePairBehaviorTransition('17-11-29L1-cl', 'bw', 'bw', 3)
]

'''
# Open all samples 
behavior_transitions = [
    SamePairBehaviorTransition(name, 'fw', 'fw', 3) for name in lm_data]
'''

found_transitions = []
for bt in behavior_transitions:
    sample_df = sample_data.get(bt.sample_id)
    if sample_df is None:
        raise ValueError('No data found for sample {}'.format(bt.sample_id))
    transitions = find_behavior_next(bt.sample_id, sample_df, bt.pre_event, bt.event, bt.max_delay)

    if transitions:
        found_transitions.append(transitions)

#print(len(transitions))
print(len(found_transitions))
print(sum([len(sample_transitions) for sample_transitions in found_transitions]))
print(found_transitions)


#%%
print(sum([len(sample_transitions) for sample_transitions in found_transitions]))


#%%
'''
#test
def test_find_behavior_next():
    data_columns = ['time', 'bw_start', 'bw_end', 'bw_overlap', 'fw_start', 'fw_end', 'fw_overlap', 'turn_start', 'turn_end', 'turn_overlap']
    data = [
        [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [1, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [2, 0, 0, 0, 1, 0, 0, 0, 0, 0],
        [3, 0, 0, 0, 0, 1, 0, 0, 0, 0],
        [4, 0, 0, 0, 0, 0, 0, 1, 0, 0],
        [5, 0, 0, 0, 0, 0, 0, 0, 1, 0],
        [6, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [7, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [8, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [9, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [10, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [11, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [12, 1, 0, 0, 0, 0, 0, 0, 0, 0],
        [13, 0, 1, 0, 0, 0, 0, 0, 0, 0],
        [14, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    ]
    toy_df = pd.DataFrame(data, columns = data_columns)

    behavior_transitions = [
    SamePairBehaviorTransition('na', 'bw', 'bw', 12)]


    found_transitions = []
    for bt in behavior_transitions:
        sample_df = toy_df
        if sample_df is None:
            raise ValueError('No data found for sample {}'.format(bt.sample_id))
        transitions = find_behavior_next(bt.sample_id, sample_df, bt.pre_event, bt.event, bt.max_delay)

        if transitions:
            found_transitions.append(transitions)

    print(len(transitions))
    print(len(found_transitions))
    print(found_transitions)
    
test_find_behavior_next()
'''


#%%
# For same_behavior_transition events get event_time and use cee_trace_config to filter by celltype and pattern.
# The results are merged and interpolated.
cell_Strans_configs = []
all_Strans_events = []

for sample in found_transitions:
    sample_ls_trans = []
    for found_transition in sample:
        #print(found_transition["sample_id"], found_transition["second_event_start"])
        #print(found_transition["sample_id"])
        sample_ls_trans.append(found_transition["second_event_start"])
        cell_Strans_configs.append(CellTransConfig(found_transition["sample_id"], "handle", 
                                                  found_transition["second_event_start"]))

# Extract for specific time window and align several events. 
# Define timepoints pre and post an event (event_df). 
# This works for single sample or multiple samples aligned 
# Note: In cell_subset_df, time was set to index, because for the previous avg calculation 
# Add index and time = column

# Set the window range left and right from the event
left_half_window_size = 200.0 #in seconds
right_half_window_size = 300.0

# trans_df defined in pargraph before 
windows = []
n_behavior = 0

for ctc in cell_Strans_configs:
    sample_df = sample_data.get(ctc.sample_id)
    if sample_df is None:
        raise ValueError('{}: could not find sample data'.format(ctc.sample_id))
        continue    
    # Extract columns matching our cell type and the optional filter pattern.
    # Pandas' filter() operations works on columns for DataFrames by default.
    cell_subset_df = sample_df.filter(regex=ctc.get_filter_regex()) #Get subset of cells 
    cell_subset_df.set_index(sample_df.time, inplace=True) #Set time to index (essential for min/max...)
    cell_subset_df.reset_index(inplace = True) # Add index and time = column
    #print(cell_subset_df)
    
    n_behavior += 1
    window_start = ctc.event_time - left_half_window_size
    window_end = ctc.event_time + right_half_window_size
        
    # Get subset of rows between window_start and window_end       
    trans = cell_subset_df[(cell_subset_df.time >= window_start) & (cell_subset_df.time <= window_end)]
    #print(trans) #ok
    # Normalizing the data to align on beginning of selected
    # behavior (event_df = Zero) by substracting events in window
    # around start of event of interest from start of event interest.
    # Note: using ":" in event.loc[] will select "all rows" in our window.
    #trans.loc[:, 'time'] = trans['time'] - row['time']
    trans.loc[:, 'time'] = trans['time'] - ctc.event_time
    #print(trans) #ok
    # Add sample_id to each column as prefix and n_behavior as suffix to distinguish events within a sample
    #trans.rename(lambda x: '{}_{}_{}'.format(ctc.sample_id, x, n_behavior), axis = 'columns', inplace = True) 

    # Rename time collum to time
    trans.rename(columns={ trans.columns[0]: 'time' }, inplace = True) 
    #print(trans) # ok
    all_Strans_events.append(trans) # Append a list with all event
#print(all_trans_events)   
        
        
# Removes first event and takes it as left_window in pd.merge_ordered and iterates than through all_events
all_trans_df = all_Strans_events.pop(0)
for right_df in all_Strans_events:
    all_trans_df = pd.merge_ordered(all_trans_df, right_df, on="time", how="outer")

# Resets the index as time and drops time column (sollte spaeter kommen)
all_trans_df.index = all_trans_df["time"]
del all_trans_df["time"]        
#print(all_trans_df)

# Index intepolation (linear interpolatione not on all_df, because index [=time] is not eaqually distributed)
int_all_Strans_df = all_trans_df.interpolate(method='index', axis=0, limit=None, inplace=False, limit_direction='both')
#print(int_all_Strans_df)


#%%
# Average and stddev, min, max, sem for same_behavior_transition events
all_Strans_avg_df = int_all_Strans_df.mean(axis=1) # Interpolated data used
all_Strans_min_df = int_all_Strans_df.min(axis=1)
all_Strans_max_df = int_all_Strans_df.max(axis=1)
# Standard deviation (distribution)
all_Strans_std_df = int_all_Strans_df.std(axis = 1)
#standard error of mean
all_Strans_sem_df = int_all_Strans_df.sem(axis = 1)
#wrong zur haelfte: Want to have avg per celltyp over time point, 
#and not avg over all cells per timepoint

# Plotting for multi-events (same_behavioral_transition)
# If a dataframe with NANs is plotted (raw-data = non interpolated), use 
# marker = '+', or 'o', since the line in the lineplot only connects 
# consecutive data points
def aligned_layout_plot(plot, tick_spacing=1, fov=(-10, 10, 0.1, 0.4), legend=False): 
    # Set fine x-axis scale
    plot.xaxis.set_major_locator(ticker.MultipleLocator(tick_spacing))

    # Set x and y limits and legend (default = False) 
    plot.axis(fov)
    plot.legend().set_visible(legend)

fig = plt.figure()

# Plot all cells from all_df, aligned at zero for event_start, specified in Cell_Trace_Config.
#sub1 = fig.add_subplot(111) #211
#all_trans_df.plot(ax=sub1, marker = '*', label = ctc.cell_type)
#aligned_layout_plot(sub1)

sub2 = fig.add_subplot(111) #212
all_Strans_avg_df.plot(ax=sub2, color = 'g', label = ctc.cell_type) #use interpolated df to calculate average...
#all_Strans_min_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
#all_Strans_max_df.plot(ax=sub2, color = 'r', linewidth=1, alpha = 0.5)
#all_Strans_avg_df.plot.line(yerr=all_Strans_std_df, ax=sub2, color = 'r', alpha = 0.1)
all_Strans_avg_df.plot.line(yerr=all_Strans_sem_df, ax=sub2, color = 'grey', alpha = 0.1)
aligned_layout_plot(sub2)

#%% [markdown]
# Overall analysis

#%%



#%%
# Correlation matrix (using merged_ordered + interpolated) from behavioral transitions
# HERE: Not aneraged!! --> int_all_Strans_df (aligned to sec_event_start)

# Before 


#%%



#%%



#%%



#%%


#%% [markdown]
# Testing 

#%%
#print(list(event_df)) #prints header
#rounded time is only visual, I still get several 'same' tp
#sample_df.round({'time' : 1})
   # Round time on 1 or 2nd decimal
    # the df.round({'time' : 1}) doesn't work if to many decimals
    #decimals = 2    
    #timestamp_df['time'] = timestamp_df['time'].apply(lambda x: round(x, decimals))


#%%



#%%
d = {'A' : [3.2, np.nan, np.nan, 5, np.nan, np.nan],
     'B' : [np.nan, 4.1, np.nan, np.nan, 6.2, np.nan], 
     'C' : [np.nan, np.nan, 1.1, np.nan, np.nan, 2.5]}
df = pd.DataFrame(data=d)
print(df)


int_df_linear = df.interpolate(method = 'linear', inplace = False)
print(int_df_linear)

int_df_index = df.interpolate(method='index', inplace=False)
print(int_df_index)
#u = df['T'].unique

# Resets the index as time and drops time column
#df.index = df["T"]
#del df["T"] 

#df.fillna(value=None, method='ffill', axis=None, inplace=False, limit=None, downcast=None) # not what I want
#df.interpolate(method='linear', axis=0, limit=None,
#                      inplace=False, limit_direction='forward', limit_area=None, downcast=None) 


# print out only rows where in'T' are duplicates
#df[df.duplicated(subset= ['T'], keep=False)]

    
    
#np.nanmean(j, axis=0)
#    for j in df.loc[df["T"]== value]:
#        np.nanmean(j, axis=0)


#df.mean(axis=1, skipna=None)
#print(df)


#%%
t = {'T' : [1,2, 3, 4, 5, 6, 7, 8], 'A' : [3.2, 5, 5.5, 5.3, 9, 8, 8, 3],
     'B' : [4.1, 6.2, 6.0, 6.2, 8, 1, 1.5, 3.7], 
     'C' : [1.1, 2.5, 2.3, 1.2, 0.9, 1.1, 1.8, 1.7]}
df1 = pd.DataFrame(data=d)
df1

d = {'A' : [3.2, 5, 5.5, 5.3, 9, 8, 8, 3],
     'B' : [4.1, 6.2, 6.0, 6.2, 8, 1, 1.5, 3.7], 
     'C' : [1.1, 2.5, 2.3, 1.2, 0.9, 1.1, 1.8, 1.7]}
df2 = pd.DataFrame(data=d)
print(df2)


#%%
# Interpolation 
# inplace  = False, since we want to keep the data sets with raw data as well

# Linear Interpolation: According to documentation, 
# because it assums the index is equally spaced.
# ‘Index’, ‘values’: use the actual numerical values of the index.

int_df2 = df2.interpolate(method = 'index', inplace = False)
print(int_df2)

#Note: First 5 values = NAN!!??!! (Method?)

#intpol_all_df = all_df.interpolate(method='index', inplace=False)
#print(intpol_all_df)
ls = [2, 1.2, 3, 8.2]
#for i in ls:
    #print(i)
a =int_df2.loc[int_df2['C'].isin(ls)]
print(a)
# find row where C = 1.2 #isin
#a =int_df2.loc[df['C'].isin(1.2)]
#a =int_df2.loc[int_df2['C'] == 1.2]

#print(a)


#%%
# Data anlysis - TODO
# Ask Ann if this kind of interpolation made sens

# Dataprocessing
# For the next step, we try two methods two normalize the data and get the timestamps 
# in synchrony between the different samples/events (1.Interpolation of some kind, 2. Binning,
# 3) fitting curve)
               
# Interpolation 
# inplace  = False, since we want to keep the data sets with raw data as well

# Linear Interpolation: According to documentation it is not correct to use, 
# because it assums the index is equally spaced.

# ‘Index’, ‘values’: use the actual numerical values of the index.

#Note: First 5 values = NAN!!??!!

#intpol_all_df = all_df.interpolate(method='index', inplace=False)
#print(intpol_all_df)


