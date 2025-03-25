import pandas as pd
import os
import re
from collections import defaultdict
import math
import numpy as np
    
print("Utils module is being imported!")

#reads the output.log file and returns a dataframe
def readlogfile(logfile_path):
    log_file_full_path = os.path.join(logfile_path, 'output.log')

    with open(log_file_full_path) as file:
        lines = file.readlines()
        start_index = 0
        data_lines = []
        for i, line in enumerate(lines):
            if 'Step E_bond TotEng Temp' in line:
                start_index = i + 1
                break

        # Collect data until a line starts with " loop"
        for line in lines[start_index:]:
            if line.strip().startswith("Loop"):
                break
            data_lines.append(line.split())

    # Now read the actual data into a DataFrame
    df = pd.DataFrame(data_lines, columns=lines[start_index-1].split())

    # Convert necessary columns to float or int as appropriate
    df['Step'] = df['Step'].astype(int)
    for col in df.columns[1:]:  # Convert all other columns assuming they are numeric
        df[col] = pd.to_numeric(df[col], errors='coerce')  # Use 'coerce' to handle any non-numeric entries safely
    return df

# average the dataframes of each path in the list "paths"
def average_runs(paths):
    dfs = []
    for path in paths:
        df = readlogfile(path)
        df.set_index('Step', inplace=True)
        #df['substrate'] = 1000 - (df['f_bondc0[2]'] + df['f_bondc1[2]'] + df['f_bondc2[2]'] + df['f_bondc3[2]'] + df['f_bondc4[2]'] + df['f_bondc5[2]'] + df['f_bondc6[2]'])
        dfs.append(df)

    # Calculate mean and variance
    combined_df = pd.concat(dfs)
    mean_df = combined_df.groupby(level=0).mean()
    var_df = combined_df.groupby(level=0).var()

    return mean_df, var_df

# reads the clusterDist.dat file and returns as list of lists
def read_distribution(logfile_path):
    #read the distribution of cluster sizes
    file_path_dist = os.path.join(logfile_path, 'cluster_results_skip_first_frame/clusterDist.dat')
    with open(file_path_dist, 'r') as f:
        # Read lines, strip newline characters, split each line into a list, and convert each item to an integer
        distribution = [[int(item) for item in line.strip().split()] for line in f.readlines()]
    return distribution 

# computes the mean size (how many enzymes it contains) of clusters >= threshold and the variance of the sizes
def calculate_mean_size(logfile_path,threshold = 5):
    #calculate the mean size of clusters > threshold 
    distribution = read_distribution(logfile_path)
    counts = [sum(1 for item in sublist if item >= threshold) for sublist in distribution] #could also return # of clusters if needed.
    mean_sizes = []
    variance_sizes = []
    for sublist in distribution:
        filtered_elements = [item for item in sublist if item >= threshold]
        if filtered_elements:  # Avoid division by zero if no elements are > 5
            mean_value = sum(filtered_elements) / len(filtered_elements)
            var_value = sum((x - mean_value) ** 2 for x in filtered_elements) / len(filtered_elements)
        else:
            mean_value = 0 
            var_value = 0
        mean_sizes.append(mean_value)
        variance_sizes.append(var_value)
    return mean_sizes, variance_sizes, counts

#calculate the percentage out of all enzymes in simulation that are in clusters
def calculate_enzyme_ratio(logfile_path,threshold = 5):
    #calculate the ratio of enzymes in the cluster
    distribution = read_distribution(logfile_path)
    enzyme_ratio_in_cluster = [sum(item for item in sublist if item >= threshold)/sum(sublist) for sublist in distribution] 
    return enzyme_ratio_in_cluster

#returns a dataframe containing cluster size max, mean, variance, enzyme ratio and step
def readcluster(logfile_path, threshold = 5):
    #read the cluster mean and max size as function of time 
    # Construct the file path
    file_path_max = os.path.join(logfile_path, 'cluster_results_skip_first_frame/clusterSizeMax.dat')
    
    with open(file_path_max, 'r') as file:
        lines = file.readlines()
        data_lines = [line.split() for line in lines]

    # Create DataFrame from read data
    df = pd.DataFrame(data_lines,columns=['clusterSizeMax']) 
    df['clusterSizeMax'] = df['clusterSizeMax'].astype(int)
    df['clusterSizeMean'], _, df['clusterNumber'] = calculate_mean_size(logfile_path, threshold)
    df['enzymeRatio'] = calculate_enzyme_ratio(logfile_path, threshold)
    df['Step'] = df.index

    return df



def read_complete_clusters(paths):
    column_names = ['ClusterNumber', 'ClusterSize', 'MeanClusterSize']
    df = pd.DataFrame(columns=column_names)

    total_cluster_numbers = []
    total_cluster_sizes = []
    total_cluster_mean_sizes = []

    for path in paths:
    # Paths to the files
        path_number = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_number.dat')
        path_sizes = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_sizes.dat')
        path_mean_sizes = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_mean_sizes.dat')

        # Initialize lists to hold file data
        cluster_numbers = []
        cluster_sizes = []
        cluster_mean_sizes = []

        # Read cluster numbers
        with open(path_number, 'r') as file:
            for line in file:
                cluster_numbers.append(int(line.strip()))

        # Read cluster sizes
        with open(path_sizes, 'r') as file:
            for line in file:
                numbers = line.split()  # This splits the string into a list of numbers as strings
                cluster_sizes.append([int(number) for number in numbers])

        with open(path_mean_sizes, 'r') as file:
            for line in file:
                mean_size = float(line.strip())
                cluster_mean_sizes.append(mean_size)

        total_cluster_numbers.append(cluster_numbers)
        total_cluster_sizes.append(cluster_sizes)
        total_cluster_mean_sizes.append(cluster_mean_sizes)


    # Create a DataFrame
    df = pd.DataFrame({
        'ClusterNumber': total_cluster_numbers,
        'ClusterSize': total_cluster_sizes,
        'MeanClusterSize' : total_cluster_mean_sizes
    })

    return df


def read_complete_clusters_one_run(path):
    column_names = ['ClusterNumber', 'ClusterSize', 'MeanClusterSize']
    df = pd.DataFrame(columns=column_names)

    path_number = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_number.dat')
    path_sizes = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_sizes.dat')
    path_mean_sizes = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_mean_sizes.dat')

    # Initialize lists to hold file data
    cluster_numbers = []
    cluster_sizes = []
    cluster_mean_sizes = []

    # Read cluster numbers
    with open(path_number, 'r') as file:
        for line in file:
            cluster_numbers.append(int(line.strip()))

    # Read cluster sizes
    with open(path_sizes, 'r') as file:
        for line in file:
            numbers = line.split()  # This splits the string into a list of numbers as strings
            cluster_sizes.append([int(number) for number in numbers])

    with open(path_mean_sizes, 'r') as file:
        for line in file:
            mean_size = float(line.strip())
            cluster_mean_sizes.append(mean_size)

    # Create a DataFrame
    df = pd.DataFrame({
        'ClusterNumber': cluster_numbers,
        'ClusterSize': cluster_sizes,
        'MeanClusterSize' : cluster_mean_sizes
    })

    return df


def average_complete_clusters(df):
    list_df = pd.DataFrame(df['ClusterNumber'].tolist())
    # Calculate the mean of each column (element-wise mean across all lists)
    averages = list_df.mean()
    # Calculate variance and standard deviation
    std_dev_cluster_number = list_df.std()
    # Convert the averages, variance, and std deviation Series back into lists
    average_complete_cluster_number = averages.tolist()
    std_dev_cluster_number_list = std_dev_cluster_number.tolist()

    # For MeanClusterSize
    list_df2 = pd.DataFrame(df['MeanClusterSize'].tolist())
    averages2 = list_df2.mean()
    # Calculate variance and standard deviation
    std_dev_cluster_size = list_df2.std()
    # Convert the averages, variance, and std deviation Series back into lists
    average_complete_cluster_size = averages2.tolist()
    std_dev_cluster_size_list = std_dev_cluster_size.tolist()

    howmany = len(df['ClusterSize'][0])
    # Initialize a list of 101 empty lists
    
    combined_complete_cluster_sizes = [[] for _ in range(howmany)]
    
    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        # row['ClusterSize'] is a list of 101 sublists
        for i, sublist in enumerate(row['ClusterSize']):
            # Extend the ith sublist in combined_sublists with the current sublist
            combined_complete_cluster_sizes[i].extend(sublist)

    average_combined_complete_cluster_sizes = [np.mean(sublist) for sublist in combined_complete_cluster_sizes]
    
    return average_complete_cluster_number, std_dev_cluster_number_list, average_complete_cluster_size, std_dev_cluster_size_list
   


# calls readcluster for each path in the list "paths" and averages the dataframes
def average_clusters(paths, threshold=5):
    dfs = []
    for path in paths:
        df = readcluster(path, threshold)
        dfs.append(df)
    
    columns = dfs[0].columns
    data = np.stack([df.values for df in dfs], axis=2)

    # Compute mean and variance for each row across the 6 dataframes
    mean_values = np.mean(data, axis=2)
    variance_values = np.var(data, axis=2)

    # Create the final dataframe with mean and variance
    final_df = pd.DataFrame(np.hstack([mean_values, variance_values]),
                            columns=[f'{col}_mean' for col in columns] + [f'{col}_var' for col in columns])

    return final_df
    # dfs = []
    # for path in paths:
    #     df = readcluster(path, threshold)
    #     dfs.append(df)
    
    # concatenated_df = pd.concat(dfs, axis=0)
    # mean_df = concatenated_df.groupby(concatenated_df.index).mean()
    # variance_df = concatenated_df.groupby(concatenated_df.index).var()
    
    # # Combine mean and variance DataFrames
    # result_df = mean_df.copy()
    # result_df['clusterSizeMaxVariance'] = variance_df['clusterSizeMax']
    # result_df['enzymeRatioVariance'] = variance_df['enzymeRatio']
    # result_df['clusterNumberVariance'] = variance_df['clusterNumber']
    
    # return result_df


#mode atom will be 1 to 1. mode molecule will be 1 to 1 for now, containing only the enzyme name, not the active site name. In the future, will convert to a list of integers
def create_index_to_name_dictionary(system_data_path, mode):
    if mode not in ["molecule", "atom"]:
        raise ValueError("Invalid mode. Choose either 'molecule' or 'atom'.")
    
    output = {}
    
    with open(system_data_path, 'r') as file:
        keyword_found = False
        keyword = "Atoms"
        
        for line in file:
            if keyword_found:
                parts = line.split()
                if len(parts) < 7:
                    continue
                
                atom_index = int(parts[0])  
                molecule_index = int(parts[1])
                name = int(parts[2])
                
                if(mode == 'molecule' and molecule_index not in output):
                    output[molecule_index] = name
                if (mode == 'atom'):
                    output[atom_index] = name
                
            elif keyword in line:
                keyword_found = True
            else:
                continue
    
    return output
            
            
#mode atom will be 1 to 1. mode molecule will have valus as list of integers
def create_name_to_list_of_indices_dictionary(system_data_path, mode):
    if mode not in ["molecule", "atom"]:
        raise ValueError("Invalid mode. Choose either 'molecule' or 'atom'.")
    
    with open(system_data_path, 'r') as file:
        keyword_found = False
        keyword = "atom types"
        atom_type = 0
        
        for line in file:
            if keyword in line:
                parts = line.split()
                atom_type = int(parts[0])
                exit
                
    output = {new_list: [] for new_list in range(1, atom_type + 1)}    
        
    with open(system_data_path, 'r') as file:
        keyword_found = False
        keyword = "Atoms"
        
        for line in file:
            if keyword_found:
                parts = line.split()
                if len(parts) < 7:
                    continue
                
                atom_index = int(parts[0])  
                molecule_index = int(parts[1])
                name = int(parts[2])
                
                if(mode == 'molecule'):
                    output[name].append(molecule_index)
                if (mode == 'atom'):
                    output[name].append(atom_index)
                
            elif keyword in line:
                keyword_found = True
            else:
                continue
    
    return output
    
#reads the composition of the biggest cluster in all frames and returns as a list of lists
def readmaxclustercomposition(logfile_path):
    # Construct the file path
    file_path = os.path.join(logfile_path, 'cluster_results_skip_first_frame/clusterMaxList_corrected.dat')
    
    # Initialize an empty list to store the data
    data_lines = []
    
    system_data_path = logfile_path + "system.data"
    index_to_name = create_index_to_name_dictionary(system_data_path, "atom")
    
    # Open the file and read lines
    with open(file_path, 'r') as file:
        for line in file:
            # Split each line, convert to integers, perform integer division by 50, and append to the list
            processed_line = [index_to_name[int(num)] for num in line.split()]
            data_lines.append(processed_line)
    
    return data_lines

# reads the composition of all clusters in all frames and returns as a list of lists
# future work: filtering based on cluster size, compeleteness
# if complete_clusters_only = false, returns all clusters, converted to type
# if complete_clusters_only = true, returns only complete clusters, NOT converted to type but is in atom indices
def read_all_cluster_compositions(logfile_path, enz_arr, complete_clusters_only, threshold = 5):
    '''Reads all cluster compositions from the clusterList.dat file and returns a list of lists of integers.'''
    file_path = os.path.join(logfile_path, 'cluster_results_skip_first_frame/clusterList_corrected.dat')
    data_lines = []
    
    complete_set = set(enz_arr)
    system_data_path = os.path.join(logfile_path, "system.data")
    index_to_name = create_index_to_name_dictionary(system_data_path, "atom")

    with open(file_path, 'r') as file:
        # Process each line to extract patterns enclosed by curly braces
        for line in file:
            data_line = []
            # Find all substrings within curly braces and split each match into lists of integers
            matches = re.findall(r'\{(.*?)\}', line)
            
            # Process each match into lists of integers
            processed_lines = [list(map(int, match.split())) for match in matches if len(match.split()) >= threshold]
            
            # Apply the transformation to each integer in the processed lines
            transformed_lines = [[index_to_name[int(value)] for value in sublist] for sublist in processed_lines]
            
            # Append the transformed lists to data_line
            if complete_clusters_only:
                for index, sublist in enumerate(transformed_lines):
                    if set(complete_set).issubset(sublist):
                        #not zero index
                        zero_index_line = [x for x in processed_lines[index]]
                        data_line.append(zero_index_line)
                data_lines.append(data_line)
            else:
                data_lines.append(transformed_lines)
        

    return data_lines

#for all paths, reads the max cluster composition and combines them into one list of lists grouped by frames
def average_maxclustercomposition(paths):
    combined_data = None
    for path in paths:
        # Read the cluster data, which is a list of lists
        data = readmaxclustercomposition(path)
        
        if combined_data is None:
            # If combined_data is not initialized, set it as the first path's data
            combined_data = data
        else:
            # Merge each corresponding sublist from the current data with combined_data
            combined_data = [x + y for x, y in zip(combined_data, data)]

    return combined_data

# for all paths, reads the composition of all clusters and combines them into one list of lists grouped by frames
def average_all_cluster_compositions(paths, enz_arr, threshold = 5):
    combined_data = None
    for path in paths:
        # Read the cluster data, which is a list of lists
        data = read_all_cluster_compositions(path, enz_arr, False, threshold)
        if combined_data is None:
            combined_data = data
        else:
            # Merge each corresponding sublist from the current data with combined_data
            combined_data = [x + y for x, y in zip(combined_data, data)]

    return combined_data


def average_all_complete_cluster_compositions(paths, enz_arr, threshold = 5):
    combined_data = None
    for path in paths:
        # Read the cluster data, which is a list of lists
        data = read_all_cluster_compositions(path, enz_arr, True, threshold)
        if combined_data is None:
            combined_data = data
        else:
            # Merge each corresponding sublist from the current data with combined_data
            combined_data = [x + y for x, y in zip(combined_data, data)]

    return combined_data



def filter_incomplete_clusters_and_save_to_file_exclude_COQ9(paths, enz_arr = [1,2,3,4,5], threshold = 5):
    for path in paths:
        # Define paths for output files
        output_path = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_excluding_COQ9.dat')
        number_output_path = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_number_excluding_COQ9.dat')
        size_output_path = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_sizes_excluding_COQ9.dat')
        mean_size_output_path = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_mean_sizes_excluding_COQ9.dat')
        
        data = read_all_cluster_compositions(path, enz_arr, True, threshold)
        
        # Open both files: one for cluster data and another for sizes
        with open(output_path, 'w') as file, open(number_output_path, 'w') as number_file, open(size_output_path, 'w') as size_file, open(mean_size_output_path, 'w') as mean_size_file:
            for sublist in data:
                if sublist == []:
                    file.write('\n')
                    number_file.write('0\n')  # Write size 0 for empty sublists
                    size_file.write('0\n')
                    mean_size_file.write('0\n')
                else:
                    file.write(' '.join(['<' + ' '.join(map(str, inner_list)) + '>' for inner_list in sublist]))
                    file.write('\n')
                    number_file.write(f"{len(sublist)}\n")  # Write the size of each sublist
                    size_file.write(' '.join([str(len(inner_list)) for inner_list in sublist]))
                    size_file.write('\n')
                    mean_size_file.write(f"{sum(len(inner_list) for inner_list in sublist) / len(sublist)}\n")
                    
        print("Filtered for only complete cluster (1-indexing) file saved for path: " + output_path)
        print(f"Number of complete clusters by frame saved for path: {number_output_path}")
        print(f"Sizes of complete clusters by frame saved for path: {size_output_path}")


def filter_incomplete_clusters_and_save_to_file(paths, enz_arr, threshold = 5):
    for path in paths:
        # Define paths for output files
        output_path = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters.dat')
        number_output_path = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_number.dat')
        size_output_path = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_sizes.dat')
        mean_size_output_path = os.path.join(path, 'cluster_results_skip_first_frame/clusterList_only_complete_clusters_mean_sizes.dat')

        
        # if os.path.exists(output_path) and os.path.exists(number_output_path) and os.path.exists(size_output_path):
        #     print(f"Filtered files already exist in {path}.")
        #     continue
        
        data = read_all_cluster_compositions(path, enz_arr, True, threshold)
        
        # Open both files: one for cluster data and another for sizes
        with open(output_path, 'w') as file, open(number_output_path, 'w') as number_file, open(size_output_path, 'w') as size_file, open(mean_size_output_path, 'w') as mean_size_file:
            for sublist in data:
                if sublist == []:
                    file.write('\n')
                    number_file.write('0\n')  # Write size 0 for empty sublists
                    size_file.write('0\n')
                    mean_size_file.write('0\n')
                else:
                    file.write(' '.join(['<' + ' '.join(map(str, inner_list)) + '>' for inner_list in sublist]))
                    file.write('\n')
                    number_file.write(f"{len(sublist)}\n")  # Write the size of each sublist
                    size_file.write(' '.join([str(len(inner_list)) for inner_list in sublist]))
                    size_file.write('\n')
                    mean_size_file.write(f"{sum(len(inner_list) for inner_list in sublist) / len(sublist)}\n")
                    
        print("Filtered for only complete cluster (0-indexing restored) file saved for path: " + output_path)
        print(f"Number of complete clusters by frame saved for path: {number_output_path}")
        print(f"Sizes of complete clusters by frame saved for path: {size_output_path}")



# given a nested list (output of average_maxclustercomposition), finds the cluster composition averaeged over specified frames 
def max_cluster_process_data(nested_lists,group_names_int,mode,frames=10):
    # Select the last 10 lists from the nested_lists
    last_n = []
    if (mode == "end"):
        last_n = nested_lists[-1*frames:]
    elif (mode == "middle"):
        temp = len(nested_lists)
        mid = int(temp // 2)
        last_n = nested_lists[mid - frames // 2: mid + frames // 2 + 1]
    else:
        print("Invalid mode. Should be either middle or end.")
      
    # Flatten the last 5 lists into a single list
    flattened_data = [item for sublist in last_n for item in sublist]
    # Total number of items in the flattened data
    total_items = len(flattened_data)
    
    # Count occurrences of numbers 1 to 5 and calculate percentage
    percentages = {i: (flattened_data.count(i) / total_items) * 100 for i in group_names_int}
    
    return percentages

# given a nested list (output of average_all_cluster_compositions), finds the cluster composition averaged over specified frames 
def all_cluster_process_data(nested_lists,group_names_int,mode,frames = 10, threshold = 5):
    last_n = []
    if (mode == "end"):
        last_n = [
            subsublist
            for sublist in nested_lists[-1*frames:]  # Access each sublist in the last 'frames' of nested_lists
            for subsublist in sublist              # Access each subsublist within a sublist
            if len(subsublist) >= threshold  # Include subsublist if it meets the size criterion
         ]
    elif (mode == "middle"):
        temp = len(nested_lists)
        mid = int(temp // 2)
        last_n = [
            subsublist
            for sublist in nested_lists[mid - frames // 2: mid + frames // 2 + 1]  # Access each sublist in the last 'frames' of nested_lists
            for subsublist in sublist              # Access each subsublist within a sublist
            if len(subsublist) >= threshold  # Include subsublist if it meets the size criterion
         ]
    else:
        print("Invalid mode. Should be either middle or end.")
        
    # Flatten the last 5 lists into a single list
    flattened_data = [item for sublist in last_n for item in sublist]
    # Total number of items in the flattened data
    total_items = len(flattened_data)

    if total_items == 0:
        return {i: 0 for i in group_names_int}
    else:
        # Count occurrences of numbers 1 to 6 and calculate percentage
        percentages = {i: (flattened_data.count(i) / total_items) * 100 for i in group_names_int}
    
    return percentages

# finds the percentage of complete clusters in all frames and returns as a list, where a complete cluster refers to a cluster that contains at least one of all types of enzymes
def percentage_of_complete_clusters_in_all_cluster_compositions(combined_list, enz_arr):
    complete_count_list = []
    
    complete_set = set(enz_arr)

    for outerlst in combined_list:
        if len(outerlst) == 0:
            complete_count_list.append(0)
            continue
        complete_count = 0
        for lst in outerlst:
            if set(lst) == complete_set:
                complete_count += 1
                print(lst)
        complete_count_list.append(complete_count / len(outerlst))
            
    return complete_count_list



# finds the average composition of ONLY complete clusters in a specified frame
def cluster_composition_only_complete_clusters(combined_list, enz_arr, mode, frames = 1):
    last_n = []

    if (mode == "end"):
        last_n = combined_list[-1*frames:]
    elif (mode == "middle"):
        temp = len(combined_list)
        mid = int(temp // 2)
        last_n = combined_list[mid - frames // 2: mid + frames // 2 + 1]
    else:
        print("Invalid mode. Should be either middle or end.")
    
    filtered_clusters = []
    complete_set = set(enz_arr)

    for outerlst in last_n:
        temp_complete_list = []
        for lst in outerlst:
            if set(lst) == complete_set:
                temp_complete_list.append(lst)
        filtered_clusters.append(temp_complete_list)
    

    # Flatten the last lists into a single list
    flattened_data = [item for sublist in filtered_clusters for item in sublist]

    percentages = []
    size = []
    for cluster in flattened_data:
        total_items = len(cluster)
        percentage_temp = {i: (cluster.count(i) / total_items) * 100 for i in enz_arr}
        percentages.append(percentage_temp)
        size.append(total_items)

    sum_dict = {key: 0 for key in enz_arr}
    count = len(percentages)

    for d in percentages:
        for key in enz_arr:
            sum_dict[key] += d.get(key, 0)
    try:
        avg_dict = {key: sum_dict[key] / count for key in enz_arr}
    except ZeroDivisionError:
        avg_dict = {key: 0 for key in enz_arr}
    try:
        variance_dict = {key: 0 for key in enz_arr}
    except ZeroDivisionError:
        variance_dict = {key: 0 for key in enz_arr}

    for d in percentages:
        for key in enz_arr:
            variance_dict[key] += (d.get(key, 0) - avg_dict[key]) ** 2
    try:
        variance_dict = {key: variance_dict[key] / count for key in enz_arr}
    except ZeroDivisionError:
        variance_dict = {key: 0 for key in enz_arr}

     # Calculate the standard deviation for each key
    std_dev_dict = {key: math.sqrt(variance_dict[key]) for key in enz_arr}

    results = {key: {"Average": avg_dict[key], "Standard Deviation": std_dev_dict[key]} for key in enz_arr}
    results_df = pd.DataFrame(results).T
   
    return results_df, size

#reads all incomplete clusters and also the enzyme type each is missing
def read_all_incomplete_cluster_compositions(logfile_path, enz_arr, threshold = 5):
    '''Reads all cluster compositions from the clusterList.dat file and returns a list of lists of integers.'''
    
    file_path = os.path.join(logfile_path, 'cluster_results_skip_first_frame/clusterList_corrected.dat')
    data_lines = []
    
    complete_set = set(enz_arr)
    system_data_path = os.path.join(logfile_path, "system.data")
    index_to_name = create_index_to_name_dictionary(system_data_path, "atom")

    containing_enzyme_type_lines = []
    missing_enzyme_type_lines = []
    with open(file_path, 'r') as file:
        # Process each line to extract patterns enclosed by curly braces
        for line in file:
            data_line = []
            missing_enzyme_type_line = []
            containing_enzyme_type_line = []
            # Find all substrings within curly braces and split each match into lists of integers
            matches = re.findall(r'\{(.*?)\}', line)
            
            # Process each match into lists of integers
            processed_lines = [list(map(int, match.split())) for match in matches if len(match.split()) >= threshold]
            
            # Apply the transformation (value // 50) + 1 to each integer in the processed lines
            transformed_lines = [[index_to_name[int(value)] for value in sublist] for sublist in processed_lines]
            
            # Append the transformed lists to data_line
            for index, sublist in enumerate(transformed_lines):
                if set(sublist) != complete_set:
                    zero_index_line = [x - 1 for x in processed_lines[index]]
                    data_line.append(zero_index_line)
                    missing_elements = complete_set - set(sublist)
                    missing_enzyme_type_line.append(missing_elements)
                    containing_enzyme_type_line.append(set(sublist))
            data_lines.append(data_line)
            missing_enzyme_type_lines.append(missing_enzyme_type_line)
            containing_enzyme_type_lines.append(containing_enzyme_type_line)

    return data_lines, missing_enzyme_type_lines, containing_enzyme_type_lines

def read_atom_type_dump(path):
    '''Reads the file from the end until the line "ITEM: ATOMS id type x y z" is encountered.'''
    
    # List all files in the given directory
    files = os.listdir(path)
    
    # Find the file that starts with "atomtype"
    atomtype_file = None
    for file in files:
        if file.startswith('atomtype'):
            atomtype_file = file
            break
    
    if atomtype_file is None:
        raise FileNotFoundError("No file starting with 'atomtype' found in the directory.")
    
    # Construct the full file path
    file_path = os.path.join(path, atomtype_file)
    
    # Initialize an empty dictionary to store the data
    
    # Open the file and read all lines
    with open(file_path, 'r') as file:
        lines = file.readlines()
    
    # Iterate through the lines in reverse order
    for i, line in enumerate(reversed(lines)):
        if line.strip() == "ITEM: ATOMS id type x y z":
            index = len(lines) - 1 - i
            break
    else:
        raise ValueError("The line 'ITEM: ATOMS id type x y z' was not found in the file.")
    
    # Return the remaining lines after the target line
    return lines[index+1:]



def read_atom_type_dump_nth_frame(path, n):
    '''Reads the file until the nth occurrence of the target line is encountered.'''

    # List all files in the given directory
    files = os.listdir(path)
    
    # Find the file that starts with "atomtype"
    atomtype_file = None
    for file in files:
        if file.startswith('atomtype'):
            atomtype_file = file
            break
    
    if atomtype_file is None:
        raise FileNotFoundError("No file starting with 'atomtype' found in the directory.")
    
    # Construct the full file path
    file_path = os.path.join(path, atomtype_file)
    
    count = 0
    
    # Open the file and read line by line
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for i, line in enumerate(lines):
            if line.strip() == "ITEM: ATOMS id type x y z" and count == n - 1:
                count += 1
                index_start = i + 1
            elif line.strip() == "ITEM: ATOMS id type x y z" and count == n:
                index_end = i - 8
                break
            elif line.strip() == "ITEM: ATOMS id type x y z" and count < n:
                count += 1
        else:
            raise ValueError("The line 'ITEM: ATOMS id type x y z' was not found in the file.")
    # Return the remaining lines after the target line
    return lines[index_start:index_end]


def create_type_dictionary(result):
    # Create a defaultdict with an empty list as the default value
    atom_type_dict = defaultdict(list)
    for i in range(1, len(result)):
        parts = result[i].split()
        atom_index = int(parts[0])
        atom_type = int(parts[1])
        atom_type_dict[atom_type].append(atom_index)
    return atom_type_dict

def create_count_dictionary(atom_type_dict):
    atom_type_count_dict = defaultdict(int)
    for key in atom_type_dict:
        atom_type_count_dict[key] = len(atom_type_dict[key])
    return atom_type_count_dict


def create_count_dictionary_faster(file_path, target_frame):
    with open(file_path, 'r') as file:
        for line in file:
            # Strip leading/trailing whitespace and split the line into parts
            parts = line.strip().split()
            if parts and parts[0] == str(target_frame):
                # Extract the desired numbers from the line
                desired_numbers = parts[4:11]
                break
    result_dict = {}
    result_dict[18] = int(desired_numbers[-1])
    result_dict[17] = int(desired_numbers[-2]) - result_dict[18]
    result_dict[16] = int(desired_numbers[-3]) - result_dict[17] - result_dict[18]
    result_dict[15] = int(desired_numbers[-4]) - result_dict[16] - result_dict[17] - result_dict[18]
    result_dict[14] = int(desired_numbers[-5]) - result_dict[15] - result_dict[16] - result_dict[17] - result_dict[18]
    result_dict[13] = int(desired_numbers[-6]) - result_dict[14] - result_dict[15] - result_dict[16] - result_dict[17] - result_dict[18]
    result_dict[12] = int(desired_numbers[-7]) - result_dict[13] - result_dict[14] - result_dict[15] - result_dict[16] - result_dict[17] - result_dict[18]
    result_dict[11] = 1000 - result_dict[12] - result_dict[13] - result_dict[14] - result_dict[15] - result_dict[16] - result_dict[17] - result_dict[18]
    return result_dict


def average_counts_atom_type_faster(paths, target_frame):
    result = []
    for path in paths:
        log_file_full_path = os.path.join(path, 'log.lammps')
        temp = create_count_dictionary_faster(log_file_full_path, target_frame)
        result.append(temp)
    
    sum_dict = {}
    sum_of_squares_dict = {}
    count_dict = {}

    # Sum the values and sum of squares for each key
    for d in result:
        for key in d:
            if key in sum_dict:
                sum_dict[key] += d[key]
                sum_of_squares_dict[key] += d[key] ** 2
                count_dict[key] += 1
            else:
                sum_dict[key] = d[key]
                sum_of_squares_dict[key] = d[key] ** 2
                count_dict[key] = 1

    avg_dict = {key: (sum_dict[key] / count_dict[key] if count_dict[key] != 0 else 0) for key in sum_dict}
    stddev_dict = {key: (np.sqrt((sum_of_squares_dict[key] / count_dict[key]) - (avg_dict[key] ** 2)) if count_dict[key] != 0 else 0) for key in sum_dict}

    return avg_dict, stddev_dict

def average_counts_atom_type_final_frame(paths, keys_to_plot):
    result = []
    for path in paths:
        temp = create_count_dictionary(create_type_dictionary(read_atom_type_dump(path)))
        result.append(temp)

    # Initialize dictionaries to store the sum and sum of squares of values
    sum_dict = {}
    sum_of_squares_dict = {}
    count_dict = {}

    # Sum the values and sum of squares for each key
    for d in result:
        for key in d:
            if key in sum_dict:
                sum_dict[key] += d[key]
                sum_of_squares_dict[key] += d[key] ** 2
                count_dict[key] += 1
            else:
                sum_dict[key] = d[key]
                sum_of_squares_dict[key] = d[key] ** 2
                count_dict[key] = 1

    # Ensure all keys from keys_to_plot are in the dictionaries
    for key in keys_to_plot:
        if key not in sum_dict:
            sum_dict[key] = 0
            sum_of_squares_dict[key] = 0
            count_dict[key] = 0

    # Calculate the average and standard deviation for each key
    avg_dict = {key: (sum_dict[key] / count_dict[key] if count_dict[key] != 0 else 0) for key in sum_dict}
    stddev_dict = {key: (np.sqrt((sum_of_squares_dict[key] / count_dict[key]) - (avg_dict[key] ** 2)) if count_dict[key] != 0 else 0) for key in sum_dict}

    return avg_dict, stddev_dict


def average_counts_atom_type_nth_frame(paths, keys_to_plot, n):
    result = []
    for path in paths:
        temp = create_count_dictionary(create_type_dictionary(read_atom_type_dump_nth_frame(path, n)))
        result.append(temp)

    # Initialize a dictionary to store the sum of values
    sum_dict = {}
    sum_of_squares_dict = {}
    count_dict = {}

    # Sum the values and sum of squares for each key
    for d in result:
        for key in d:
            if key in sum_dict:
                sum_dict[key] += d[key]
                sum_of_squares_dict[key] += d[key] ** 2
                count_dict[key] += 1
            else:
                sum_dict[key] = d[key]
                sum_of_squares_dict[key] = d[key] ** 2
                count_dict[key] = 1

    # Ensure all keys from keys_to_plot are in the dictionaries
    for key in keys_to_plot:
        if key not in sum_dict:
            sum_dict[key] = 0
            sum_of_squares_dict[key] = 0
            count_dict[key] = 0

    # Calculate the average and standard deviation for each key
    avg_dict = {key: (sum_dict[key] / count_dict[key] if count_dict[key] != 0 else 0) for key in sum_dict}
    stddev_dict = {key: (np.sqrt((sum_of_squares_dict[key] / count_dict[key]) - (avg_dict[key] ** 2)) if count_dict[key] != 0 else 0) for key in sum_dict}

    return avg_dict, stddev_dict