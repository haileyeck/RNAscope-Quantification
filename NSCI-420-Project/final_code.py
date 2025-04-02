# This is a python program for quantifying the intensity of mRNA in cells. 
# It assumes Ch2 contains a viral marker to indicate the virus was delivered to a cell (e.g., Sacas9),
# and Ch3 contains a probe marker which is the measured mRNA intensity (e.g., Syngap1). 
# It splits the cells into Sacas9+ and Sacas9-. 

import pandas as pd
import os
import random

# fill in the variables below with the path to each folder
raw_measurements_folder = '/Users/haileyeckersberg/Desktop/syngap project'
ctcf_folder = '/Users/haileyeckersberg/syngap project test'

###
###
# CHANGE NAMES OF THESE
###
###
output_folder = '/Users/haileyeckersberg/syngap project test output'
output_folder_2 = '/Users/haileyeckersberg/syngap project test/greatest_ctcfs'

os.makedirs(output_folder, exist_ok=True)
os.makedirs(output_folder_2, exist_ok=True)


# checks if your input folder exists
os.makedirs(raw_measurements_folder, exist_ok=True)

# for this code to work, your files should be named based on the following scheme:
# 1. each animal has a representative letter 
# 2. there is a number representing the slice you measured
# 3. there is a number representing each image 
# e.g. the first image you take of the first slice corresponding to animal A will have the identifier:
#      A1_1
# the background measurements for each image follow a similar scheme where you'd have: 
#      BG_A1_1

# based on set naming system:
# 1. add a letter representing each animal to the animals string  
# 2. set a range for the number of slices each animal has
animals = 'ABCDEFGH'
slices = range(1, 4)

# indicate which animals are controls and which are treated
control = 'BDFH'
treated = 'ACEG'

# list to hold measurements for each animal
all_cells_data = []
#all_cells_data_wo_neg = []


# FIRST, 
# check for file formatting errors, calculate CTCFs, and save CSVs with CTCFs
for animal in animals:
    for slice_num in slices:
        image_num = 1

        while True: 
            image_id = f"{animal}{slice_num}_{image_num}"

            bg_filepath = f"{raw_measurements_folder}/BG_{image_id}.csv"
            main_filepath = f"{raw_measurements_folder}/{image_id}.csv"

            try:
                bg_data = pd.read_csv(bg_filepath)

                # calculates mean bg from bg csv 
                mean_bg_ch2 = (bg_data[bg_data['Ch'] == 2]['IntDen'] / bg_data[bg_data['Ch'] == 2]['Area']).mean()
                mean_bg_ch3 = (bg_data[bg_data['Ch'] == 3]['IntDen'] / bg_data[bg_data['Ch'] == 3]['Area']).mean()

                main_data = pd.read_csv(main_filepath)

                # Initialize a list to store CTCF calculations
                ctcf_ch2_data = []
                ctcf_ch2_data_without_neg = []
                ctcf_ch3_data = []
                ctcf_ch3_data_without_neg = []

                ctcf_both_data = []
                ctcf_both_data_without_neg = []

                # iterate over rows in pairs (Ch2 and Ch3 for each region)
                for i in range(0, len(main_data), 2):
                    # next bit will check for errors in the csv files
                    # if there aren't measurements for both channels for a partiular region
                    if i + 1 >= len(main_data):
                        print(f"Incomplete pair for region in {image_id}. Region number: {(i // 2) + 1}")
                        exit(1)
                    
                    # if two consecutive rows do not have Ch = 2 and Ch = 3
                    ch2_row = main_data.iloc[i]
                    ch3_row = main_data.iloc[i + 1]
                    if ch2_row['Ch'] != 2 or ch3_row['Ch'] != 3:
                        print(f"Pattern error in {image_id}: Expected Ch2-Ch3 pair, got Ch{ch2_row['Ch']}-Ch{ch3_row['Ch']}. Region number: {(i // 2) + 1}")
                        exit(1)
                    
                    # if the region areas don't match between two consecutive rows 
                    if ch2_row['Area'] != ch3_row['Area']:
                        print(f"Area mismatch in {image_id} for region {(i // 2) + 1}")
                        exit(1)
                    
                    region_num = (i // 2) + 1
                    
                    # calculate CTCF for Ch2 and Ch3
                    ctcf_ch2 = ch2_row['IntDen'] - (ch2_row['Area'] * mean_bg_ch2)
                    ctcf_ch3 = ch3_row['IntDen'] - (ch3_row['Area'] * mean_bg_ch3)
                    
                    # append ctcf results to master lists
                    ctcf_ch2_data.append([region_num, ctcf_ch2])
                    if ctcf_ch2 > 0: 
                        ctcf_ch2_data_without_neg.append([region_num, ctcf_ch2])
                    
                    ctcf_ch3_data.append([region_num, ctcf_ch3])
                    if ctcf_ch3 > 0: 
                        ctcf_ch3_data_without_neg.append([region_num, ctcf_ch2])

                    # add to this whether neg or not
                    ctcf_both_data.append([region_num, ctcf_ch2, ctcf_ch3])

                    # check if CTCF values are non-negative for both channels
                    if ctcf_ch2 >= 0 and ctcf_ch3 >= 0:
                        # append result to the final list only if both values are non-negative
                        ctcf_both_data_without_neg.append([region_num, ctcf_ch2, ctcf_ch3])
                    
                    all_cells_data.append({
                        'Identifier': image_id,'Region': (i // 2) + 1,'CTCF_Ch2': ctcf_ch2,'CTCF_Ch3': ctcf_ch3          
                    })

                    #if ctcf_ch2 > 0 and ctcf_ch3 > 0:
                    #    all_cells_data_wo_neg.append({
                    #        'Identifier': image_id,'Region': (i // 2) + 1,'CTCF_Ch2': ctcf_ch2,'CTCF_Ch3': ctcf_ch3          
                    #    })

                # convert results to dataframes
                ctcf_df2 = pd.DataFrame(ctcf_ch2_data, columns=['Region', 'CTCF_Ch2'])
                ctcf_df2_no_neg = pd.DataFrame(ctcf_ch2_data_without_neg, columns=['Region', 'CTCF_Ch2'])
                ctcf_df3 = pd.DataFrame(ctcf_ch3_data, columns=['Region', 'CTCF_Ch3'])
                ctcf_df3_no_neg = pd.DataFrame(ctcf_ch3_data_without_neg, columns=['Region', 'CTCF_Ch2'])
                
                ctcf_df23 = pd.DataFrame(ctcf_both_data, columns=['Region','CTCF_Ch2','CTCF_Ch3'])
                ctcf_df23_no_neg = pd.DataFrame(ctcf_both_data_without_neg, columns=['Region','CTCF_Ch2','CTCF_Ch3'])
                

                # save CTCF values to new csv
                # CURRENTLY SET TO FILES WITH POS AND NEG VALUES
                output_file2 = f"{ctcf_folder}/CTCF_{image_id}_Sacas9.csv"
                output_file3 = f"{ctcf_folder}/CTCF_{image_id}_SynGAP.csv"
                output_file23 = f"{ctcf_folder}/CTCF_{image_id}_both.csv"
                ctcf_df2.to_csv(output_file2, index=False)
                ctcf_df3.to_csv(output_file3, index=False)
                ctcf_df23.to_csv(output_file23, index=False)
                print(f"Processed and saved CTCF data for {image_id}.")

                # IF WANT ONLY POSITIVE, USE THIS CODE
                #output_file2_no_neg = f"{ctcf_folder}/CTCF_{identifier}_Ch2_no_neg.csv"
                #output_file3_no_neg = f"{ctcf_folder}/CTCF_{identifier}_Ch3_no_neg.csv"
                #output_file23_no_neg = f"{ctcf_folder}/CTCF_{identifier}_Ch23_no_neg.csv"
                #ctcf_df2_no_neg.to_csv(output_file2_no_neg, index=False)
                #ctcf_df3_no_neg.to_csv(output_file3_no_neg, index=False)
                #ctcf_df23_no_neg.to_csv(output_file23_no_neg, index=False)
                #print(f"Processed and saved non-negative CTCF data for {image_id}.")


            # will check if file exists and stop iteration if no more images for this slice
            except FileNotFoundError:
                break
            
            # go to next image
            image_num += 1


#IF ONLY WANT POS CTCFs, CHANGE THIS
all_cells_df = pd.DataFrame(all_cells_data)
#all_cells_df = pd.DataFrame(all_cells_data_wo_neg)

ctcf_csv_path = os.path.join(ctcf_folder, 'All_Cells_CTCF.csv')
all_cells_df.to_csv(ctcf_csv_path, index=False)

print(f"All cells' CTCF data saved to {ctcf_csv_path}.")

            
greatest_ctcf_ch2 = []
greatest_ctcf_ch2_for_avg = 0
total_images = 0

# list that will hold Ch3 CTCF for all controls
# will be used later in random sampling
ctcf_ch3_controls = []

# NEXT:
# finds greatest CTCF for Ch2 (Sacas9) in control animals
for animal in control:
    for slice_num in slices:
        image_num = 1
       
        print(f"On animal: {animal}.")
        
        while True:
            print(f"On image: {image_num}.")

            identifier = f"{animal}{slice_num}_{image_num}"

            both_file_path = f"{ctcf_folder}/CTCF_{identifier}_both.csv"

            try:
                df_both = pd.read_csv(both_file_path)
                ctcf_ch3_controls.extend(df_both['CTCF_Ch3'].tolist())
               
                # sort and append the greatest CH2 CTCF value
                sorted_df_both = df_both.sort_values(by=["CTCF_Ch2"], ascending=False)
                greatest_ctcf_ch2.append([identifier, sorted_df_both['CTCF_Ch2'].iloc[0]])
                greatest_ctcf_ch2_for_avg += sorted_df_both['CTCF_Ch2'].iloc[0]
                total_images += 1

                # do i need to save this? 
                # save sorted DataFrame
                output_file2 = f"{output_folder}/CTCF_{identifier}_sorted_by_sacas9.csv"
                sorted_df_both.to_csv(output_file2, index=False)
                
                print(f"Processed and saved Ch2 CTCF sort data for {identifier}.")

            except FileNotFoundError:
                break
            
            image_num += 1

#####
#####
#####
# do i need this??? i don't think so...
# save greatest CTCF values for channel 2
greatest_ch2_ctcfs = pd.DataFrame(greatest_ctcf_ch2, columns=['Identifier', 'CTCF_Ch2'])
greatest_ch2_ctcfs_file = f"{output_folder_2}/greatest_CTCF_Ch2.csv"
greatest_ch2_ctcfs.to_csv(greatest_ch2_ctcfs_file, index=False)

# compute average greatest Ch2 CTCF and set as Sacas9+ threshold
avg_greatest_ch2_ctcf = greatest_ctcf_ch2_for_avg / total_images
print(f"Avg of all greatest Sacas9 CTCFs/threshold for positive Sacas9: {avg_greatest_ch2_ctcf}.")

sacas9_threshold = avg_greatest_ch2_ctcf

# NEXT:
# sort treated animals to Sacas9+ and Sacas9-

ctcf_ch3_sacas9_positive = []
ctcf_ch3_sacas9_negative = []

for animal in treated:
    for slice_num in slices:
        image_num = 1

        print(f"On animal: {animal}.")

        while True:
            identifier = f"{animal}{slice_num}_{image_num}"
            both_file_path = f"{ctcf_folder}/CTCF_{identifier}_both.csv"

            try:
                df_both = pd.read_csv(both_file_path)

                # create new row to indicate if cell is Sacas9+(=1) or Sacas9-(=0)
                df_both['Sacas9 (+ or -)'] = df_both['CTCF_Ch2'].apply(lambda x: 1 if x >= sacas9_threshold else 0)

                ctcf_ch3_sacas9_positive.extend(df_both[df_both['Sacas9 (+ or -)'] == 1]['CTCF_Ch3'].tolist())
                ctcf_ch3_sacas9_negative.extend(df_both[df_both['Sacas9 (+ or -)'] == 0]['CTCF_Ch3'].tolist())


                output_file_labeled = f"{output_folder}/CTCF_{identifier}_Sacas9_thresholded.csv"
                df_both.to_csv(output_file_labeled, index=False)

                print(f"Processed and labeled cells for {identifier} based on Sacas9 threshold.")

            except FileNotFoundError: 
                break

            image_num += 1

# NEXT: 
# need to do random sampling and then fix final_data below...

sample_size = min(len(ctcf_ch3_sacas9_positive), len(ctcf_ch3_sacas9_negative), len(ctcf_ch3_controls))

sampled_ctcf_ch3_sacas9_positive = random.sample(ctcf_ch3_sacas9_positive, sample_size)
sampled_ctcf_ch3_sacas9_negative = random.sample(ctcf_ch3_sacas9_negative, sample_size)
sampled_ctcf_ch3_controls = random.sample(ctcf_ch3_controls, sample_size)



final_data = {
    'SaCas9+': sampled_ctcf_ch3_sacas9_positive,
    'SaCas9-': sampled_ctcf_ch3_sacas9_negative,
    'Untreated': sampled_ctcf_ch3_controls
}
final_df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in final_data.items()]))

# save to CSV
output_file = f"{output_folder}/categorized_CTCF_Ch3.csv"
final_df.to_csv(output_file, index=False)

print("Finished organizing and saving categorized CTCF_Ch3 data.")
