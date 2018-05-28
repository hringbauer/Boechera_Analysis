'''
Created on May 28th, 2018
Removes Gaps from cm Folder! .ibd File has to be in cM units

Input: Name of .ibd file in IBD_cm/.

Output: Save cleaned .ibd file in IBD_merged/ with name+.merged.ib

@ Author Harald Ringbauer
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import _pickle as pickle
import itertools
import bisect
import os
import sys

max_gap = 10.0  # Maximum Length of Gaps to merge (Up to 1 cM)

############################
############################ Define File Names
arguments = sys.argv
if len(arguments)!= 2:
	print("Nr Arguments given: %i" % len(arguments))
	raise ValueError("Please give ONE argument")

# IBD Input File
ibd_file_path = str(arguments[1]) 

# Output File to save to
merge_file = "./Data/IBD_merged/" + os.path.basename(ibd_file_path)

############################
############################ Load the Data

### Load the IBD data
df = pd.read_csv(ibd_file_path)
# Set to new column Values
columns = ["Ind1", "HapIndex1", "Ind2", "HapIndex2","Scaffold", "IBDStart", "IBDEnd", "Lod", "IBDlen"] 
df.columns = columns

# Calculate IBD Len column
df["IBDlen"] = df["IBDEnd"] - df["IBDStart"]  # Lenght column redundant but faster analysis
assert(np.min(df["IBDlen"])>0) # Sanity Check

print("\nSuccessfully loaded IBD Dataset from:\n %s" % ibd_file_path)
print("\nNr of Blocks: %i" % len(df))
print("Nr of Blocks > 5 cm: %i " % len(df[df.IBDlen>5]))
print("Nr of Blocks > 10 cm: %i " % len(df[df.IBDlen>10]))
print("Nr of Blocks > 20 cm: %i " % len(df[df.IBDlen>20]))

############################
############################ Do the main Job
print("\nStarting the script.")
print("Maximum Gap that is closed: %i cm" % max_gap)

# Order raws (after Ind1, then Ind2 then Chromosome)!!
df.sort_values(['Ind1', 'Ind2', "Scaffold", "IBDStart"], ascending=[True, True, True, True], inplace=True) # Sort Dataframe

df_d = df # Do start with the original dataframe

print("Starting filling Gaps...")

for i in range(10):
    print("Doing Cycle: %i" %i)
    print("Length of Dataframe: %i" % len(df_d))
    same = df_d.duplicated(subset=['Ind1','Ind2',"Scaffold"], keep='first') # True if same as previous (but first)
    same = np.where(same==True)[0]
    print("Nr of blocks on same Pair of Individuals and Scaffold: %i" % len(same))

    # Flag raws where gap could be filled (index is shifted to second raw - since always comparison to previous one)
    start = df_d["IBDStart"][same]
    end = df_d["IBDEnd"][same-1]
    block_bb = np.maximum(df_d["IBDlen"][same].values, df_d["IBDlen"][same-1].values) # Extract the longest block
    gaps = start.values-end.values

    # Whether gap smaller than block or previous block
    # Gap is smaller than Maximum Gap. And also smaller than max. Block. Catches also overlap (gap<0)
    fill = (gaps < max_gap) * (gaps < block_bb)  

    print("Blocks with small enough gaps to merge: %i" % np.sum(fill))
    if np.sum(fill)==0:
        print("Stoppinig after %i repetitions." %i)
        break

    # Extract indices of all starting points:
    block_uq = np.ones(len(df_d), dtype=bool) # The Boolean array of unique Blocks
    block_uq[same[fill]]=0  # Set the filler blocks to 0
    unique = np.where(block_uq==True)[0]

    df_m  = df_d.iloc[unique,:].reset_index(drop=True)  # To avoid carrying on previous index

    # Fill end points with end point of block before next starting point
    # Create Starting and End Values for overlapping blocks:
    starts = np.split(df_d["IBDStart"], unique)[1:] # Extract the Start of Filler blocks
    min_starts = np.array([np.min(i) for i in starts])


    ends =  np.split(df_d["IBDEnd"], unique)[1:]  # Extract the End points of Filler blocks
    max_ends = np.array([np.max(i) for i in ends])

    df_m["IBDStart"] = min_starts
    df_m["IBDEnd"] = max_ends
    #df_m.loc[:len(df_m)-2, "IBDEnd"] = df_d["IBDEnd"].iloc[unique-1][1:].values


    df_m["IBDlen"]=df_m["IBDEnd"]-df_m["IBDStart"]
    assert(np.min(df_m["IBDlen"])>0) # Sanity check
    print("\nLength of merged Dataframe: %i" % len(df_m))
    df_d = df_m # Prepare next Cycle   


# Do the saving of IBD blocks:
print("\nNr of Blocks: %i" % len(df_m))
print("Nr of Blocks > 5 cm: %i " % len(df_m[df_m.IBDlen>5]))
print("Nr of Blocks > 10 cm: %i " % len(df_m[df_m.IBDlen>10]))
print("Nr of Blocks > 20 cm: %i " % len(df_m[df_m.IBDlen>20]))



df_m.to_csv(merge_file, float_format='%.4f')
print("Successfully saved to %s" %merge_file)

assert(np.min(df_m["IBDlen"])>0) # Sanity Check


