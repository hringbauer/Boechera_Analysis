'''
Created on May 28th, 2018
Transform .ibd File to Genetic Map

Input: Name of .ibd file in IBD_raw/
Requires linkage map in "../Data/LinkageMap/"	 

Output: Save transformed .ibd file to IBD_cm/

@ Author Harald Ringbauer
'''

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os
import sys

############################
############################ Define the File Names
arguments = sys.argv
if len(arguments)!= 3:
	print("Nr Arguments given: %i" % len(arguments))
	raise ValueError("Please give TWO arguments. LM path and .ibd path")

ld_file = str(arguments[1]) # Which Linkage Map to use
ibd_file = str(arguments[2]) # Which IBD file to analyze

## Output File:  
cm_file = "./Data/IBD_cm/" + os.path.basename(ibd_file)

############################
############################ Load the Data

### Load the genetic map
l_df = pd.read_csv(ld_file, sep="\t") 
ld_columns = ["LG", "Scaffold", "Bin", "cm", "Bin_End", "Method", "Start", "End"]
l_df.columns = ld_columns

### Load the IBD file
df_ibd = pd.read_csv(ibd_file, sep="\t", header=None)
columns_ibd = ["Ind1", "HapIndex1", "Ind2", "HapIndex2","Scaffold", "IBDStart", "IBDEnd", "Lod", "IBDlen"] 
df_ibd.columns = columns_ibd


############################
############################ Run the actual Code
# For every Linkage Group make mapping df. I.e. make list of dataframe.

mapping=[] # List of mapping functions!! (i.e. scipy interpolation object)
max_bps = [] # The upper limits (mainly for plotting)

## For the moment only try out first linkage group:

for lg in range(1,8):

    l_df1 = l_df.loc[l_df["LG"]==lg]
    assert(len(l_df1)>0)
    dup = l_df1["cm"].duplicated(keep="first")
    data = l_df1[~dup][["cm", "Bin_End"]].values
    bp = data[:, 1]
    cm = data[:, 0]

    ### Create the Interpolation object:
    f = interp1d(bp, cm, fill_value="extrapolate")
    mapping.append(f)  # add to mapping
    max_bps.append(np.max(bp))
    #f2 = interp1d(bp, cm, kind='cubic', fill_value="extrapolate")



### Do a plot
fig, axarr = plt.subplots(4, 2, figsize=(10,12))
axarr = axarr.flatten() # Flatten the array for better Looping


for i in range(len(mapping)):
    f = mapping[i]
    ax = axarr[i]

    xnew = np.linspace(0, max_bps[i], num=100, endpoint=True)

    # Extract the Data from linkage map:
    l_df1 = l_df.loc[l_df["LG"]==(i+1)]
    assert(len(l_df1)>0)
    dup = l_df1["cm"].duplicated(keep="first")
    data = l_df1[~dup][["cm", "Bin_End"]].values
    bp = data[:, 1]
    cm = data[:, 0]

    ax.plot(bp, cm, 'o', xnew, f(xnew), '-', linewidth=2)
    ax.set_title("LG: %i" % (i+1), fontsize=12, pad=2)
axarr[0].legend(['Linkage Map', 'Interpolation'], loc='upper left', fontsize=16)

# Set the Labels
fig.text(0.5, 0.04, 'bp', ha='center', va='center', fontsize=12)
fig.text(0.06, 0.5, 'cM', ha='center', va='center', rotation='vertical', fontsize=14)

#plt.show()

fig.savefig('l_map.png', bbox_inches='tight', pad_inches=0)
print("Mapping Figure saved!")



####### Print some Summary Statistics
tot_mb = np.sum(max_bps) / (1e6)
print("Total lengt (mb): %.3f" % tot_mb)  # Only 180 mega bases???

tot_cm = np.sum([mapping[i](max_bps[i]) for i in range(7)])
print("Total length (cm): %.3f" % tot_cm)

print("mb per cm: %.4f" % (tot_mb/tot_cm))

####### Do the actual transformation

df_cm = df_ibd.copy() # Create a new copy

for i, df_scaffold in df_ibd.groupby('Scaffold'):
    print("\nDoing saffold %i ..." % i)
    print("Nr of blocks: %i" % len(df_scaffold))
    
    # Calculate the new Coordinates.
    starts = mapping[i-1](df_scaffold["IBDStart"].values)
    ends = mapping[i-1](df_scaffold["IBDEnd"].values)
    
    print("Nr of Blocks negative length: %i " % np.sum((ends-starts)<0))
    corr_ends = np.maximum(starts+0.05, ends) # Set the end values to be at least 0.05 cM after start values
    
    inds = df_scaffold.index
    df_cm.loc[inds, "IBDStart"] = starts
    df_cm.loc[inds, "IBDEnd"] = corr_ends
    
df_cm["IBDlen"] = df_cm["IBDEnd"]-df_cm["IBDStart"]
print("Finished Mapping.")

df_cm.to_csv(cm_file, index=False, float_format='%.3f')
print("Successfully saved to %s" % cm_file)
assert(len(df_cm)>0)










