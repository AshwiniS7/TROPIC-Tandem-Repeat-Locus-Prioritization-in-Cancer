'''
This program filters repeat loci from the eh-denovo output based on the criteria described in the README. (key difference from the filter_v097 version is that it incorporates donors with files labeled with PCAWG)

Author: Ashwini Suriyaprakash
Modified by:
Version: 5/10/21
Revised: 

'''

import sys
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

## INPUTS
repeat_file = sys.argv[1]
patient_file = sys.argv[2]
manifest_file = sys.argv[3]
filter_constants_file = sys.argv[4]

## OUTPUTS
single_file = sys.argv[5]
paired_file = sys.argv[6]
filtered_file = sys.argv[7]
almostfiltered_file = sys.argv[8]
fi_file = sys.argv[9]

repeats = pd.read_csv(repeat_file,sep='\t',header= None, names=['chr' , 'start', 'stop', 'motif', 'p-value',
                                                                 'corrected_p-value', 'raw_data'] )
## sorts repeats by p-value (helpful for FDR)
repeats = repeats.sort_values(by=['p-value'])
repeats.reset_index(inplace=True, drop = True)
print("eh-denovo RL Count: ",len(repeats))
# print(repeats)

## applies 1st filter based on chromosome
print("applying chromosome filter...")
filtered1 = repeats[(repeats['chr'].str.isnumeric()) | (repeats['chr'] == "X") | (repeats['chr'] == "Y")].copy()
filtered1.reset_index(inplace=True, drop = True)
filtered1['keep'] = 'yes'
print("RL Count: ",len(filtered1))
filtered1.to_csv("debug/f1", sep='\t', index = False, header = True)

patients = pd.read_csv(patient_file,sep='\t',header= None, names=['name' , 'type', 'long_file'] )
patients = patients.drop(['long_file'], axis=1)
patients

manifest = pd.read_csv(manifest_file)
manifest

## collects and formats the normal samples
normal = manifest[manifest['sample_type'].str.contains('Normal', case = False, na = False, regex=False)][['case_id', 'name','sample_type', 'aliquot_id']]
normal = normal.reset_index(drop=True)
normal.replace("\.bam(\.bai)?",".bam",regex=True, inplace=True)
normal['name'] = 'control_' + normal['name']
normal['sample_type'] = 'normal'
normal = normal.drop_duplicates() 
normal

## collects and formats the tumor samples
tumor = manifest[manifest['sample_type'].str.contains('tumour', case = False, na = False, regex= False)][['case_id', 'name','sample_type', 'aliquot_id']]
tumor = tumor.append(manifest[manifest['sample_type'].str.contains('tumor', case = False, na = False, regex= False)][['case_id', 'name','sample_type', 'aliquot_id']])
tumor = tumor.reset_index(drop=True)
tumor.replace("\.bam(\.bai)?",".bam",regex=True, inplace=True)
tumor['name'] = 'case_' + tumor['name']
tumor['sample_type'] = 'tumor'        
tumor = tumor.drop_duplicates() 
tumor

## combines normal and tumor samples for each patient (from the manifest) into one dataframe
all_samples = normal.append(tumor)
all_samples
print(all_samples)

## reading constants
filter_constants = pd.read_csv(filter_constants_file,sep='\t')
MAXIRR = filter_constants.loc[0, 'MAXIRR']
MAXIRR_PERCENT = filter_constants.loc[0, 'MAXIRR_PERCENT']
FDR_CUTOFF = filter_constants.loc[0, 'FDR_CUTOFF']
EXP_IRRDIFF = filter_constants.loc[0, 'EXP_IRRDIFF']
EXP_PERCENT = filter_constants.loc[0, 'EXP_PERCENT'] 

## applies the 2nd filter based on number of anchored IRRs
print("applying max anchored IRR filter...")
for i in range(0,len(filtered1)):
# for i in range(0,2):
    print("Repeat "+ str(i))
    rawdata = filtered1.loc[i,'raw_data']
    chrom = filtered1.loc[i,'chr']
    p_val = filtered1.loc[i,'p-value']
    motif = filtered1.loc[i,'motif']
    start = filtered1.loc[i,'start']
    stop = filtered1.loc[i,'stop']
    
    x = rawdata.split(",")
    
    main_data = pd.DataFrame(columns=['case_id', 'sample_type','anchored_IRRs', 'sample_id'])
    
    ## checks other patients that were input to eh-denovo but did not show up in the output file
    ## assigns 0.0 to their anchored IRRs
    for p in range(len(patients)):
        flag = False
        for file in x:
            if (file.find(patients.loc[p,"name"]) != -1):
                flag = True
        if not (flag):
            x.append(patients.loc[p,"name"] + ":0.0")
    
    
    abnormal_IRR_count = 0
    
    ## maps patient file and anchored IRRs to case_id
    for w in x:
        w = w.lower().strip()
        anchored_IRRs = float(w.split(":")[1])
        fname_orig = w.split(":")[0]
        if "pcawg" in fname_orig:
            fname = fname_orig.replace('pcawg.','PCAWG.')
        else:
            fname = fname_orig
        print("Filename: " + fname)
        df = all_samples[all_samples['name'] == fname].reset_index(drop=True) 
        # print(df)
        if (len(df) > 0):
            print("Found filename in manifest")
            case_id = df.at[0,'case_id']
            sample_type = df.at[0,'sample_type']
            sample_id = df.at[0, 'aliquot_id']
            print(sample_id)
            # print("checking for IRR_max")
            if (anchored_IRRs > MAXIRR):
                abnormal_IRR_count+=1
        
            main_data = main_data.append({'case_id':case_id, 'sample_type':sample_type, 'anchored_IRRs': anchored_IRRs, 'sample_id': sample_id},
                                     ignore_index=True)

    # if there is no patient data, abort and move to the next repeat loci
    if (len(main_data) != 0):
	# print(main_data)
        # print("Sample count with abnormal IRRs:", abnormal_IRR_count)
        # print("Total sample count: ", len(main_data))
        # print("Sample count with expansions: ",exp_count)
        # print("Total paired count: ", len(paired))
        ## generates a file containing individual samples and anchored IRRs
        if (abnormal_IRR_count < MAXIRR_PERCENT/100*len(main_data)):
            main_data.to_csv(single_file +"_" + str(motif)+"_"+str(chrom)+"_" + str(start)+"_"+str(stop) + 
                         ".txt", sep='\t', index = False, header = True)
            filtered1.loc[i, "keep"] = "Yes"
        else:
            print("DROPPING expansion locus")
            filtered1.loc[i, "keep"] = "No"

filtered2 = filtered1[(filtered1['keep'] == "Yes")].copy()
filtered2.reset_index(inplace=True, drop = True)
print("RL Count: ",len(filtered2))
filtered2.to_csv("debug/f12_" + str(MAXIRR), sep='\t', index = False, header = True)
# print(filtered2)

## applies the 3rd filter (FDR q-value)
q_vals = []
print("Applying FDR adjustment...")
for i in range(0,len(filtered2)):
   p_val = (filtered2.loc[i,'p-value'])
   # print("p_val: ", p_val)
   q_val = p_val*len(filtered2)/(i+1)
   # print("Hypothesis:", len(filtered2))
   #print("q_val: ", q_val)
   #print("i+1: ", (i+1))
   q_vals.append(q_val)

filtered2["q-value"] = q_vals

for i in range(len(filtered2)-2, -1, -1):
   q_val = (filtered2.loc[i,'q-value'])
   next_q_val = (filtered2.loc[i+1,'q-value'])
   filtered2.loc[i,'q-value'] = min(q_val, next_q_val)   

filtered2.to_csv("debug/f12_" + str(MAXIRR) + "_3", sep='\t', index = False, header = True)
filtered3 = filtered2[(filtered2['q-value'] < FDR_CUTOFF)].copy()
filtered3.reset_index(inplace=True, drop = True)
# print(filtered3)
print("RL Count: ",len(filtered3))
filtered3.to_csv("debug/f12_" + str(MAXIRR) + "_3_" + str(FDR_CUTOFF), sep='\t', index = False, header = True)


count_list = []
exp_count_list = []
exp_percent_list = []

## applies the 4th filter based on difference of anchored IRRs between tumor and normal samples
print("Applying expansion filter...")    
for i in range(0,len(filtered3)):
# for i in range(0,2):
    # print("Repeat "+ str(i))
    rawdata = filtered3.loc[i,'raw_data']
    chrom = (filtered3.loc[i,'chr'])
    p_val = (filtered3.loc[i,'p-value'])
    motif = (filtered3.loc[i,'motif'])
    start = (filtered3.loc[i,'start'])
    stop = (filtered3.loc[i,'stop'])
    
    orig_file = single_file + "_" + motif + "_" + str(chrom) + "_" + str(start) + "_" + str(stop) + ".txt"
    orig_data = pd.read_csv(orig_file,sep='\t')
    
    sort_samples = orig_data.sort_values(by=['case_id','sample_type'])
    sort_samples.reset_index(inplace = True,drop=True)
    paired = pd.DataFrame(columns=['case_id', 'normal_anchored_IRRs', 'tumor_anchored_IRRs','anchored_IRRs_diff', 'exp', 'normal_sample_id', 'tumor_sample_id'])

    exp_count = 0
    ## Creates tumor and normal pairings
    for j in range(len(sort_samples)-1):
        if ((sort_samples.loc[j,'case_id'] == sort_samples.loc[j+1,'case_id']) and (sort_samples.loc[j, 'sample_type'] == 'normal') and (sort_samples.loc[j+1,'sample_type'] == 'tumor')):
            case_id = sort_samples.loc[j,'case_id']
            normal_IRRs = sort_samples.loc[j,'anchored_IRRs']
            tumor_IRRs = sort_samples.loc[j+1,'anchored_IRRs']
            normal_sample_id = sort_samples.loc[j, 'sample_id']
            tumor_sample_id = sort_samples.loc[j+1, 'sample_id']
            metric = (tumor_IRRs-normal_IRRs)/(normal_IRRs+1)
            if (metric >= EXP_IRRDIFF):
                exp = "Yes"
                exp_count+=1
            else:
                exp = "No"

            paired = paired.append({'case_id':case_id, 'normal_anchored_IRRs':normal_IRRs,
                                'tumor_anchored_IRRs':tumor_IRRs,'anchored_IRRs_diff': metric,
                                   'exp': exp, 'normal_sample_id': normal_sample_id, 'tumor_sample_id': tumor_sample_id},ignore_index=True)

    
    count_list.append(len(paired))
    exp_count_list.append(exp_count)
    exp_percent = exp_count/len(paired) * 100
    exp_percent_list.append(exp_percent)
    # print("Expansion count: ",exp_count)
    # print("Total count: ", len(paired))
    
    ## generates a file containing paired tumor and normal samples and number of anchored IRRs
    if (exp_percent >= EXP_PERCENT):
        filtered3.loc[i, "keep"] = "Yes"
        paired.to_csv(paired_file +"_"+ str(motif)+"_"+str(chrom)+"_" + str(start)+"_"+str(stop) + ".txt", sep='\t', index = False, header = True)
    else:
        filtered3.loc[i, "keep"] = "No"
        os.remove(single_file +"_" + str(motif)+"_"+str(chrom)+"_" + str(start)+"_"+str(stop) +".txt")
        # filtered2.loc[i, "keep"] = "No"

filtered3['paired_count'] = count_list
filtered3['expansion_count'] = exp_count_list
filtered3['expansion_percent'] = exp_percent_list
filtered4 = filtered3[(filtered3['keep'] == "Yes")].copy()
filtered4.reset_index(inplace=True, drop = True)
## storing loci that did not survive the last step
filtered5 = filtered3[(filtered3['keep'] == "No")].copy()
filtered5.reset_index(inplace=True, drop = True)
print("RL Count: ",len(filtered4))
filtered4.to_csv("debug/f12_" + str(MAXIRR) + "_3_" + str(FDR_CUTOFF) + "_4_" + str(EXP_IRRDIFF), sep='\t', index = False, header = True)

## generates a file containing all the repeat loci which passed the filters
filtered4.to_csv(filtered_file, sep='\t', index = False, header = True)
filtered5.to_csv(almostfiltered_file, sep='\t', index = False, header = True)
fi = pd.DataFrame(columns=['eh-denovo_count', 'chrfilter_count','maxIRRfilter_count', 'FDR', 'expfilter_count'])
fi = fi.append({'eh-denovo_count':len(repeats), 'chrfilter_count':len(filtered1), 'maxIRRfilter_count':len(filtered2), 'FDR':len(filtered3), 'expfilter_count':len(filtered4)}, ignore_index=True)
fi.to_csv(fi_file, sep='\t', index = False, header = True)
