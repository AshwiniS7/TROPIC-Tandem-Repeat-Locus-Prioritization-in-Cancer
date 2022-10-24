# TROPIC: Tandem-Repeat-Locus-Prioritization-in-Cancer

Overview
--------

TROPIC prioritizes tandem repeat loci that are associated with different cancers.

The filtering step (implemented by filter_v100.py) eliminates irrelevant loci 
from a list of loci identified by eh-Denovo. It consists of four different filters each of which 
incrementally select loci based on different criterion and pass them to the next filter. The first 
filter selects loci located in chromosomes 1-22, X and Y. The second filter
selects loci which have less samples containing a high number of anchored IRRs. The third filter 
selects loci whose false discovery rate adjusted p-value is less than a certain threshold.
The fourth filter selects loci which have more paired normal and tumor samples having a greater
difference in anchored IRRs. 

{sarea} = source directory containing filter_v100.py

{c} = cancer name (eg: RECA)

{tf_fold} = patientdata

INPUTS:
1. (file) {c}_repeatloci: BED file provided as output from eh-denovo (contains repeat loci that possibly play a role in cancer c) 
2. (file) {c}_patients: tsv file containing the list of filenames and sample IDs for all patients with cancer {c} that were given to eh-denovo
3. (file) manifest: file containing data on all patients from the TCGA and ICGC
4. (file) filter_constants: file containing the thresholds and constants used in filter.py

OUTPUTS:
1. patientdata: folder containing the output data
2. (file) {c}_filtered_repeatloci: BED file containing repeat loci that survived all of the filtering steps in filter.py
3. (file) {c}_filterinfo: text file containing how many loci survived each filtering step executed in filter.py

Execution:
python {sarea}/filter_v100.py input/{c}_repeatloci input/{c}_patients input/manifest input/filter_constants output/{tf_fold}/{c}_single output/{tf_fold}/{c}_paired output/{c}_filtered_repeatloci output/{c}_almostfiltered_repeatloci output/{c}_filterinfo
