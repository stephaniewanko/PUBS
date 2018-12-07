#!/bin/python

#Author: Christina Stephens and Stephanie Wankowicz
#Date: October 10, 2018

"""

(1) Calculate all wildtype (WT) barcodes freq. normalized to total sum of barcode counts for each timepoint
(2) Remove outlier barcode counts
(3) Log transform normalized frequencies
(4) Fit points (linear regression plugin)
(5) Select well-fitted slopes
(6) Select the medium slope and associated barcode


"""

import csv
import argparse
import collections
import itertools
import numpy as np
from scipy.stats import linregress
import statistics
import pickle

def variant_counter_from_fastqs(csv_files,barcode_variant_dict):
    """generate a dictionary of variants and their wt-normalized frequencies for 
        three time points from an input csv file organized by barcode    
    """

    #variant_counter = collections.defaultdict(lambda: [[],[]])
    #variant_counter_raw = collections.defaultdict(lambda: [[],[]])

    variant_counter = collections.defaultdict(lambda: [[],[],[]])
    variant_counter_raw = collections.defaultdict(lambda: [[],[],[]])
    for i, csv_fn in enumerate(csv_files):
        sum_t = 0
        with open(csv_fn,'r') as f:
            read_csv = csv.reader(f, delimiter=',')
            for row in read_csv:
                if row[1] != 'count':
                    sum_t += float(row[1])
   
        
        dict_keys = list(barcode_variant_dict.keys()) 
        with open(csv_fn,'r') as f:
            read_csv = csv.reader(f, delimiter=',')
            for row1 in read_csv:
                if row1[1] != 'count':
                    barcode = row1[0][:20]
                    count=row1[1]
                    if barcode in dict_keys:
                        if barcode_variant_dict[barcode][1] == 'WT':
                            variant_counter[barcode][i].append(float(count)/sum_t)
                            variant_counter_raw[barcode][i].append(float(count))
    return variant_counter, variant_counter_raw

def calculate_variant_fitness(variant_timepoint_counter,variant_counter_raw):
    """uses method described in Doug Folwer's Enrich2 paper - fitness is the slope of linear regression line"""
    s = {}
    count = 0
    count1 = 0
    r_squared = []
    only_good_r_squared = []
    all_counts = [[],[],[]]
    for barcode, count_list in variant_timepoint_counter.items():
        
        all_counts[0] += variant_counter_raw[barcode][0]
        all_counts[1] += variant_counter_raw[barcode][1]
        all_counts[2] += variant_counter_raw[barcode][2]
        
    #median_count = []
    #for i,time in enumerate(all_counts):
    #    l = len(list(time))
    #    if l % 2 == 0:
    #        list1 = sorted(time, key=int)
    #        middle = float(l/2)
    #        median_count.append(list1[int(middle - .5)])
    #    else:
    #        median_count.append(np.median(time)) 
    
        t0 = np.array(all_counts[0])
        t0_mean = np.mean(t0, axis=0)
        t0_std = np.std(t0, axis=0)
    
        t1 = np.array(all_counts[1])
        t1_mean = np.mean(t1, axis=0)
        t1_std = np.std(t1, axis=0)
    
        t2 = np.array(all_counts[2])
        t2_mean = np.mean(t2, axis=0)
        t2_std = np.std(t2, axis=0)

    for barcode, count_list in variant_timepoint_counter.items():
        
        count+=1
        
        # if the array is empty add a single value of 0.5
        if len(count_list[0]) == 0:
            count_list[0].append(1)
            variant_counter_raw[barcode][0].append(0)
        if len(count_list[1]) == 0:
            count_list[1].append(1)
            variant_counter_raw[barcode][1].append(0)
        if len(count_list[2]) == 0:
            count_list[2].append(1)
            variant_counter_raw[barcode][2].append(0)

        get_out=0    
        #if median_count[0] <= variant_counter_raw[barcode][0][0]:   
        #    get_out = 1
        #if median_count[1] <= variant_counter_raw[barcode][1][0]:   
        #    get_out = 1
        #if median_count[2] <= variant_counter_raw[barcode][2][0]:   
        #    get_out = 1
  

        if (t0_mean - t0_std) > variant_counter_raw[barcode][0][0] or variant_counter_raw[barcode][0][0] > (t0_mean + t0_std):   
            get_out = 1
        if (t1_mean - t1_std) > variant_counter_raw[barcode][1][0] or variant_counter_raw[barcode][1][0] > (t1_mean + t1_std):   
            get_out = 1
        if (t2_mean - t2_std) > variant_counter_raw[barcode][2][0] or variant_counter_raw[barcode][2][0] > (t2_mean + t2_std):   
            get_out = 1

        if get_out == 0:
            l0 = [0] * len(count_list[0])
            l1 = [1] * len(count_list[1])
            l2 = [2] * len(count_list[2])

            f = count_list[0] + count_list[1] + count_list[2]
            x = l0 + l1 + l2

            freq = np.array(f)
            m_vt = np.log(freq)
            for i,item in enumerate(m_vt):
                if item == np.log(1):
                    m_vt[i] = 0
            slope, intercept, r_value, p_value, std_err = linregress(x, m_vt)
            r_squared.append(r_value*r_value)

            if r_value*r_value > 0.9:
                s[barcode] = slope
                only_good_r_squared.append(r_value*r_value)
        
    l = len(list(s.values()))
    if l % 2 == 0:
        list1 = sorted(list(s.values()), key=float)
        middle = float(l/2)
        medium_score = list1[int(middle - .5)]
    else:
        medium_score = np.median(list(s.values()))

    for i,v in s.items():
        if s[i] == medium_score:
            medium_barcode = i

    print("Total number of WT barcodes: " + str(count))
    print("Number of selected WT barcodes by count: " + str(len(r_squared)))
    print("Number of selected WT barcodes by r_squared and count: " + str(len(only_good_r_squared)))

    return (medium_barcode, variant_counter_raw[medium_barcode],r_squared)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""script to generate a pickle file containing counts of library 
    variants""")
    required = parser.add_argument_group('required')
    required.add_argument('-c', '--barcode_csv',  nargs='*',required=True,
                          help='3 preprocessed csv files containing sequenced barcodes and raw counts for single timepoints '
                               'it is assumed that the order of fastq files corresponds to the order of time points')
    required.add_argument('-b', '--barcode_pickle', required=True,
                          help='pickle containing a dictionary with barcode as keys and library variants as values. '
                               'reads in fastq file will be searched for exact matches to dictionary keys')
    
    parser.add_argument('-n', '--name_suffix')
    args = parser.parse_args()
    with open(args.barcode_pickle, 'rb') as f:
        barcode_variant_dict = pickle.load(f)

    variant_timepoint_counter,variant_counter_raw = variant_counter_from_fastqs(args.barcode_csv, barcode_variant_dict)
    
    medium_barcode,raw_counts,r_squared = calculate_variant_fitness(variant_timepoint_counter,variant_counter_raw)
    
    with open('r_squared_file.txt', 'w') as f:
        for item in r_squared:
            f.write("%s\n" % item)

    print("")
    print("Reference barcode: " + str(medium_barcode))
    print("raw count for 0hr: " + str(raw_counts[0]) + "  12hr: " + str(raw_counts[1]) + "  24hr: " + str(raw_counts[2]))
    print("")
