#!/bin/python

#Author: Christina Stephens and Stephanie Wankowicz
#Date: October 10, 2018

"""

Average the deviations of all variant fitness scores 

"""

import csv
import argparse
import collections
import itertools
import numpy as np
import sys

def pick_variants(csv_file):
    
    with open(csv_file,'r') as f:
        read_csv = csv.reader(f, delimiter=',')
        dev_t = []
        dev_c = []
        for row1 in read_csv:
            if row1[1] != 'fitness_x':
                dev_t.append(float(row1[2]))
                dev_c.append(float(row1[5]))


    arr_dev_t = np.array(dev_t)
    dev_mean_t = np.mean(arr_dev_t, axis=0)

    arr_dev_c = np.array(dev_c)
    dev_mean_c = np.mean(arr_dev_c, axis=0)

    neg_fitness = {}
    pos_fitness = {}

    with open(csv_file,'r') as f:
        read_csv = csv.reader(f, delimiter=',')
        
        for row in read_csv:
            if row[1] != 'fitness_x':
                if float(row[2]) <= dev_mean_t and float(row[5]) <= dev_mean_c: 
                    if float(row[6]) < 0:
                        neg_fitness[row[3]] = row[6]
                    elif float(row[6]) > 0:
                        pos_fitness[row[3]] = row[6]
    return(neg_fitness,pos_fitness)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""script to generate a pickle file containing counts of library 
    variants""")
    required = parser.add_argument_group('required')
    required.add_argument('-c', '--variant_csv',required=True,
            help="submit csv file containinf fitness, stdv, and variant")
    parser.add_argument('-n', '--name_suffix')
    args = parser.parse_args()

    neg_fitness,pos_fitness = pick_variants(args.variant_csv)

    with open('pos_fitness.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in pos_fitness.items():
            writer.writerow([key, value])

    with open('neg_fitness.csv', 'w') as csv_file:
        writer = csv.writer(csv_file)
        for key, value in neg_fitness.items():
            writer.writerow([key, value])
    
