#Author: Elissa Fink
#Date: October 19, 2018

"""

Step 4 of Fitness Calculation Protocol for PUBS 2018
(1) Calculate fitness for each barcode
(2) Calculate average fitness for each variant

python3 step4_collapse_variant_fitness.py <input.csv>

Intermediate output: 'intermediate_csv_bc_fitness_' + input_name.csv
Output is 'collapsed_fitness_' + input_name.csv

"""

import numpy as np
import pandas as pd
from scipy import stats
import math
import sys

"""Calculates average slope and standard deviation for all barcodes of unique variant"""
def weighted_calc(slopes):
    slope_mean = np.mean(slopes)
    slope_sd = np.std(slopes)
    return (slope_mean,slope_sd)

"""Generates fitness, R^2, standard error for each unique barcode"""
def barcode_fitness(df_in):
    time = []
    ratio = []
    var_num = []
    var_aa = []
    x_ind = []
    bc_info = []
    for item in df_in.iterrows():
        row = list(item[1])
        time.append(int(row[5]))
        ratio.append(row[6])
        var_num.append(int(row[3]))
        var_aa.append(row[4])
    for t in time:
        if t == 0:
            x_ind.append(0)
        if t == 12:
            x_ind.append(1)
        if t == 24:
            x_ind.append(2)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x_ind,ratio)
    bc_info.append(slope)        
    bc_info.append(r_value**2)
    bc_info.append(std_err)
    bc_info.append(var_aa[0])
    bc_info.append(var_num[0])
    return (bc_info)

"""Main"""
#open file from command line
f = open(sys.argv[1],"r")
df = pd.read_csv(f,sep=',')
#df1 = df.drop(['Unnamed: 0'],axis=1)
df2 = df.sort_values(by=['barcode'])
unique_bc = (df2.barcode.unique())

#build new df with bc slopes and stats
df_slope = []
df_rsq = []
df_sterr = []
df_var_aa = []
df_var_num = []
df_variant = []
print ("Starting barcode fitness calculations...")
for i in range(0,len(unique_bc)):
    temp_df = df2.loc[df2['barcode'] == unique_bc[i]]
    bc_info = barcode_fitness(temp_df)
    df_slope.append(bc_info[0])
    df_rsq.append(bc_info[1])
    df_sterr.append(bc_info[2])
    df_var_aa.append(bc_info[3])
    df_var_num.append(bc_info[4])
    combo = "".join([str(bc_info[4]),bc_info[3]])
    df_variant.append(combo)
df3 = pd.DataFrame({'barcode':unique_bc,'full_variant':df_variant,'variant_aa':df_var_aa,'variant_num':df_var_num,'fitness_slope':df_slope,'r_squared':df_rsq,'st_err':df_sterr})
df4 = df3.sort_values(['variant_num','variant_aa'],ascending=[True,True])
name_in = str(sys.argv[1])
intermediate_name = 'intermediate_csv_bc_fitness_'
int_save = "".join([intermediate_name,name_in])
df4.to_csv(int_save,sep=',')
print ("Barcode fitness calculations completed and saved in intermediate csv.")

#to collapse fitness to variant
print ("Starting average fitness calculation for variants...this will take a while!")
unique_variants = list(df4.full_variant.unique())
final_fitness = {}
for item in unique_variants:
    all_slopes = []
    for row in df4.iterrows():
        row_list = list(row[1])
        row_var = row_list[2]
        if item == row_var:
            all_slopes.append(row_list[1])
    collapsed_info = weighted_calc(all_slopes)
    final_fitness[item]=collapsed_info
print ("Finished average fitness calculation for variants...")

#make df of variant fitness
df_dict = pd.DataFrame.from_dict(final_fitness, orient = 'index')
df_dict_1 = df_dict.reset_index()
df_dict_1.columns = ['variant', 'fitness','stdev']

#add WT aa to variant
var_pre_wt = list(df_dict_1['variant'])
full_variant_wt = []
wt = '*MDVFMKGLSKAKEGVVAAAEKTKQGVAEAAGKTKEGVLYVGSKTKEGVVHGVATVAEKTKEQVTNVGGAVVTGVTAVAQKTVEGAGSIAAATGFVKKDQLGKNEEGAPQEGILEDMPVDPDNEAYEMPSEEGYQDYEPEA'
for var in var_pre_wt:
    if int(var[0]) == 0:
        full_variant_wt.append(var)
    else:
        temp = int(var[:-1])
        wt_aa = wt[temp]
        combo = "".join([wt_aa,var])
        full_variant_wt.append(combo)

#new df with full variant
df_new = df_dict_1.drop('variant', axis=1)
df_new['variant'] = full_variant_wt

#save fitness output
gen_name = 'collapsed_fitness_'
file_save_name = "".join([gen_name,name_in])
df_new.to_csv(file_save_name,sep=',')
