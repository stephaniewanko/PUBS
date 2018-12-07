#Author: Stephanie Wankowicz
#Date: 2018/10/16
'''
This code should take in three different csvs with the columns: barcode, count. This script should be run for each replicate, for each enviornment variable.

The input of the script should be:

python 20181015_normalization_lnvalues.py Representative_WT_Barcode_Seq csv_timepoint0 csv_timepoint12 csv_timepoint24

'''

def csv_pickle_check(df):
    df['Variant'] = df["barcode"].map(barcode_pickle)
    df = df.dropna(subset=['Variant'])
    df=df.reset_index(drop=True)
    return(df)

def clean_input(barcode_count,timepoint):
    barcode_count['Variant_Position']=barcode_count['Variant'].apply(pd.Series)[0]
    barcode_count['Variant_AA_Change']=barcode_count['Variant'].apply(pd.Series)[1]
    barcode_count['Timepoint']=timepoint
    #print(barcode_count.head())
    return(barcode_count)

def get_WT_barcode_count(WT_barcode,tp):
    WT_barcode_counts=barcode_all_tp.loc[barcode_all_tp['barcode']==WT_barcode]
    tp_WT_count1=((WT_barcode_counts.loc[WT_barcode_counts['Timepoint']==tp])['count']).iloc[0]
    tp_WT_count = tp_WT_count1.tolist()
    return(tp_WT_count)


def frequency_norm(df):
    tot_bar_count=df['count'].sum()
    print(tot_bar_count)
    df['Freq_Norm']=(df['count']/totalbarcode_count)
    return(df)

def median_norm(df,tp1,tp2,tp3):
    df_tp1=df.loc[(df['Timepoint'] == tp1)]
    med_df_tp1=df_tp1['count'].median()
    print(med_df_tp1)
    df_tp1['Med_Norm']=df_tp1['count']/med_df_tp1
    df_tp2=df.loc[(df['Timepoint'] == tp2)]
    med_df_tp2=df_tp2['count'].median()
    print(med_df_tp2)
    df_tp2['Med_Norm']=df_tp2['count']/med_df_tp2
    df_tp3=df.loc[(df['Timepoint'] == tp3)]
    med_df_tp3=df_tp3['count'].median()
    print(med_df_tp3)
    df_tp3['Med_Norm']=df_tp3['count']/med_df_tp3
    tp1_tp2=df_tp1.append(pd.DataFrame(data=df_tp2),ignore_index=True)
    all_tp=tp1_tp2.append(pd.DataFrame(data=df_tp3),ignore_index=True)
    return(all_tp)


def ln_barcode_cal(df, tp1, tp2, tp3, wt1,wt2, wt3):
    print(df.tail())
    #print(type(df.loc[1,'count']))
    df_zero_count=df.loc[df['count']==0]
    print('Zero length index: ')
    print(len(df_zero_count.index))
    print(df_zero_count.head())
    df_zero_count['ln(freq_BC/freq_med_WT)']=0
    df2=df.loc[df['count']>=1]
    df_tp1=df2.loc[df2['Timepoint'] == tp1]
    #print(df_tp1.loc[:,'count'])
    df_tp1['ln(freq_BC/freq_med_WT)']=np.log2(df_tp1.loc[:,'count']/wt1)
    df_tp2=df2.loc[(df2['Timepoint'] == tp2)]
    df_tp2['ln(freq_BC/freq_med_WT)']=np.log2(df_tp2.loc[:,'count']/wt2)
    df_tp3=df2.loc[(df2['Timepoint'] == tp3)]
    df_tp3['ln(freq_BC/freq_med_WT)']=np.log2(df_tp3.loc[:,'count']/wt3)
    tp1_tp2=df_tp1.append(pd.DataFrame(data=df_tp2),ignore_index=True)
    tp1_tp2_tp3=tp1_tp2.append(pd.DataFrame(data=df_tp3),ignore_index=True)
    all_tp=tp1_tp2_tp3.append(pd.DataFrame(data=df_zero_count),ignore_index=True)
    print(df_tp1.head())
    return(all_tp)

def quantile_normalize(df, tp1, tp2, tp3):
    df_tp1=df.loc[(df['Timepoint'] == tp1)]
    ordered_df_tp1=(df_tp1.sort_values(['count'], ascending=False)).reset_index(drop=True)
    print(len(ordered_df_tp1.index))
    df_tp2=df.loc[(df['Timepoint'] == tp2)]
    ordered_df_tp2=(df_tp2.sort_values(['count'], ascending=False)).reset_index(drop=True)
    print(len(ordered_df_tp2.index))
    df_tp3=df.loc[(df['Timepoint'] == tp3)]
    ordered_df_tp3=(df_tp3.sort_values(['count'], ascending=False)).reset_index(drop=True)
    print(len(ordered_df_tp3.index))
    for i in ordered_df_tp1.index:
     if i in (ordered_df_tp2.index & ordered_df_tp3.index):
       mean_value=((ordered_df_tp1.loc[i,'count']+ordered_df_tp2.loc[i,'count']+ordered_df_tp3.loc[i,'count'])/3)
     elif i in (ordered_df_tp2.index):
       mean_value=((ordered_df_tp1.loc[i,'count']+ordered_df_tp2.loc[i,'count'])/3)
     elif i in (ordered_df_tp3.index):
       mean_value=((ordered_df_tp1.loc[i,'count']+ordered_df_tp3.loc[i,'count'])/3)
     else:
       mean_value=((ordered_df_tp1.loc[i,'count'])/3)
     ordered_df_tp1.loc[i,'quant_norm']=mean_value
     if i in ordered_df_tp2.index:
       ordered_df_tp2.loc[i,'quant_norm']=mean_value
     if i in ordered_df_tp3.index:
       ordered_df_tp3.loc[i,'quant_norm']=mean_value
    #merge df back together
    tp1_tp2=ordered_df_tp1.append(pd.DataFrame(data=ordered_df_tp2),ignore_index=True)
    all_tp=tp1_tp2.append(pd.DataFrame(data=ordered_df_tp3),ignore_index=True)
    return(all_tp)

def barcode_dropout_df_mani(df,tp):
    df=df.drop(columns=['count', 'Timepoint'])
    df['count']=0
    df['Timepoint']=tp
    df=df[['barcode','count','Variant', 'Variant_Position','Variant_AA_Change', 'Timepoint']]
    return(df)

def barcode_dropout_tp0(df,tp1,tp2,tp3):
    barcode_tp0=df.loc[df['Timepoint']==tp1]
    barcode_tp12=df.loc[df['Timepoint']==tp2]
    barcode_tp24=df.loc[df['Timepoint']==tp3]
    barcode_tp12_list=barcode_tp12['barcode'].tolist()
    barcode_tp24_list=barcode_tp24['barcode'].tolist()
    barcode_tp0_list=barcode_tp0['barcode'].tolist()
    not_in_tp0_1=list(set(barcode_tp12_list)-set(barcode_tp0_list))
    not_in_tp0_2=list(set(barcode_tp24_list)-set(barcode_tp0_list))
    missing_tp0_df_1=(barcode_tp12.loc[barcode_tp12['barcode'].isin(not_in_tp0_1)]).reset_index()
    missing_tp0_df_1=barcode_dropout_df_mani(missing_tp0_df_1,0)
    missing_tp0_df_2=(barcode_tp24.loc[barcode_tp24['barcode'].isin(not_in_tp0_2)]).reset_index()
    missing_tp0_df_2=barcode_dropout_df_mani(missing_tp0_df_2,0)
    missing_tp0_df=missing_tp0_df_1.append(missing_tp0_df_2, ignore_index=True)
    missing_tp0_df=missing_tp0_df.drop_duplicates()
    df=df.append(missing_tp0_df, ignore_index=True)
    return(df)

def barcode_dropout_tp12_24(df,tp1,tp2,tp3):
    barcode_tp0=df.loc[df['Timepoint']==tp1]
    barcode_tp12=df.loc[df['Timepoint']==tp2]
    barcode_tp24=df.loc[df['Timepoint']==tp3]
    barcode_tp12_list=barcode_tp12['barcode'].tolist()
    barcode_tp24_list=barcode_tp24['barcode'].tolist()
    barcode_tp0_list=barcode_tp0['barcode'].tolist()
    missing_tp12=list(set(barcode_tp0_list)-set(barcode_tp12_list))
    missing_tp24=list(set(barcode_tp0_list)-set(barcode_tp24_list))
    missing_tp12_df=(barcode_tp0.loc[barcode_tp0['barcode'].isin(missing_tp12)]).reset_index()
    missing_tp12_df=barcode_dropout_df_mani(missing_tp12_df,12)
    missing_tp24_df=(barcode_tp0.loc[barcode_tp0['barcode'].isin(missing_tp24)]).reset_index()
    missing_tp24_df=barcode_dropout_df_mani(missing_tp24_df,24)
    missing=missing_tp24_df.append(missing_tp12_df, ignore_index=True)
    df=df.append(missing, ignore_index=True)
    return(df)

if __name__ == "__main__":

#packages needed
    import pandas as pd
    import numpy as np
    import pickle
    import sys
##Importing CSVs
    WT_representative_barcode=sys.argv[4]
    barcode_values_tp0=pd.read_csv(sys.argv[1], engine='python')
    barcode_values_tp12=pd.read_csv(sys.argv[2], engine='python')
    barcode_values_tp24=pd.read_csv(sys.argv[3], engine='python')
    print ('All inputs are in!')
#Importing pickle dictionary
    print ('Loading the Pickle!')
    f=open('2018_10_01_mismatch0.pkl', "rb")
    barcode_pickle = pickle.load(f)


#subset barcodes from CSVs to only those in the dictionary
    print ('Checking if variants are in each pickle')
    barcode_values_tp0=csv_pickle_check(barcode_values_tp0)
    barcode_values_tp12=csv_pickle_check(barcode_values_tp12)
    barcode_values_tp24=csv_pickle_check(barcode_values_tp24)
    #print(barcode_values_tp0.head())

#cleaning input CSV and seperating Variant and Position
    print('Cleaning up the input files.')
    barcode_values_tp0=clean_input(barcode_values_tp0,0)
    barcode_values_tp12=clean_input(barcode_values_tp12,12)
    barcode_values_tp24=clean_input(barcode_values_tp24,24)

#binding all of the barcodes together
    barcode_0_12 = barcode_values_tp0.append(pd.DataFrame(data = barcode_values_tp12), ignore_index=True)
    barcode_all_tp=barcode_0_12.append(pd.DataFrame(data=barcode_values_tp24),ignore_index=True)

#get the counts from the WT representative barcode (argv[4])
    WT_count_tp0=get_WT_barcode_count(WT_representative_barcode, 0)
    WT_count_tp12=get_WT_barcode_count(WT_representative_barcode, 12)
    WT_count_tp24=get_WT_barcode_count(WT_representative_barcode, 24)
    #print(WT_count_tp24)
    #print(type(WT_count_tp24))
#look at barcode dropout
    print("Looking at the barcode dropout rate!")
    barcode_all_tp=barcode_dropout_tp0(barcode_all_tp,0,12,24)
    barcode_all_tp=barcode_dropout_tp12_24(barcode_all_tp,0,12,24)
#calculate ln value
    print('Calculating the LN values!')
    barcode_all_tp=ln_barcode_cal(barcode_all_tp, 0, 12, 24,WT_count_tp0,WT_count_tp12,WT_count_tp24)
#output csv
    output_file_name=sys.argv[5]
    barcode_all_tp.to_csv(output_file_name+'step_3.csv',index=False)
