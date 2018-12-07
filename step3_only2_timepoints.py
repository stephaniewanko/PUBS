#Author: Stephanie Wankowicz
#Date: 2018/10/16
'''
This code should take in three different csvs with the columns: barcode, count. This script should be run for each replicate, for each enviornment variable.

The input of the script should be:

python 20181015_normalization_lnvalues.py  csv_timepoint0 csv_timepoint12 csv_timepoint24 Representative_WT_Barcode_Seq output name

Where to get each
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
    print(ordered_df_tp3.head())
    tp1_tp2=ordered_df_tp1.append(pd.DataFrame(data=ordered_df_tp2),ignore_index=True)
    all_tp=tp1_tp2.append(pd.DataFrame(data=ordered_df_tp3),ignore_index=True)
    return(all_tp)


def ln_barcode_cal(df, tp1, tp2, wt1,wt2):
    df_zero_count=df.loc[df['count']==0]
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
    all_tp=tp1_tp2_tp3.append(pd.DataFrame(data=tp1_tp2_tp3),ignore_index=True)
    #print(df_tp1.head())
    return(all_tp)

#determining the number of 'dropout barcodes'
def barcode_dropout(df,tp1,tp2):
    barcode_tp0=barcode_all_tp.loc[(barcode_all_tp['Timepoint']==tp1)]
    barcode_tp2=barcode_all_tp.loc[(barcode_all_tp['Timepoint']==tp2)]
    #barcode_tp24=barcode_all_tp.loc[(barcode_all_tp['Timepoint']==tp3)]
    notinbarcode_2 = barcode_tp0[(barcode_tp0['barcode'].isin(barcode_tp2['barcode']) == False)]
    #notinbarcode_24 = barcode_tp0[(barcode_tp0['barcode'].isin(barcode_tp24['barcode']) == False)]
    print("The number of barcodes that dropped out at timepoint 12:")
    print(len(notinbarcode_2))
    #print("The number of barcodes that dropped out at timepoint 12:")
    #print(len(notinbarcode_24))
    n=1
    idx_2=len(barcode_tp2.index)
    print(idx_2)
    for i in notinbarcode_12['barcode']:
    		barcode_tp2.loc[idx_2+n,'barcode']=i
    		barcode_tp2.loc[idx_2+n,'count']=0
    		n+=1
    #idx_24=len(barcode_tp24.index)
    #n=1
    #for i in notinbarcode_24['barcode']:
#    		idx=len(barcode_tp24.index)
#    		barcode_tp24.loc[idx_24+n,'barcode']=i
#    		barcode_tp24.loc[idx_24+n,'count']=0
#    		n+=1
    barcode_tp0['count']+= 1
    barcode_tp2['count']+= 1
    all_tp=barcode_tp0.append(pd.DataFrame(data=barcode_tp2),ignore_index=True))
    print(head(all_tp))
    return(all_tp)
#    barcode_tp24['count']+= 1


if __name__ == "__main__":

#packages needed
    import pandas as pd
    import numpy as np
    import pickle
    import sys
##FILL IN THIS VALUE FROM THE WT SCRIPT
    WT_representative_barcode=sys.argv[4]
    barcode_values_tp0=pd.read_csv(sys.argv[1])
    barcode_values_tp2=pd.read_csv(sys.argv[2])
    print ('All inputs are in!')
#Importing pickle dictionary
    print ('Loading the Pickle!')
    f=open('2018_10_01_mismatch0.pkl', "rb")
    barcode_pickle = pickle.load(f)

#subset barcodes from CSVs to only those in the dictionary
    print ('Checking if variants are in each pickle')
    barcode_values_tp0=csv_pickle_check(barcode_values_tp0)
    barcode_values_tp2=csv_pickle_check(barcode_values_tp2)
    #print(barcode_values_tp0.head())

#cleaning input CSV and seperating Variant and Position
    print('Cleaning up the input files.')
    barcode_values_tp0=clean_input(barcode_values_tp0,0)
    barcode_values_tp2=clean_input(barcode_values_tp12,24)

#binding all of the barcodes together
    barcode_all_tp = barcode_values_tp0.append(pd.DataFrame(data = barcode_values_tp2), ignore_index=True)

#get the counts from the WT representative barcode (argv[4])
    WT_count_tp0=get_WT_barcode_count(WT_representative_barcode, 0)
    WT_count_tp2=get_WT_barcode_count(WT_representative_barcode, 24)
    print(WT_count_tp24)
    print(type(WT_count_tp24))
#look at barcode dropout
    print("Looking at the barcode dropout rate!")
    barcode_dropout(barcode_all_tp,0,24)
#calculate ln value
    print('Calculating the LN values!')
    barcode_all_tp=ln_barcode_cal(barcode_all_tp, 0,24,WT_count_tp0,WT_count_tp2)
#output csv
    output_file_name=sys.argv[5]
    barcode_all_tp.to_csv(output_file_name+'.csv')
