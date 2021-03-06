#Author: Stephanie Wankowicz
#Date: 10/21/2018
#script to (1)calculate the alternative fitness values (2) compare the alternative fitness value to the slope calculated fitness values
'''
The inputs for this script:
1)Regular Fitness Value (ouptut of script name: )
2)Collapse barcode CSV (output of script: fitness_preprocessing2.py)
3)WT Representative Barcode; string input; (output for script name:)
4)Output name
'''
def alt_fit_value(count1,count2,wt1,wt2):
    #get count from wt sequence
    if count1==0 or count2==0:
        return(0)
    else:
        afvalue=np.log(count2/wt2)-np.log(count1/wt1)
        return(afvalue)

def get_WT_barcode_count(WT_barcode,tp):
    WT_barcode_counts=barcode_all_tp.loc[barcode_all_tp['barcode']==WT_barcode]
    tp_WT_count1=((WT_barcode_counts.loc[WT_barcode_counts['Timepoint']==tp])['count']).iloc[0]
    tp_WT_count = tp_WT_count1.tolist()
    return(tp_WT_count)

if __name__ == "__main__":

#packages needed
    import pandas as pd
    import numpy as np
    import pickle
    import sys
    import matplotlib as plt
    import re
    #import matplotlib.pyplot as plt
    plt.use('Agg')

#import CSVs
    fitness_values=pd.read_csv(sys.argv[1])
    barcode_all_tp=pd.read_csv(sys.argv[2])
    WT_representative_barcode=sys.argv[3]
    output = str(sys.argv[4])
    print ('All inputs are loaded!')

#get the counts from the WT representative barcode (argv[4])
    print('Getting the representative WT barcode count!')
    WT_count_tp0=get_WT_barcode_count(WT_representative_barcode, 0)
    WT_count_tp12=get_WT_barcode_count(WT_representative_barcode, 12)
    WT_count_tp24=get_WT_barcode_count(WT_representative_barcode, 24)

#calculate alterntative fitness values
    print('Calculating the alternative fitenss value.')
    df_fitness_values=pd.DataFrame(columns=('Barcode', 'Count_TP0', 'Count_TP12', 'Count_TP24', 'Variant', 'Variant_Position', 'Variant_AA_Change'))
    n=0
    barcode_tp0=barcode_all_tp.loc[(barcode_all_tp['Timepoint']==0)]
    barcode_tp12=barcode_all_tp.loc[(barcode_all_tp['Timepoint']==12)]
    barcode_tp24=barcode_all_tp.loc[(barcode_all_tp['Timepoint']==24)]
    for i in barcode_all_tp['barcode'].unique():
      df_fitness_values.loc[n,'Barcode']=i
      df_fitness_values.loc[n,'Count_TP0']=(barcode_tp0.loc[(barcode_tp0['barcode']==i)])['count'].iloc[0]
      df_fitness_values.loc[n,'Variant']=(barcode_tp0.loc[(barcode_tp0['barcode']==i)])['Variant'].iloc[0]
      df_fitness_values.loc[n,'Variant_Position']=(barcode_tp0.loc[(barcode_tp0['barcode']==i)])['Variant_Position'].iloc[0]
      df_fitness_values.loc[n,'Variant_AA_Change']=(barcode_tp0.loc[(barcode_tp0['barcode']==i)])['Variant_AA_Change'].iloc[0]
      df_fitness_values.loc[n,'Count_TP12']=(barcode_tp12.loc[(barcode_tp12['barcode']==i)])['count'].iloc[0]
      df_fitness_values.loc[n,'Count_TP24']=(barcode_tp24.loc[(barcode_tp24['barcode']==i)])['count'].iloc[0]
      n+=1

    df_fitness_values['alt_fit_value_0_12'] = df_fitness_values.apply(lambda x: alt_fit_value(x['Count_TP12'], x['Count_TP0'],WT_count_tp0,WT_count_tp12), axis=1)
    df_fitness_values['alt_fit_value_12_24'] = df_fitness_values.apply(lambda x: alt_fit_value(x['Count_TP24'], x['Count_TP12'],WT_count_tp12,WT_count_tp24), axis=1)
    df_fitness_values['Variant2']=df_fitness_values['Variant_Position'].map(str)+df_fitness_values['Variant_AA_Change']
    fitness_values['Variant'] = fitness_values['variant'].map(lambda x: str(x)[-2:])
    fitness_values['fitness']=fitness_values[['fitness']].apply(pd.to_numeric)
    #df_fitness_values.to_csv('df_fitness_values.csv', ignore_index=True)

#collapsing barcode
    print('Collapsing Barcodes and finding Median Alternative Fitness Values.')
    alt_fit_value_collapse=pd.DataFrame()
    n=0
    for i in df_fitness_values['Variant2']:
        alt_fit_value_collapse.loc[n,'Variant']=i
        alt_fit_value_collapse.loc[n,'Variant_Position']=(df_fitness_values.loc[(df_fitness_values['Variant2']==i)])['Variant_Position'].iloc[0]
        alt_fit_value_collapse.loc[n,'Variant_AA_Change']=(df_fitness_values.loc[(df_fitness_values['Variant2']==i)])['Variant_AA_Change'].iloc[0]
        alt_fit_value_collapse.loc[n,'Alt_Fit_Value_0_12']=np.median(df_fitness_values.loc[(df_fitness_values['Variant2']==i)]['alt_fit_value_0_12'])
        alt_fit_value_collapse.loc[n,'Alt_Fit_Value_12_24']=np.median(df_fitness_values.loc[(df_fitness_values['Variant2']==i)]['alt_fit_value_12_24'])
        alt_fit_value_collapse.loc[n,'Alt_Fit_Value_12_24_SD']=np.std(df_fitness_values.loc[(df_fitness_values['Variant2']==i)]['alt_fit_value_12_24'])
        alt_fit_value_collapse.loc[n,'Alt_Fit_Value_0_12_SD']=np.std(df_fitness_values.loc[(df_fitness_values['Variant2']==i)]['alt_fit_value_0_12'])
        n+=1
    #print(type(fitness_values.loc[1,'fitness']))
    alt_fit_value_collapse['Variant']= alt_fit_value_collapse['Variant'].replace('.(WT)', 'WT', regex=True)
    print('Alt Fitness Value DF:')
    print(alt_fit_value_collapse.head(10))
    print('Fitness Value DF:')
    print(fitness_values.head())
#Compare Fitness values and output values
    #merge alt and fitness value df
    print('Creating plots to show the difference between fitness values.')
    both_fit_value=alt_fit_value_collapse.merge(fitness_values,how='outer',on='Variant')
    both_fit_value.to_csv(output+'_step6_output.csv', index=False)
#Figure 1
    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    ax.scatter(both_fit_value['Alt_Fit_Value_0_12'],both_fit_value['Alt_Fit_Value_12_24'])
    x_limit=both_fit_value['Alt_Fit_Value_0_12'].max()+0.5
    y_limit=both_fit_value['Alt_Fit_Value_12_24'].max()+0.5
    ax.set_xlim(0,x_limit)
    ax.set_ylim(0,y_limit)
    fig1.savefig('Alt_Fit_Value1vAlt_Fit_Value1.pdf')

    fig2=plt.figure()
    ax=fig2.add_subplot(111)
    ax.scatter(both_fit_value['Alt_Fit_Value_0_12'],both_fit_value['fitness'])
    x_limit=both_fit_value['Alt_Fit_Value_0_12'].max()+0.5
    y_limit=both_fit_value['fitness'].max()+0.5
    ax.set_xlim(0,x_limit)
    ax.set_ylim(0,y_limit)
    fig2.savefig('OverallFitnessvAlt_Fit_Value1.pdf')

    fig3=plt.figure()
    ax=fig3.add_subplot(111)
    ax.scatter(both_fit_value['Alt_Fit_Value_0_12'],both_fit_value['fitness'])
    x_limit=both_fit_value['Alt_Fit_Value_0_12'].max()+0.5
    y_limit=both_fit_value['fitness'].max()+0.5
    ax.set_xlim(0,x_limit)
    ax.set_ylim(0,y_limit)
    fig3.savefig('OverallFitnessvAlt_Fit_Value2.pdf')
