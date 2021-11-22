import numpy as np
import pandas as pd
import re
import sys
import os

#####Arguments setting#####
filename = sys.argv[1]
outputfile_url = sys.argv[2]
csvfilename = os.listdir(filename)[0].split('.')[0]
process_AreaMean =sys.argv[3]

#####Parameters setting#####
#ramp_sep (default:ramping)
if 'ramping' in csvfilename:
    ramp_sep = 'ramping'
elif 'separate' in csvfilename:
    ramp_sep = 'separate'
else:
    ramp_sep = 'ramping'

#####open file#####
if os.listdir(filename)[0].split('.')[-1] == 'csv':
    f1 = pd.read_csv('%s/%s.csv'%(filename,csvfilename))
    f1 = f1.dropna(subset=['Replicate Name'])
elif os.listdir(filename)[0].split('.')[-1] == 'xlsx':
    f1 = pd.read_excel('%s/%s.xlsx'%(filename,csvfilename))
    f1 = f1.dropna(subset=['Replicate Name'])
#####處理 replicate name#####
string = f1['Replicate Name']
string = ';'.join(string)
#先判斷 asc desc
regex_ad_1 = re.compile('desc|asc|des')
regex_ad_2 = re.compile('-[0-9]+|_[0-9]+')
match_ad_1 = regex_ad_1.findall(string)
match_ad_2 = regex_ad_2.findall(string)
#再判斷 -[0-9]
regex_09 = re.compile('-[0-9]+|_[0-9]+')
match_09 = regex_09.findall(string)
#其他電壓Optimization
regex_opt = re.compile('Reg|SG\+B')
match_opt = regex_opt.findall(string)
#剩下沒重複的
regex_noRep = re.compile('desc|asc|des|ramp')
match_noRep = regex_noRep.findall(string)
#改寫Replicate Name
if (len(match_ad_1) == len(f1['Replicate Name'])) & (len(match_ad_2) == len(f1['Replicate Name'])):
    f1['Replicate Name'] = [i+j for i,j in zip(match_ad_1,match_ad_2)]
elif len(match_09) == len(f1['Replicate Name']):
    f1['Replicate Name'] = [i[1:] for i in match_09]
elif len(match_opt) == len(f1['Replicate Name']):
    f1['Replicate Name'] = [i for i in match_opt]
else: #沒有重複的就不用寫出Area Mean
    f1['Replicate Name'] = [i for i in match_noRep]
#處理separate 的Opt Step
if ramp_sep == 'separate':
    OptStep = [0]
    for i in range(1,1+int((len(set(f1['SampleGroup']))-1)/2)):
        OptStep.append(i)
        OptStep.append(-i)
    OptStep.sort()
    OptStep_dict = dict(zip(sorted(list(set(f1['SampleGroup']))),OptStep))
    f1['Opt Step'] = [OptStep_dict[i] for i in f1['SampleGroup']]

#####寫出一個新的檔案_preprocessed#####
PeptideSequence = np.unique(f1['Peptide Sequence'])
#處理asc/desc的合併
if (len(match_ad_1) == len(f1['Replicate Name'])) & (len(match_ad_2) == len(f1['Replicate Name'])):
    f_write = open('%s/%s.csv'%(outputfile_url,csvfilename.split('/')[-1]),'w')
    f_write.write('Peptide Sequence,Precursor Charge,Fragment Ion,Product Charge,Replicate Name,Opt Step,Area,Background\n')
    f1['Replicate Name2'] = match_ad_1
    if process_AreaMean == "TRUE":
        for ps in PeptideSequence:
            f2 = f1[f1['Peptide Sequence']==ps]
            PrecursorCharge = np.unique(f2['Precursor Charge'])
            for pc1 in PrecursorCharge:
                f3 = f2[f2['Precursor Charge']==pc1]
                FragmentIon = np.unique(f3['Fragment Ion'])
                for fi in FragmentIon:
                    f4 = f3[f3['Fragment Ion']==fi]
                    ProductCharge = np.unique(f4['Product Charge'])
                    for pc2 in ProductCharge:
                        f5 = f4[f4['Product Charge']==pc2]
                        for rn in np.unique(match_ad_1):
                            f6 = f5[f5['Replicate Name2']==rn]
                            for os in [-5,-4,-3,-2,-1,0,1,2,3,4,5]:
                                f7 = f6[f6['Opt Step']==os]
                                Area = sum(f7['Area'])
                                Background = sum(f7['Background'])
                                f_write.write('%s,%s,%s,%s,(AreaMean)(%s),%s,%s,%s\n'%(ps,pc1,fi,pc2,rn,os,Area,Background))
        if len(np.unique(match_ad_1)) > 1:
            for ps in PeptideSequence:
                f2 = f1[f1['Peptide Sequence']==ps]
                PrecursorCharge = np.unique(f2['Precursor Charge'])
                for pc1 in PrecursorCharge:
                    f3 = f2[f2['Precursor Charge']==pc1]
                    FragmentIon = np.unique(f3['Fragment Ion'])
                    for fi in FragmentIon:
                        f4 = f3[f3['Fragment Ion']==fi]
                        ProductCharge = np.unique(f4['Product Charge'])
                        for pc2 in ProductCharge:
                            f5 = f4[f4['Product Charge']==pc2]
                            for os in [-5,-4,-3,-2,-1,0,1,2,3,4,5]:
                                Area = 0
                                Background = 0
                                for rn in np.unique(match_ad_1):
                                    if rn == 'asc':
                                        f6 = f5[(f5['Replicate Name2']==rn)&(f5['Opt Step']==-os)]
                                    else:
                                        f6 = f5[(f5['Replicate Name2']==rn)&(f5['Opt Step']==os)]
                                    Area += sum(f6['Area'])
                                    Background += sum(f6['Background'])
                                asc = np.unique(match_ad_1)[0]

                                f_write.write('%s,%s,%s,%s,(AreaMean)(%s),%s,%s,%s\n'%(ps,pc1,fi,pc2,'&'.join(np.unique(match_ad_1)),os,Area,Background))
    #把不重複的寫入
    columns = ['Peptide Sequence','Precursor Charge','Fragment Ion','Product Charge','Replicate Name','Opt Step','Area','Background']
    f2 = f1[columns]
    f2['Precursor Charge'] = list(map(str,f2['Precursor Charge']))
    f2['Product Charge'] = list(map(str,f2['Product Charge']))
    f2['Opt Step'] = list(map(int,f2['Opt Step']))
    f2['Opt Step'] = list(map(str,f2['Opt Step']))
    f2['Area'] = list(map(int,f2['Area']))
    f2['Area'] = list(map(str,f2['Area']))
    f2['Background'] = list(map(int,f2['Background']))
    f2['Background'] = list(map(str,f2['Background']))
    for i in f2.index:
        f_write.write(','.join(f2.loc[i])+'\n')
    f_write.close()    

#單純處理重複的
elif len(match_09) == len(f1['Replicate Name']):
    f_write = open('%s/%s.csv'%(outputfile_url,csvfilename.split('/')[-1]),'w')
    f_write.write('Peptide Sequence,Precursor Charge,Fragment Ion,Product Charge,Replicate Name,Opt Step,Area,Background\n')
    if process_AreaMean == "TRUE":
        for ps in PeptideSequence:
            f2 = f1[f1['Peptide Sequence']==ps]
            PrecursorCharge = np.unique(f2['Precursor Charge'])
            for pc1 in PrecursorCharge:
                f3 = f2[f2['Precursor Charge']==pc1]
                FragmentIon = np.unique(f3['Fragment Ion'])
                for fi in FragmentIon:
                    f4 = f3[f3['Fragment Ion']==fi]
                    ProductCharge = np.unique(f4['Product Charge'])
                    for pc2 in ProductCharge:
                        f5 = f4[f4['Product Charge']==pc2]
                        for os in [-5,-4,-3,-2,-1,0,1,2,3,4,5]:
                            Area = 0
                            Background = 0
                            for rn in np.unique(f1['Replicate Name']):
                                f6 = f5[(f5['Replicate Name']==rn)&(f5['Opt Step']==os)]
                                Area += sum(f6['Area'])
                                Background += sum(f6['Background'])
                            f_write.write('%s,%s,%s,%s,(AreaMean),%s,%s,%s\n'%(ps,pc1,fi,pc2,os,Area,Background))
    #把不重複的寫入
    columns = ['Peptide Sequence','Precursor Charge','Fragment Ion','Product Charge','Replicate Name','Opt Step','Area','Background']
    f2 = f1[columns]
    f2['Precursor Charge'] = list(map(str,f2['Precursor Charge']))
    f2['Product Charge'] = list(map(str,f2['Product Charge']))
    f2['Opt Step'] = list(map(int,f2['Opt Step']))
    f2['Opt Step'] = list(map(str,f2['Opt Step']))
    f2['Area'] = list(map(int,f2['Area']))
    f2['Area'] = list(map(str,f2['Area']))
    f2['Background'] = list(map(int,f2['Background']))
    f2['Background'] = list(map(str,f2['Background']))
    for i in f2.index:
        f_write.write(','.join(f2.loc[i])+'\n')
    f_write.close()
    
#不用處理的   
else:
    f_write = open('%s/%s.csv'%(outputfile_url,csvfilename.split('/')[-1]),'w')
    f_write.write('Peptide Sequence,Precursor Charge,Fragment Ion,Product Charge,Replicate Name,Opt Step,Area,Background\n')
    columns = ['Peptide Sequence','Precursor Charge','Fragment Ion','Product Charge','Replicate Name','Opt Step','Area','Background']
    f2 = f1[columns]
    f2['Precursor Charge'] = list(map(str,f2['Precursor Charge']))
    f2['Product Charge'] = list(map(str,f2['Product Charge']))
    f2['Opt Step'] = list(map(int,f2['Opt Step']))
    f2['Opt Step'] = list(map(str,f2['Opt Step']))
    f2['Area'] = list(map(int,f2['Area']))
    f2['Area'] = list(map(str,f2['Area']))
    f2['Background'] = list(map(int,f2['Background']))
    f2['Background'] = list(map(str,f2['Background']))
    for i in f2.index:
        f_write.write(','.join(f2.loc[i])+'\n')
    f_write.close()

