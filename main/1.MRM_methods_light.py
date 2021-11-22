import pandas as pd
import numpy as np
import os
import sys
import re

#input url
Peptide_Method_File_url = sys.argv[1] #.xlsx
aaIndex_url = sys.argv[2] #.xlsx/csv
formula_url = sys.argv[3] #.xlsx/csv
outputFile_url = sys.argv[4]
Time = sys.argv[5]
DoPrecursorCharge4 = float(sys.argv[6])
Q1MzMinimum = float(sys.argv[7])
Q1MzMaximum = float(sys.argv[8])
Q3MzMinimum = float(sys.argv[9])
Q3MzMaximum = float(sys.argv[10])
DefaultFile_url_exist = sys.argv[11] 
if DefaultFile_url_exist == 'non-existent':
    DefaultFile_url = 'None'
elif DefaultFile_url_exist == 'existent':
    DefaultFile_url = sys.argv[12] 


#建立output出來的資料夾
if 'methods' not in os.listdir(outputFile_url):
    os.mkdir('%s/methods_light'%outputFile_url)
    
#建立 1.Peptide Name + Sequence資料 2.Methods資料
PeptideFile = pd.read_excel(Peptide_Method_File_url,sheet_name='Peptide list') 
PeptideFile.columns = [i.lower() for i in PeptideFile.columns] #資料的column names都要變小寫
MethodFile =pd.read_excel(Peptide_Method_File_url,sheet_name='Ramping step')
MethodFile.columns = [i.lower() for i in MethodFile.columns] #資料的column names都要變小寫


#建立 1.amino acid:Mass資料 2.Mass-to-charge Ratio公式資料
if aaIndex_url.split('.')[-1] == 'xlsx':
    aaIndex = pd.read_excel(aaIndex_url)
elif aaIndex_url.split('.')[-1] == 'csv':
    aaIndex = pd.read_csv(aaIndex_url)

aa_string = '|'.join(aaIndex['aa'])
aa_string = aa_string.replace('[','\[')
aa_string = aa_string.replace(']','\]')
aaIndex = dict(zip(aaIndex['aa'],aaIndex['MW']))

if formula_url.split('.')[-1] == 'xlsx':
    formula = pd.read_excel(formula_url)
elif formula_url.split('.')[-1] == 'csv':
    formula = pd.read_csv(formula_url)
formula = dict(zip(formula['Peptide'],formula['Formula']))



#建立 1.DefaultFile
if DefaultFile_url != 'None':
    if DefaultFile_url.split('.')[-1] == 'xlsx':
        DefaultFile = pd.read_excel(DefaultFile_url)
    elif DefaultFile_url.split('.')[-1] == 'csv':
        DefaultFile = pd.read_csv(DefaultFile_url)
    
#################下面開始處理計算#################
#計算full peptide質量 (Dict)
aaMass_dict = dict()
for i in PeptideFile['peptide sequence']:
    aa_list = re.findall(aa_string,i)
    aaMass_dict[i] = sum([aaIndex[j] for j in aa_list])


#處理OutputFile
Names = []
yPep_dict = dict()
bPep_dict = dict()
uniprot_dict = dict(zip(PeptideFile['peptide sequence'],PeptideFile['uniprot']))
for pepseq in PeptideFile['peptide sequence']:
    if pepseq[-1] == ']':
        pepseq_del_tail = pepseq[:list(re.finditer('\[',pepseq))[-1].span()[0]]
    else:
        pepseq_del_tail = pepseq
    '''    
    #切出b-peptide
    b_peptide = []
    position = 2 #b2開始
    while True:
        if position < len(pepseq_del_tail):
            if pepseq_del_tail[0:position][-1] == '[':
                while True:
                    position += 1
                    if pepseq_del_tail[0:position][-1] == ']':
                        b_peptide.append(pepseq_del_tail[0:position])
                        position += 1
                        break
            else:
                b_peptide.append(pepseq_del_tail[0:position])
                position += 1
        else:
            break
    '''        
    #切出b-peptide
    b_aa_element = []
    b_aa = ''
    for ps in pepseq:
        b_aa += ps
        if ('[' in b_aa)&(']' not in b_aa):
            continue
        elif (']' in b_aa)&(b_aa[-1] != ']'):
            b_aa_element.append(b_aa[:-1])
            b_aa = b_aa[-1]      
        elif (len(b_aa)==2)&(b_aa[-1]!='['):
            b_aa_element.append(b_aa[:-1])
            b_aa = b_aa[-1]
    b_aa_element.append(b_aa)
    b_peptide = []
    for position in range(2,len(b_aa_element)):
        b_peptide.append(''.join(b_aa_element[:position]))
    #切出y-peptide
    y_aa_element = []
    y_aa = ''
    for ps in pepseq:
        y_aa += ps
        if ('[' in y_aa)&(']' not in y_aa):
            continue
        elif (']' in y_aa)&(y_aa[-1] != ']'):
            y_aa_element.append(y_aa[:-1])
            y_aa = y_aa[-1]      
        elif (len(y_aa)==2)&(y_aa[-1]!='['):
            y_aa_element.append(y_aa[:-1])
            y_aa = y_aa[-1]
    y_aa_element.append(y_aa)
    y_peptide = []
    for position in range(1,len(y_aa_element)-1):
        y_peptide.append(''.join(y_aa_element[position:]))


    #計算y/b-peptide質量
    for i in y_peptide:
        aa_list = re.findall(aa_string,i)
        yPep_dict[i] = sum([aaIndex[j] for j in aa_list])
    for i in b_peptide:
        aa_list = re.findall(aa_string,i)
        bPep_dict[i] = sum([aaIndex[j] for j in aa_list])

    #配對階層
    PreCharge = ['2','3']
    ProductCharge = ['1','2']
    for i in PreCharge:
        name = uniprot_dict[pepseq] + '_' + pepseq + '_' + i + '_Q1' + '+'*int(i)
        for j in y_peptide:
            aa_list = re.findall(aa_string,j)
            j_rm_modification = [aa for aa in aa_list if aa[0]!='[']          
            name_y = name + '_y' + str(len(j_rm_modification)) + '_' + j
            for k in ProductCharge:
                Names.append(name_y + '_+' + k + '_' + 'y' + '+'*int(k))
        for j in b_peptide:
            aa_list = re.findall(aa_string,j)
            j_rm_modification = [aa for aa in aa_list if aa[0]!='[']         
            name_b = name + '_b' + str(len(j_rm_modification)) + '_' + j
            for k in ProductCharge:
                Names.append(name_b + '_+' + k + '_' + 'b' + '+'*int(k))
OutputFile = pd.DataFrame([i.split('_') for i in Names])
OutputFile.columns = ['uniprot','species','Q1 peptide sequence','precursor charge','Q1 charge','split name','Q3 peptide sequence',
                          'product charge','Q3 charge']

#計算Q1和Q3
Q1 = []
Q3 = []
for i in OutputFile.index:
    mass = aaMass_dict[OutputFile['Q1 peptide sequence'][i]]
    Q1.append(eval(formula[OutputFile['Q1 charge'][i]]))
    if 'y' in OutputFile['split name'][i]:
        mass = yPep_dict[OutputFile['Q3 peptide sequence'][i]]
    elif 'b' in OutputFile['split name'][i]:
        mass = bPep_dict[OutputFile['Q3 peptide sequence'][i]]
    Q3.append(eval(formula[OutputFile['Q3 charge'][i]]))
OutputFile['Q1'] = Q1
OutputFile['Q3'] = Q3








#把 Q1+++<1250 的製造Q1++++
for i in OutputFile.index:
    if (OutputFile['Q1 charge'][i]=='Q1+++')&(OutputFile['Q1'][i]>=DoPrecursorCharge4):
        OutputFile = OutputFile.append(OutputFile.loc[i],ignore_index=True)
        OutputFile.at[OutputFile.index[-1],'Q1 charge'] = 'Q1++++'
        OutputFile.at[OutputFile.index[-1],'precursor charge'] = '4'
        mass = aaMass_dict[OutputFile['Q1 peptide sequence'][OutputFile.index[-1]]]
        OutputFile.at[OutputFile.index[-1],'Q1'] = eval(formula[OutputFile.at[OutputFile.index[-1],'Q1 charge']])
 

#保留 Q1++和Q1+++ 300 < Q1 < 1250
OutputFile = OutputFile.drop(OutputFile[(OutputFile['Q1 charge']!='Q1++++')&((OutputFile['Q1']<=Q1MzMinimum)|(OutputFile['Q1']>=Q1MzMaximum))].index)
OutputFile = OutputFile.reset_index(drop=True)
#保留 300 < Q3 < 1250
OutputFile = OutputFile[(OutputFile['Q3']>Q3MzMinimum) & (OutputFile['Q3']<Q3MzMaximum)]
#SORT 
OutputFile = OutputFile.sort_values(by=['Q1 peptide sequence','precursor charge','split name', 'product charge'])
OutputFile = OutputFile.reset_index(drop=True)
OutputFile.index = list(range(len(OutputFile.index)))



#Name
Name = []
for i in OutputFile.index:
    if OutputFile['product charge'][i] == '+1':
        Name.append('%s_%s.%s.+%s%s.light'%(OutputFile['uniprot'][i],OutputFile['species'][i],OutputFile['Q1 peptide sequence'][i],
                                           OutputFile['precursor charge'][i],OutputFile['split name'][i]))
    else:
        Name.append('%s_%s.%s.+%s%s%s.light'%(OutputFile['uniprot'][i],OutputFile['species'][i],OutputFile['Q1 peptide sequence'][i],
                                             OutputFile['precursor charge'][i],OutputFile['split name'][i],OutputFile['product charge'][i]))
OutputFile['Name'] = Name
 
#name_default 
name_default = []
for i in OutputFile.index:
    if OutputFile['product charge'][i] == '+1':
        name_default.append('%s.+%s%s'%(OutputFile['Q1 peptide sequence'][i],
                                      OutputFile['precursor charge'][i],
                                      OutputFile['split name'][i]))
    else:
        name_default.append('%s.+%s%s%s'%(OutputFile['Q1 peptide sequence'][i],
                                        OutputFile['precursor charge'][i],
                                        OutputFile['split name'][i],
                                        OutputFile['product charge'][i]))
OutputFile['name_default'] = name_default
    
#計算Overlap
Q1_unique = []
[Q1_unique.append(x) for x in OutputFile['Q1'] if x not in Q1_unique]

Overlap_Mass_list = []
for q1_mass in Q1_unique:
    OutputFile2 = OutputFile[OutputFile['Q1']==q1_mass]
    for q3 in OutputFile2.index:
        diff = abs(OutputFile2['Q3'] - OutputFile2['Q3'][q3]) <= 0.1
        
        if sum(diff) > 1:
            Overlap_Mass_list.append(list(diff[diff==True].index))
Overlap_Mass_set = []
[Overlap_Mass_set.append(x) for x in Overlap_Mass_list if x not in Overlap_Mass_set]

OutputFile['Mass overlap'] = ''
label = 1
multi_overlap = 1
for i in Overlap_Mass_set:
    for j in i:
        if OutputFile['Mass overlap'][j] == '':
            OutputFile['Mass overlap'][j] = label
        else:
            OutputFile['Mass overlap'][j] = 'multi-group overlap(%s)'%multi_overlap
    multi_overlap +=1
    label += 1

    
#處理default數值
if DefaultFile_url != 'None':
    DefaultFile.columns = [x.lower() for x in DefaultFile.columns]
    DefaultFile['name_default'] = ['%s.%s'%(DefaultFile['name'][x].split('.')[1],
                                            DefaultFile['name'][x].split('.')[2]) for x in range(len(DefaultFile['name']))]

    DP_default = []
    EP_default = []
    CE_default = []
    CXP_default = []
    for i in OutputFile['name_default']:
        if len(DefaultFile[DefaultFile['name_default']==i]) != 0:
            DP_default.append(list(DefaultFile[DefaultFile['name_default']==i]['dp'])[0])
            EP_default.append(list(DefaultFile[DefaultFile['name_default']==i]['ep'])[0])
            CE_default.append(list(DefaultFile[DefaultFile['name_default']==i]['ce'])[0])
            CXP_default.append(list(DefaultFile[DefaultFile['name_default']==i]['cxp'])[0])
        else:
            if len(MethodFile['dp_default']) != 0:
                DP_default.append(MethodFile['dp_default'][0])
            else:
                DP_default.append('NaN')
            if len(MethodFile['ep_default']) != 0:
                EP_default.append(MethodFile['ep_default'][0])
            else:
                EP_default.append('NaN')
            if len(MethodFile['ce_default']) != 0:
                CE_default.append(MethodFile['ce_default'][0])
            else:
                CE_default.append('NaN')
            if len(MethodFile['cxp_default']) != 0:
                CXP_default.append(MethodFile['cxp_default'][0])
            else:
                CXP_default.append('NaN')
    OutputFile['dp_default'] = DP_default
    OutputFile['ep_default'] = EP_default
    OutputFile['ce_default'] = CE_default
    OutputFile['cxp_default'] = CXP_default
    
elif DefaultFile_url == 'None':
    OutputFile['dp_default'] = MethodFile['dp_default'][0]
    OutputFile['ep_default'] = MethodFile['ep_default'][0]
    OutputFile['ce_default'] = MethodFile['ce_default'][0]
    OutputFile['cxp_default'] = MethodFile['cxp_default'][0]
    
#固定參數
Parameters = ['dp','ep','ce','cxp']
Sep_columns = [i + '_seperate' for i in Parameters]
Ramping_columns = [i + '_ramping' for i in Parameters]
Default_columns = [i + '_default' for i in Parameters]


#處理Seperate
#如果某column中的cell都不為'Na'，就要處理
for sepcol in Sep_columns:
    if sum(pd.isna(MethodFile[sepcol])) != len(MethodFile[sepcol]):
        for j in [i for i in MethodFile[sepcol] if not pd.isna(i)]:
            OutputFile2 = OutputFile.copy()
            OutputFile2['Time'] = Time
            OutputFile2['%s'%sepcol.split('_')[0].upper()] = j
            if DefaultFile_url == 'None':
                Default_columns2 = [x for x in Default_columns if sepcol.split('_')[0]+'_default' != x]
                for defcol in Default_columns2:
                    OutputFile2['%s'%defcol.split('_')[0].upper()] = MethodFile[defcol][0]
            else:
                Default_columns2 = [x for x in Default_columns  if sepcol.split('_')[0]+'_default' != x]
                for defcol in Default_columns2:
                    Parameter_content = []
                    for k in OutputFile2['name_default']:
                        Parameter_content.append(list(OutputFile2[OutputFile2['name_default'] == k][defcol])[0])
                    OutputFile2['%s'%defcol.split('_')[0].upper()] = Parameter_content
  
            #儲存資料
            output_columns = ['Q1','Q3','Time','Name','DP','EP','CE','CXP','Mass overlap']
            
            if 'seperated' not in os.listdir('%s/methods_light'%outputFile_url):
                os.mkdir('%s/methods_light/seperated'%outputFile_url)
            if DefaultFile_url == 'None':
                OutputFile2[output_columns].to_csv('%s/methods_light/seperated/seperated_%s%d.csv'%(outputFile_url,sepcol.split('_')[0].upper(),j),index=False)
            else:
                OutputFile2[output_columns].to_csv('%s/methods_light/seperated/seperated_%s%d(best_default=true).csv'%(outputFile_url,sepcol.split('_')[0].upper(),j),index=False)
                
                
                
#處理Ramping
#如果某column中的cell都不為'Na'，就要處理
for rampcol in Ramping_columns:
    if sum(pd.isna(MethodFile[rampcol])) != len(MethodFile[rampcol]):
        #展開OutputFile2
        OutputFile2 = OutputFile.copy()
        OutputFile2['Name'] = [str(OutputFile2['Q1'][i])+'_'+str(OutputFile2['Q3'][i])+'_'+OutputFile2['Name'][i]+
                               '_'+str(OutputFile2['Mass overlap'][i])+
                               '_'+str(OutputFile2['dp_default'][i])+
                               '_'+str(OutputFile2['ep_default'][i])+
                               '_'+str(OutputFile2['ce_default'][i])+
                               '_'+str(OutputFile2['cxp_default'][i])
                               for i in range(len(OutputFile2['Q1']))]
        Names = []
        for name in OutputFile2['Name']:
            Ramp_center = MethodFile[rampcol][0]
            Ramp_interval = MethodFile[rampcol][1]
            Ramp_Points = MethodFile[rampcol][2]
            Ramp_list = [Ramp_center]
            for i in range(int((Ramp_Points-1)/2)):
                Ramp_list.append(Ramp_center-Ramp_interval*(i+1))
            Ramp_list.reverse()
            for i in range(int((Ramp_Points-1)/2)):
                Ramp_list.append(Ramp_center+Ramp_interval*(i+1))

            for ramp in Ramp_list:
                if len(str(max(Ramp_list))) == 2:
                    Names.append('%s_%s%02d_%d'%(name,rampcol.split('_')[0].upper(),ramp,ramp))
                elif len(str(max(Ramp_list))) == 3:
                    Names.append('%s_%s%03d_%d'%(name,rampcol.split('_')[0].upper(),ramp,ramp))
                else:
                    Names.append('%s_%s%d_%d'%(name,rampcol.split('_')[0].upper(),ramp,ramp))
        OutputFile2 = pd.DataFrame([i.split('_') for i in Names])
        OutputFile2.columns = ['Q1','Q3','uniprot','species','Mass overlap','dp_default','ep_default','ce_default',
                               'cxp_default','ramping',rampcol.split('_')[0].upper()]
        Default_columns2 = [x for x in Default_columns if rampcol.split('_')[0]+'_default' != x]
        for defcol in Default_columns2:
            OutputFile2['%s'%defcol.split('_')[0].upper()] = OutputFile2[defcol]
        OutputFile2['Time'] = Time
        OutputFile2['Name'] = [OutputFile2['uniprot'][x]+'_'+OutputFile2['species'][x]+'_'+OutputFile2['ramping'][x] for x in range(len(OutputFile2['uniprot']))]


        #儲存資料
        output_columns = ['Q1','Q3','Time','Name','DP','EP','CE','CXP','Mass overlap']
        
        if 'ramping' not in os.listdir('%s/methods_light'%outputFile_url):
            os.mkdir('%s/methods_light/ramping'%outputFile_url)
        if DefaultFile_url == 'None':
            OutputFile2[output_columns].to_csv('%s/methods_light/ramping/%s_ramping.csv'%(outputFile_url,rampcol.split('_')[0].upper()),index=False)
        else:
            OutputFile2[output_columns].to_csv('%s/methods_light/ramping/%s_ramping(best_default=true).csv'%(outputFile_url,rampcol.split('_')[0].upper()),index=False)