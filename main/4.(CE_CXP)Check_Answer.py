import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import statsmodels.api as sm
import time
from scipy.stats import t
import sys
import shutil
plt.rcParams['figure.figsize'] = (5.0, 5.0)
plt.rcParams['figure.dpi'] = 100

#arguments
standard_folder = sys.argv[1]
prediction_folder = sys.argv[2]
DoDistance = sys.argv[3]
DoMSE = sys.argv[4]

#parameters
filename = os.listdir(prediction_folder)
filename = [i for i in filename if ('CE' in i)|('CXP' in i)]
filename = filename[0]

if 'CE' in filename:
    CHARGE = 'CE'
    CHARGE_RANGE = [5,10,15,20,25,30,35,40,45,50,55]
    CHARGE_SPAN = CHARGE_RANGE[1]-CHARGE_RANGE[0]
elif 'CXP' in filename:
    CHARGE = 'CXP'
    CHARGE_RANGE = [5,10,15,20,25,30,35,40,45,50,55]
    CHARGE_SPAN = CHARGE_RANGE[1]-CHARGE_RANGE[0]
if 'ramping' in filename:
    ramp_sep_type = 'ramping'
elif ('separate' in filename)|('seperate' in filename):
    ramp_sep_type = 'separate'
else:
    ramp_sep_type = 'ramping'

#Standard charge
standard_file = os.listdir(standard_folder)[0]
if standard_file.split('.')[-1] == 'xlsx':
    f1 = pd.read_excel('%s/%s'%(standard_folder,standard_file))
else:
    f1 = pd.read_csv('%s/%s'%(standard_folder,standard_file))
Name = list(f1['name'])
f1['name2'] = [i.split('.')[1]+'_'+i.split('.')[2][1]+'.'+i.split('.')[2][2:] for i in Name]
Name = []
for i in f1['name2']:
    if i[-2] == '+':
        Name.append(i.split('+')[0]+'.'+i.split('+')[1])
    else:
        Name.append(i+'.1')
f1['name2'] = Name       
Answer = dict(zip(f1['name2'],f1['CE']))
#Create methods
Methods = [ i for i in os.listdir('%s/%s'%(prediction_folder,filename)) if i[-4:] not in ['.png','xlsx','.csv']]

#####################Distance########################
if DoDistance == 'TRUE':
    for method in Methods:
        FileList = [i for i in os.listdir('%s/%s/%s'%(prediction_folder,filename,method)) if i[-3:]=='csv' ]
        FileList = [i for i in FileList if i[:4]!='(Top']
        PointNumber_list = []
        Sum_Difference_list = []
        Mean_Difference_list = []
        Method_list = []
        for fn in FileList:
            f2 = pd.read_csv('%s/%s/%s/%s'%(prediction_folder,filename,method,fn))
            if method in ['Regression','S-Gsmoothing+Regression']:
                Modeled = dict(zip(f2['Name']+';'+f2['%s%s'%(CHARGE,ramp_sep_type)],f2['Best %s'%CHARGE]))
            else:
                Modeled = dict(zip(f2['Name']+';'+f2['%s%s'%(CHARGE,ramp_sep_type)],f2['Best %s (marginRm)'%CHARGE]))
            Modeled_name = list(Modeled.keys())
            answer_value = []
            modeled_value = []
            for n in Modeled_name:
                if n.split(';')[0] in Answer:
                    answer_value.append(Answer[n.split(';')[0]])
                    modeled_value.append(Modeled[n])

            ColorPlotValue = []
            difference = 0
            for i in range(len(answer_value)):
                ColorPlotValue.append(np.ceil(abs(answer_value[i]-modeled_value[i])/CHARGE_SPAN)*CHARGE_SPAN)
                difference += abs(answer_value[i]-modeled_value[i])
            f3 = pd.DataFrame(list(zip(answer_value,modeled_value,ColorPlotValue)))
            f3.columns = ['answer_value','modeled_value','colorplotvalue']

            PointNumber_list.append(len(answer_value))
            Sum_Difference_list.append(difference)
            Mean_Difference_list.append(difference/len(answer_value))    
            Method_list.append(fn.split('.')[0])

            Diff = np.unique(f3['colorplotvalue'])
            for d in Diff:
                plt.scatter(f3[f3['colorplotvalue']==d]['modeled_value'],f3[f3['colorplotvalue']==d]['answer_value'],label='diff≤'+str(d))
            plt.title('%s\n'%filename + fn.split('.')[0],fontsize=16)
            plt.xlabel('Estimated %s'%CHARGE,fontsize=14)
            plt.ylabel('Standard %s'%CHARGE,fontsize=14)
            #plt.xlim(min(answer_value+modeled_value)-5,max(answer_value+modeled_value)+5)
            #plt.ylim(min(answer_value+modeled_value)-5,max(answer_value+modeled_value)+5)
            plt.xticks(CHARGE_RANGE)
            plt.yticks(CHARGE_RANGE)
            plt.grid(linestyle='--')
            plt.plot(CHARGE_RANGE,CHARGE_RANGE,color='black')
            plt.legend()
            plt.savefig('%s/%s/%s/(Diff)%s.png'%(prediction_folder,filename,method,fn.split('.')[0]))
            plt.clf()

        df = pd.DataFrame(list(zip(Method_list,Sum_Difference_list,PointNumber_list,Mean_Difference_list)))
        df.columns = ['Method','sum of difference','point number','mean of difference']
        df.to_csv('%s/%s/(Diff)%s_%s.csv'%(prediction_folder,filename,filename,method),index=False)
    
#####################################MSE#######################################
if DoMSE == 'TRUE':
    for method in Methods:
        FileList = [i for i in os.listdir('%s/%s/%s'%(prediction_folder,filename,method)) if i[-3:]=='csv']
        FileList = [i for i in FileList if i[:4]!='(Top']

        Method_list = []
        PointNumber_list = []
        R_list = []
        R_squared_list = []
        MSE_residuals = []
        for fn in FileList:
            f2 = pd.read_csv('%s/%s/%s/%s'%(prediction_folder,filename,method,fn))
            if method in ['Regression','S-Gsmoothing+Regression']:
                Modeled = dict(zip(f2['Name']+';'+f2['%s%s'%(CHARGE,ramp_sep_type)],f2['Best %s'%CHARGE]))
            else:
                Modeled = dict(zip(f2['Name']+';'+f2['%s%s'%(CHARGE,ramp_sep_type)],f2['Best %s (marginRm)'%CHARGE]))
            Modeled_name = list(Modeled.keys())
            answer_value = []
            modeled_value = []
            for n in Modeled_name:
                if n.split(';')[0] in Answer:
                    answer_value.append(Answer[n.split(';')[0]])
                    modeled_value.append(Modeled[n])

            X = sm.add_constant(modeled_value)
            Y = answer_value
            model = sm.OLS(Y,X)
            results = model.fit()

            Method_list.append(fn.split('.')[0])
            PointNumber_list.append(len(answer_value))
            R_squared_list.append(results.rsquared)
            R_list.append(np.sqrt(results.rsquared))
            MSE_residuals.append(results.mse_resid)

            x_line = np.arange(min(modeled_value),max(modeled_value),0.1)
            y_line = results.params[0] + results.params[1]*x_line
            r_squared = round(results.rsquared*100,4)
            plt.plot(x_line,y_line,color='red',label='$R^2$=%s%%'%r_squared)
            plt.scatter(modeled_value,answer_value,color='black')
            plt.xticks(CHARGE_RANGE)
            plt.yticks(CHARGE_RANGE)
            plt.grid(linestyle='--')
            #plt.xlim(min(answer_value+modeled_value)-5,max(answer_value+modeled_value)+5)
            #plt.ylim(min(answer_value+modeled_value)-5,max(answer_value+modeled_value)+5)
            plt.title('%s\n'%filename + fn.split('.')[0],fontsize=16)
            plt.xlabel('Estimated %s'%CHARGE,fontsize=14)
            plt.ylabel('Standard %s'%CHARGE,fontsize=14)
            plt.legend()

            plt.savefig('%s/%s/%s/(MSE)%s.png'%(prediction_folder,filename,method,fn.split('.')[0]))
            plt.clf()

        df = pd.DataFrame(list(zip(Method_list,PointNumber_list,R_squared_list,R_list,MSE_residuals)))
        df.columns = ['Method','point number','R-squared','R','MSE_residuals']
        df.to_csv('%s/%s/(MSE)%s_%s.csv'%(prediction_folder,filename,filename,method),index=False)    
    
################################合併####################################
#(Diff)
if DoDistance == 'TRUE':
    csvFileName = [i for i in os.listdir('%s/%s'%(prediction_folder,filename)) if i[0:6] == '(Diff)']
    Method_list = []
    Sum_Difference_list = []
    PointNumber_list = []
    Mean_Difference_list = []
    for csv in csvFileName:
        f = pd.read_csv('%s/%s/%s'%(prediction_folder,filename,csv))
        Method_list.extend(list(f['Method']))
        Sum_Difference_list.extend(list(f['sum of difference']))
        PointNumber_list.extend(list(f['point number']))
        Mean_Difference_list.extend(list(f['mean of difference']))
    df = pd.DataFrame(list(zip(Method_list,Sum_Difference_list,PointNumber_list,Mean_Difference_list)))
    df.columns = ['Method','sum of difference','point number','mean of difference']
    df.to_excel('%s/%s/(Diff,Comparison)%s.xlsx'%(prediction_folder,filename,filename),na_rep='NA',index=False)

#(MSE)
if DoMSE == 'TRUE':
    csvFileName = [i for i in os.listdir('%s/%s'%(prediction_folder,filename)) if i[0:5] == '(MSE)']
    Method_list = []
    PointNumber_list = []
    R_squared_list = []
    R_list = []
    MSE_residuals_list = []
    for csv in csvFileName:
        f = pd.read_csv('%s/%s/%s'%(prediction_folder,filename,csv))
        Method_list.extend(list(f['Method']))
        PointNumber_list.extend(list(f['point number']))
        R_squared_list.extend(list(f['R-squared']))
        R_list.extend(list(f['R']))
        MSE_residuals_list.extend(list(f['MSE_residuals']))
    df = pd.DataFrame(list(zip(Method_list,PointNumber_list,R_squared_list,R_list,MSE_residuals_list)))
    df.columns = ['Method','point number','R-squared','R','MSE_residuals']
    df.to_excel('%s/%s/(MSE,Comparison)%s.xlsx'%(prediction_folder,filename,filename),na_rep='NA',index=False)

#Delete csv file
csvFileName = [i for i in os.listdir('%s/%s'%(prediction_folder,filename)) if i[-3:] == 'csv']
for csv in csvFileName:
    os.remove('%s/%s/%s'%(prediction_folder,filename,csv))    
for method in Methods:
    empty_folder = os.listdir('%s/%s/%s'%(prediction_folder,filename,method))
    empty_folder = [i for i in empty_folder if (i[-3:]!='csv')&(i[-3:]!='png')]
    for ef in empty_folder:
        if len(os.listdir('%s/%s/%s/%s'%(prediction_folder,filename,method,ef))) == 0:
            shutil.rmtree('%s/%s/%s/%s'%(prediction_folder,filename,method,ef))