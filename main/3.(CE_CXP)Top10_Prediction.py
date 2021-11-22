import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from Prediction_Formula import *
import statsmodels.api as sm
import sys
plt.rcParams['figure.figsize'] = (9.0, 5.0)
plt.rcParams['figure.dpi'] = 80

filename = sys.argv[1]
outputfile_url = sys.argv[2]
csvfilename = os.listdir('%s'%filename)[0].split('.')[0]
plot_pic = sys.argv[3].split('/')[-1]


#Read file
f1 = pd.read_csv('%s/%s.csv'%(filename,csvfilename))
f1['Precursor Charge'] = list(map(str,f1['Precursor Charge']))
f1['Product Charge'] = list(map(str,f1['Product Charge']))
f1['Replicate Name'] = list(map(str,f1['Replicate Name']))
f1['Name'] = f1[['Precursor Charge','Fragment Ion','Product Charge']].apply(lambda x: '.'.join(x), axis = 1) 
f1['Name'] = f1[['Peptide Sequence','Name']].apply(lambda x: '_'.join(x), axis = 1) 
f1['Area_Background'] = f1['Area'] + f1['Background']
#RepName
f1['RepName'] = f1[['Name','Replicate Name']].apply(lambda x: ';'.join(x), axis = 1)
f1.head()


###Setting Parameters
if 'CE' in csvfilename:
    CHARGE = 'CE'
    CHARGE_RANGE = [5,10,15,20,25,30,35,40,45,50,55]
elif 'CXP' in csvfilename:
    CHARGE = 'CXP'   
    CHARGE_RANGE = [5,10,15,20,25,30,35,40,45,50,55]
if 'ramping' in csvfilename:
    ramp_sep_type = 'ramping'
elif 'separate' in csvfilename:
    ramp_sep_type = 'separate'
else:
    ramp_sep_type = 'ramping'
    
Name = []
for i in f1['RepName']:
    if i not in Name:
        Name.append(i)    
################################Bezier################################
################################Bezier################################
################################Bezier################################
if csvfilename not in os.listdir('%s'%outputfile_url):
    os.mkdir('%s/%s'%(outputfile_url,csvfilename))
if 'Bezier' not in os.listdir('%s/%s'%(outputfile_url,csvfilename)):
    os.mkdir('%s/%s/Bezier'%(outputfile_url,csvfilename))
if 'Cubic_Bezier' not in os.listdir('%s/%s/Bezier'%(outputfile_url,csvfilename)):
    os.mkdir('%s/%s/Bezier/Cubic_Bezier'%(outputfile_url,csvfilename))

CEramping_forcsv = []
CEbest_forcsv = []
CEbest_marginRm_forcsv = []
Name_forcsv = []    
PredArea_forcsv = []

for name in Name:
    f2 = f1[f1['RepName']==name]
    y_Area=list(f2['Area_Background'])
    if 'des' in name:
        y_Area.reverse()
    
    x_ReplicateName = CHARGE_RANGE
    interval = x_ReplicateName[1]- x_ReplicateName[0]

    if len(x_ReplicateName) != len(y_Area):
        print('length error: %s'%name)
        continue

    X_plot = []
    Y_plot = []
    for t in np.arange(0,1+0.1/interval,0.1/interval):
        Bezier_point = FEvaluate_Bezier(x_ReplicateName[:4],y_Area[:4],0,0,1,2,t)
        X_plot.append(Bezier_point[0])
        Y_plot.append(Bezier_point[1])

    for i in range(len(x_ReplicateName)-3):
        for t in np.arange(0,1+0.1/interval,0.1/interval):
            Bezier_point = FEvaluate_Bezier(x_ReplicateName[i:4+i],y_Area[i:4+i],0,1,2,3,t)
            X_plot.append(Bezier_point[0])
            Y_plot.append(Bezier_point[1])

    for t in np.arange(0,1+0.1/interval,0.1/interval):
        Bezier_point = FEvaluate_Bezier(x_ReplicateName[-4:],y_Area[-4:],1,2,3,3,t)
        X_plot.append(Bezier_point[0])
        Y_plot.append(Bezier_point[1])

    Name_forcsv.append(name.split(';')[0])
    CEramping_forcsv.append(list(f2['Replicate Name'])[0])
    x_max = X_plot[Y_plot.index(max(Y_plot))]
    y_max = max(Y_plot)
    if y_max != 0:
        CEbest_forcsv.append(x_max)
    else:
        CEbest_forcsv.append('NA')


    X_plot_rmMargin = list(np.array(X_plot)[(np.array(X_plot)>=x_ReplicateName[1])&(np.array(X_plot)<=x_ReplicateName[-2])])
    Y_plot_rmMargin = list(np.array(Y_plot)[(np.array(X_plot)>=x_ReplicateName[1])&(np.array(X_plot)<=x_ReplicateName[-2])])   
    x_max_rmMargin = X_plot_rmMargin[Y_plot_rmMargin.index(max(Y_plot_rmMargin))]
    y_max_rmMargin = max(Y_plot_rmMargin)        
    if y_max_rmMargin != 0:
        CEbest_marginRm_forcsv.append(x_max_rmMargin)
    else:
        CEbest_marginRm_forcsv.append('NA')    
        
        
    PredArea_forcsv.append(y_max_rmMargin)

    if plot_pic == 'TRUE':
        plt.scatter(x_ReplicateName,y_Area)
        plt.plot(X_plot,Y_plot,label='Cubic Bézier curves')
        plt.xticks(x_ReplicateName,x_ReplicateName,fontsize=12)
        plt.yticks(fontsize=12)
        if y_max_rmMargin != 0:
            plt.scatter(x_max_rmMargin,y_max_rmMargin,color='red')
            plt.annotate('%s=%.1f'%(CHARGE,x_max_rmMargin), (x_max_rmMargin,y_max_rmMargin),color='red')
        plt.grid(linestyle='--',color='gray')
        plt.xlabel('%s'%CHARGE,fontsize=20)
        plt.ylabel('Area + Background',fontsize=20)
        plt.title('%s(Bezier)\n%s'%(csvfilename,name))
        plt.legend()
        plt.savefig('%s/%s/Bezier/Cubic_Bezier/%s.png'%(outputfile_url,csvfilename,name),dpi=80,bbox_inches='tight')
        plt.clf()

output = pd.DataFrame(list(zip(Name_forcsv,CEbest_forcsv,CEbest_marginRm_forcsv,CEramping_forcsv,PredArea_forcsv)))
output.columns = ['Name','Best %s'%CHARGE,'Best %s (marginRm)'%CHARGE,'%s%s'%(CHARGE,ramp_sep_type),'Predicted Area']
output.to_csv('%s/%s/Bezier/Cubic_Bezier.csv'%(outputfile_url,csvfilename),index=False)   
################################Bezier################################
################################Bezier################################
################################Bezier################################


################################Regression################################
################################Regression################################
################################Regression################################
if 'Regression' not in os.listdir('%s/%s'%(outputfile_url,csvfilename)):
    os.mkdir('%s/%s/Regression'%(outputfile_url,csvfilename))
    

def FitPolynomial(f1,points_number,pol_level,margineRm=False):
    Name_forcsv = []
    CEbest_forcsv = []
    CEramping_forcsv = []
    PredArea_forcsv = []
    
    Name = []
    [Name.append(x) for x in f1['RepName'] if x not in Name]
    
    if '%sPoints%sDegree'%(points_number,pol_level) not in os.listdir('%s/%s/Regression'%(outputfile_url,csvfilename)):
        os.mkdir('%s/%s/Regression/%sPoints%sDegree'%(outputfile_url,csvfilename,points_number,pol_level))
    
    for n in Name:
        Name_forcsv.append(n.split(';')[0])   
        f2 = f1[(f1['RepName']==n)]
        f2 = f2.sort_values(by=['Opt Step'])
        x_axis = CHARGE_RANGE
        y_axis = list(f2['Area_Background'])
        if 'des' in n:
            y_axis.reverse()    
        if len(y_axis) != len(x_axis):
            print('length error: %s'%n)
            continue        
        
            
        CEramping_forcsv.append(list(f2['Replicate Name'])[0])        
        #fit model
        #polynomial和point自訂
        Parameter_set = []
        x_forline_set = []
        y_forline_set = []
        x_perpoint_set = []
        y_perpoint_set = []
        start = 0
        end = points_number
        while True:
            x_axis_shift = x_axis[start:end]
            y_axis_shift = y_axis[start:end]
            Parameter_set.append(regression3to6_param(pol_level,x_axis_shift,y_axis_shift))

            x_forline = np.arange(min(x_axis_shift),max(x_axis_shift),0.1)
            x_forline_set.append(x_forline)      
            parameters_number_minus1 = pol_level
            y_forline = regression3to6_param(pol_level,x_axis_shift,y_axis_shift)[parameters_number_minus1]
            while True:
                parameters_number_minus1 -= 1
                y_forline += regression3to6_param(pol_level,x_axis_shift,y_axis_shift)[parameters_number_minus1]*x_forline**(pol_level-parameters_number_minus1)
                if parameters_number_minus1 == 0:
                    break
            y_forline_set.append(y_forline)
                    
            #找最大y值所用的
            if margineRm == False:
                x_perpoint = np.arange(min(x_axis_shift),max(x_axis_shift)+0.1,0.1)
            elif margineRm == True:
                x_perpoint = np.arange(x_axis_shift[1],x_axis_shift[-2]+0.1,0.1)
            parameters_number_minus1 = pol_level
            y_perpoint = regression3to6_param(pol_level,x_axis_shift,y_axis_shift)[parameters_number_minus1]
            while True:
                parameters_number_minus1 -= 1
                y_perpoint += regression3to6_param(pol_level,x_axis_shift,y_axis_shift)[parameters_number_minus1]*x_perpoint**(pol_level-parameters_number_minus1)
                if parameters_number_minus1 == 0:
                    break
            x_perpoint_set.append(x_perpoint)
            y_perpoint_set.append(y_perpoint)

            start += 1
            end += 1
            if end > len(x_axis):
                break
        #Find peak
        x_max_global = 0
        y_max_global = 0
        for j in range(len(x_perpoint_set)):
            xy_perpoint = dict(zip(x_perpoint_set[j],y_perpoint_set[j]))
            x_max = max(xy_perpoint, key=xy_perpoint.get)
            y_max = xy_perpoint[x_max]       
            if y_max >= y_max_global:
                x_max_global = x_max
                y_max_global = y_max
        if max(y_axis)==0:
            CEbest_forcsv.append('NA')
        else:    
            CEbest_forcsv.append('%.1f'%x_max_global)
        
        #預測的Area
        PredArea_forcsv.append(y_max_global)
        
        if plot_pic == 'TRUE':
            #plot
            fig = plt.figure()
            ax1 = fig.add_subplot(1,1,1)
            ax1.xaxis.set_ticks_position('bottom')
            ax1.yaxis.set_ticks_position('left')
            ax1.scatter(x_axis,y_axis,color='black')
            plt.xticks(x_axis,x_axis,rotation=0,fontsize=12)
            plt.grid(linestyle='--', which='major',color='gray')
            plt.xlabel('%s'%CHARGE)
            plt.ylabel('Area+Background')
            plt.title('%s(%sP,%sD)\n%s'%(csvfilename,points_number,pol_level,n))

            #plot 回歸線
            for j in range(len(x_forline_set)):
                ax1.plot(x_forline_set[j],y_forline_set[j],linestyle='--')
            #plot最高點
            if max(y_axis)!=0:
                ax1.scatter(x_max_global,y_max_global,color='red')
                ax1.annotate('%s=%.1f'%(CHARGE,x_max_global), (x_max_global,y_max_global),color='red')

            plt.savefig('%s/%s/Regression/%sPoints%sDegree/%s.png'%(outputfile_url,csvfilename,points_number,pol_level,n),dpi=100,bbox_inches='tight')
            plt.close(fig)                        

    output = pd.DataFrame(list(zip(Name_forcsv,CEbest_forcsv,CEramping_forcsv,PredArea_forcsv)))
    output.columns = ['Name','Best %s'%CHARGE,'%s%s'%(CHARGE,ramp_sep_type),'Predicted Area']
    output.to_csv('%s/%s/Regression/%sP%sD.csv'%(outputfile_url,csvfilename,points_number,pol_level),index=False)
FitPolynomial(f1,6,3,margineRm=True)
FitPolynomial(f1,6,4,margineRm=True)
FitPolynomial(f1,7,3,margineRm=True)
FitPolynomial(f1,7,4,margineRm=True)
FitPolynomial(f1,7,5,margineRm=True)
FitPolynomial(f1,8,3,margineRm=True)
FitPolynomial(f1,8,4,margineRm=True)
FitPolynomial(f1,8,5,margineRm=True)
FitPolynomial(f1,8,6,margineRm=True)
FitPolynomial(f1,9,3,margineRm=True)
FitPolynomial(f1,9,4,margineRm=True)
FitPolynomial(f1,9,5,margineRm=True)
FitPolynomial(f1,9,6,margineRm=True)
FitPolynomial(f1,10,3,margineRm=True)
FitPolynomial(f1,10,4,margineRm=True)
FitPolynomial(f1,10,5,margineRm=True)
FitPolynomial(f1,10,6,margineRm=True)
FitPolynomial(f1,11,3,margineRm=True)
FitPolynomial(f1,11,4,margineRm=True)
FitPolynomial(f1,11,5,margineRm=True)
FitPolynomial(f1,11,6,margineRm=True)
################################Regression################################
################################Regression################################
################################Regression################################



################################SG-smooth################################
################################SG-smooth################################
################################SG-smooth################################
if 'S-Gsmoothing_w=5_p=2' not in os.listdir('%s/%s'%(outputfile_url,csvfilename)):
    os.mkdir('%s/%s/S-Gsmoothing_w=5_p=2'%(outputfile_url,csvfilename))
m = 2
k = 3

Methods = ['Original','RegressionPadding','VicinalPointPadding','ZeroPadding','MinimumPadding','MarginExtension']

for method in Methods:
    if method not in os.listdir('%s/%s/S-Gsmoothing_w=5_p=2'%(outputfile_url,csvfilename)):
        os.mkdir('%s/%s/S-Gsmoothing_w=5_p=2/%s'%(outputfile_url,csvfilename,method))
        
    CEramping_forcsv = []
    CEbest_forcsv = []
    CEbest_marginRm_forcsv = []
    Name_forcsv = []    
    PredArea_forcsv = []
    
    for name in Name:
        f2 = f1[f1['RepName']==name]
        y_Area=list(f2['Area_Background']) 
        if 'des' in name:
            y_Area.reverse()
            
        x_ReplicateName = CHARGE_RANGE

        if len(x_ReplicateName) != len(y_Area):
            print('length error: %s'%name)
            continue        

        if method == 'Original':
            y_smooth = SavitzkyGolaySmooth(y_Area,m,k)      
        elif method == 'RegressionPadding':
            y_smooth = SavitzkyGolaySmooth_RegressionPadding(y_Area,m,k) 
        elif method == 'VicinalPointPadding':
            y_smooth = SavitzkyGolaySmooth_VicinalPointPadding(y_Area,m,k) 
        elif method == 'ZeroPadding':
            y_smooth = SavitzkyGolaySmooth_ZeroPadding(y_Area,m,k) 
        elif method == 'MinimumPadding':
            y_smooth = SavitzkyGolaySmooth_MinimumPadding(y_Area,m,k) 
        elif method == 'MarginExtension':
            y_smooth = SavitzkyGolaySmooth_MarginExtension(y_Area,m,k) 


        CEramping_forcsv.append(list(f2['Replicate Name'])[0])
        Name_forcsv.append(name.split(';')[0])

        #marginRm = False
        if max(y_smooth) == 0:
            CEbest_forcsv.append('NA')
        else:
            if y_smooth.count(max(y_smooth)) >= 1:
                max_index = []
                pos = 0
                while True:
                    try:
                        max_index.append(y_smooth.index(max(y_smooth),pos))
                        pos = y_smooth.index(max(y_smooth),pos)+1
                    except:
                        break
                CEbest_forcsv.append(np.mean(np.array(x_ReplicateName)[max_index]))
            else:        
                CEbest_forcsv.append(x_ReplicateName[y_smooth.index(max(y_smooth))])

        #marginRm = True
        x_Rm = x_ReplicateName[1:-1]
        y_Rm = y_smooth[1:-1]
        if max(y_Rm) == 0:
            CEbest_marginRm_forcsv.append('NA')
        else:
            if y_Rm.count(max(y_Rm)) >= 1:
                max_index = []
                pos = 0
                while True:
                    try:
                        max_index.append(y_Rm.index(max(y_Rm),pos))
                        pos = y_Rm.index(max(y_Rm),pos)+1
                    except:
                        break
                CEbest_marginRm_forcsv.append(np.mean(np.array(x_Rm)[max_index]))
            else:        
                CEbest_marginRm_forcsv.append(x_Rm[y_Rm.index(max(y_Rm))])    

        PredArea_forcsv.append(max(y_Rm))        

        if plot_pic == 'TRUE':
            plt.bar(x_ReplicateName,y_Area, edgecolor='white',label='raw',width=3,zorder=-1)
            plt.plot(x_ReplicateName,y_smooth,color='red',label='S-G smoothing')
            plt.xticks(x_ReplicateName,x_ReplicateName,fontsize=12)
            plt.yticks(fontsize=12)
            plt.xlabel('%s'%CHARGE,fontsize=20)
            plt.ylabel('Area + Background',fontsize=20)
            plt.title('%s(%s)\n%s'%(csvfilename,method,name))
            plt.legend()
            plt.grid(linestyle='--',color='gray')
            plt.savefig('%s/%s/S-Gsmoothing_w=5_p=2/%s/%s.png'%(outputfile_url,csvfilename,method,name),dpi=80,bbox_inches='tight') 
            plt.clf()

    output = pd.DataFrame(list(zip(Name_forcsv,CEbest_forcsv,CEbest_marginRm_forcsv,CEramping_forcsv,PredArea_forcsv)))
    output.columns = ['Name','Best %s'%CHARGE,'Best %s (marginRm)'%CHARGE,'%s%s'%(CHARGE,ramp_sep_type),'Predicted Area']
    output.to_csv('%s/%s/S-Gsmoothing_w=5_p=2/%s.csv'%(outputfile_url,csvfilename,method),index=False)  
################################SG-smooth################################
################################SG-smooth################################
################################SG-smooth################################


################################SG-smooth+Bezier################################
################################SG-smooth+Bezier################################
################################SG-smooth+Bezier################################
if 'S-Gsmoothing+Bezier' not in os.listdir('%s/%s'%(outputfile_url,csvfilename)):
    os.mkdir('%s/%s/S-Gsmoothing+Bezier'%(outputfile_url,csvfilename))
m = 2
k = 3

Methods = ['Original','RegressionPadding','VicinalPointPadding','ZeroPadding','MinimumPadding','MarginExtension']

for method in Methods:
    if '%s+Bezier'%method not in os.listdir('%s/%s/S-Gsmoothing+Bezier'%(outputfile_url,csvfilename)):
        os.mkdir('%s/%s/S-Gsmoothing+Bezier/%s+Bezier'%(outputfile_url,csvfilename,method))
        
    CEramping_forcsv = []
    CEbest_forcsv = []
    CEbest_marginRm_forcsv = []
    Name_forcsv = []   
    PredArea_forcsv = []
    
    for name in Name:
        f2 = f1[f1['RepName']==name]
        y_Area=list(f2['Area_Background']) 
        if 'des' in name:
            y_Area.reverse()
            
        x_ReplicateName = CHARGE_RANGE

        if len(x_ReplicateName) != len(y_Area):
            print('length error: %s'%name)
            continue


        if method == 'Original':
            y_smooth = SavitzkyGolaySmooth(y_Area,m,k)      
        elif method == 'RegressionPadding':
            y_smooth = SavitzkyGolaySmooth_RegressionPadding(y_Area,m,k) 
        elif method == 'VicinalPointPadding':
            y_smooth = SavitzkyGolaySmooth_VicinalPointPadding(y_Area,m,k) 
        elif method == 'ZeroPadding':
            y_smooth = SavitzkyGolaySmooth_ZeroPadding(y_Area,m,k) 
        elif method == 'MinimumPadding':
            y_smooth = SavitzkyGolaySmooth_MinimumPadding(y_Area,m,k) 
        elif method == 'MarginExtension':
            y_smooth = SavitzkyGolaySmooth_MarginExtension(y_Area,m,k) 

        CEramping_forcsv.append(list(f2['Replicate Name'])[0])
        Name_forcsv.append(name.split(';')[0])

        ########################################
        X_plot = []
        Y_plot = []
        interval = x_ReplicateName[1]- x_ReplicateName[0]
        for t in np.arange(0,1+0.1/interval,0.1/interval):
            Bezier_point = FEvaluate_Bezier(x_ReplicateName[:4],y_smooth[:4],0,0,1,2,t)
            X_plot.append(Bezier_point[0])
            Y_plot.append(Bezier_point[1])

        for i in range(len(x_ReplicateName)-3):
            for t in np.arange(0,1+0.1/interval,0.1/interval):
                Bezier_point = FEvaluate_Bezier(x_ReplicateName[i:4+i],y_smooth[i:4+i],0,1,2,3,t)
                X_plot.append(Bezier_point[0])
                Y_plot.append(Bezier_point[1])

        for t in np.arange(0,1+0.1/interval,0.1/interval):
            Bezier_point = FEvaluate_Bezier(x_ReplicateName[-4:],y_smooth[-4:],1,2,3,3,t)
            X_plot.append(Bezier_point[0])
            Y_plot.append(Bezier_point[1])


        x_max = X_plot[Y_plot.index(max(Y_plot))]
        y_max = max(Y_plot)
        if y_max != 0:
            CEbest_forcsv.append(x_max)
        else:
            CEbest_forcsv.append('NA')


        X_plot_rmMargin = list(np.array(X_plot)[(np.array(X_plot)>=x_ReplicateName[1])&(np.array(X_plot)<=x_ReplicateName[-2])])
        Y_plot_rmMargin = list(np.array(Y_plot)[(np.array(X_plot)>=x_ReplicateName[1])&(np.array(X_plot)<=x_ReplicateName[-2])])   
        x_max_rmMargin = X_plot_rmMargin[Y_plot_rmMargin.index(max(Y_plot_rmMargin))]
        y_max_rmMargin = max(Y_plot_rmMargin)        
        if y_max_rmMargin != 0:
            CEbest_marginRm_forcsv.append(x_max_rmMargin)
        else:
            CEbest_marginRm_forcsv.append('NA')    
            
        PredArea_forcsv.append(y_max_rmMargin)                   
        ##########################################


        if plot_pic == 'TRUE':
            plt.bar(x_ReplicateName,y_Area, edgecolor='white',label='raw',width=3,zorder=-1)
            plt.scatter(x_ReplicateName,y_smooth,color='black',label='S-G smoothing',zorder=1)
            plt.plot(X_plot,Y_plot,color='black',label='Bezier')
            plt.xticks(x_ReplicateName,x_ReplicateName,fontsize=12)
            plt.yticks(fontsize=12)
            plt.xlabel('%s'%CHARGE,fontsize=20)
            plt.ylabel('Area + Background',fontsize=20)
            plt.title('%s(%s+Bezier)\n%s'%(csvfilename,method,name))
            plt.legend()
            plt.grid(linestyle='--',color='gray')
            if y_max_rmMargin != 0:
                plt.scatter(x_max_rmMargin,y_max_rmMargin,color='red')
                plt.annotate('%s=%.1f'%(CHARGE,x_max_rmMargin), (x_max_rmMargin,y_max_rmMargin),color='red')            
            plt.savefig('%s/%s/S-Gsmoothing+Bezier/%s+Bezier/%s.png'%(outputfile_url,csvfilename,method,name),dpi=80,bbox_inches='tight') 
            plt.clf()

    output = pd.DataFrame(list(zip(Name_forcsv,CEbest_forcsv,CEbest_marginRm_forcsv,CEramping_forcsv,PredArea_forcsv)))
    output.columns = ['Name','Best %s'%CHARGE,'Best %s (marginRm)'%CHARGE,'%s%s'%(CHARGE,ramp_sep_type),'Predicted Area']
    output.to_csv('%s/%s/S-Gsmoothing+Bezier/%s+Bezier.csv'%(outputfile_url,csvfilename,method),index=False)          
################################SG-smooth+Bezier################################
################################SG-smooth+Bezier################################
################################SG-smooth+Bezier################################   


##########################分開Replicate Name##############################
Methods = [ i for i in os.listdir('%s/%s'%(outputfile_url,csvfilename)) if i[-4:] not in ['.png','xlsx','.csv']]
for method in Methods:
    datafile = os.listdir('%s/%s/%s'%(outputfile_url,csvfilename,method))
    datafile = [i for i in datafile if i[-4:]=='.csv']
    os.mkdir('%s/%s/%s/replicate_name_combination'%(outputfile_url,csvfilename,method))
    for df in datafile:
        f1 = pd.read_csv('%s/%s/%s/%s'%(outputfile_url,csvfilename,method,df))

        rep = np.unique(f1['%s%s'%(CHARGE,ramp_sep_type)])
        for i in rep:
            f2 = f1.loc[f1['%s%s'%(CHARGE,ramp_sep_type)] == i]
            if str(i)[0]=='(':
                f2.to_csv('%s/%s/%s/%s%s'%(outputfile_url,csvfilename,method,i,df),index=False)
            else:
                f2.to_csv('%s/%s/%s/(%s)%s'%(outputfile_url,csvfilename,method,i,df),index=False)
        os.replace('%s/%s/%s/%s'%(outputfile_url,csvfilename,method,df),'%s/%s/%s/replicate_name_combination/%s'%(outputfile_url,csvfilename,method,df))    
        
        
######################Top 10 #######################################
#單一Replication name過濾出top10
MethFold = [ i for i in os.listdir('%s/%s'%(outputfile_url,csvfilename)) if i[-4:] not in ['.png','xlsx','.csv']]
for mf in MethFold:
    Meth = os.listdir('%s/%s/%s'%(outputfile_url,csvfilename,mf))
    Meth = [i for i in Meth if i[-4:]=='.csv']
    os.mkdir('%s/%s/%s/non-Top10_results'%(outputfile_url,csvfilename,mf))
    for m in Meth:
        f1 = pd.read_csv('%s/%s/%s/%s'%(outputfile_url,csvfilename,mf,m))
        f1['Peptide Sequence'] = [i.split('_')[0] for i in f1['Name']]
        f1 = f1.sort_values(by=['Peptide Sequence', 'Predicted Area'], ascending=False)
        
        f_write = open('%s/%s/%s/(Top10)%s'%(outputfile_url,csvfilename,mf,m),'w')
        f_write.write(','.join(list(f1.columns))+',Rank\n')
        for pep in np.unique(f1['Peptide Sequence']):
            f2 = f1[f1['Peptide Sequence']==pep].reset_index(drop=True).loc[0:9]
            f2['Rank']= list(range(1,len(f2)+1))
            for i in list(range(len(f2))):
                f_write.write(','.join(list(map(str,f2.loc[i])))+'\n')            
        f_write.close()     
        os.replace('%s/%s/%s/%s'%(outputfile_url,csvfilename,mf,m),'%s/%s/%s/non-Top10_results/%s'%(outputfile_url,csvfilename,mf,m))    
        
        
        
        