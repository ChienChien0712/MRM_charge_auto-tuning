####################################################################################################
import numpy as np
import pandas as pd
import statsmodels.api as sm
#Bezier
def FEvaluate_Bezier(xarr,yarr,i0,i1,i2,i3,t):
    pts0 = [xarr[i0],yarr[i0]]
    pts1 = [xarr[i1],yarr[i1]]
    pts2 = [xarr[i2],yarr[i2]]
    pts3 = [xarr[i3],yarr[i3]]
    
    d01 = np.sqrt((pts0[0]-pts1[0])**2 + (pts0[1]-pts1[1])**2)
    d12 = np.sqrt((pts1[0]-pts2[0])**2 + (pts1[1]-pts2[1])**2)
    d23 = np.sqrt((pts2[0]-pts3[0])**2 + (pts2[1]-pts3[1])**2) 
    d02 = np.sqrt((pts0[0]-pts2[0])**2 + (pts0[1]-pts2[1])**2)
    d13 = np.sqrt((pts1[0]-pts3[0])**2 + (pts1[1]-pts3[1])**2)
    
    #bz0
    bz0 = pts1
    
    #bz1
    if i0 != i1:
        f = 1/6
    else:
        f = 1/3
    bz1 = [pts1[0]+(pts2[0]-pts0[0])*f, pts1[1]+(pts2[1]-pts0[1])*f]
    
    #bz2
    if i2 != i3:
        f = 1/6
    else:
        f = 1/3
    bz2 = [pts2[0]+(pts1[0]-pts3[0])*f, pts2[1]+(pts1[1]-pts3[1])*f]
    
    #bz3
    bz3 = pts2
    
    #Bezier4
    x_bezier4 = (1-t)**3*bz0[0] + 3*((1-t)**2)*t*bz1[0] + 3*(1-t)*(t**2)*bz2[0] + t**3*bz3[0]
    y_bezier4 = (1-t)**3*bz0[1] + 3*((1-t)**2)*t*bz1[1] + 3*(1-t)*(t**2)*bz2[1] + t**3*bz3[1]
  
    return [x_bezier4,y_bezier4]
####################################################################################################
#Bacic SG smooth
def SGsmooth_coeff(w,k):
    X = []
    for i in range(k):
        X.append(list(np.arange(-w,w+1,1)**i))
    X = np.matrix(X).T
    B = (X*(X.T*X).I*X.T)
    return B[w,:].tolist()[0]

def SGsmooth_lag_coeff(w,k,lag):
    X = []
    for i in range(k):
        X.append(list(np.arange(-w,w+1,1)**i))
    X = np.matrix(X).T
    B = (X*(X.T*X).I*X.T)
    return B[w+lag,:].tolist()[0]

def SavitzkyGolaySmooth(series,SG_w,SG_k):
    coeff = np.array(SGsmooth_coeff(SG_w,SG_k))
    
    seriesRaw = series.copy()
    seriesSmooth = series.copy()
    tmp = []
    for i in range(SG_w,len(seriesRaw)-SG_w):
        sub_series = []
        for j in np.arange(-SG_w,SG_w+1,1):
            sub_series.append(seriesRaw[i+j])
        sub_series = np.array(sub_series)
        tmp.append(sum(sub_series*coeff))
    seriesSmooth[SG_w:len(seriesRaw)-SG_w] = tmp
    return seriesSmooth
####################################################################################################
#Padding Methods of SG smooth
def SavitzkyGolaySmooth_RegressionPadding(series,SG_w,SG_k):
    coeff = np.array(SGsmooth_coeff(SG_w,SG_k))
    
    seriesRaw = series.copy()
    tmp = []
    for i in range(SG_w,len(seriesRaw)-SG_w):
        sub_series = []
        for j in np.arange(-SG_w,SG_w+1,1):
            sub_series.append(seriesRaw[i+j])
        sub_series = np.array(sub_series)
        if i == SG_w:
            for l in range(SG_w):
                tmp.append(sum(sub_series*np.array(SGsmooth_lag_coeff(SG_w,SG_k,-SG_w+l))))        
        tmp.append(sum(sub_series*coeff))
        if i == len(seriesRaw)-SG_w-1:
            for l in range(SG_w):
                tmp.append(sum(sub_series*np.array(SGsmooth_lag_coeff(SG_w,SG_k,1+l))))            
        
    return tmp

def SavitzkyGolaySmooth_VicinalPointPadding(series,SG_w,SG_k):
    coeff = np.array(SGsmooth_coeff(SG_w,SG_k))
    
    seriesRaw = series.copy()
    tmp = []
    for i in range(SG_w,len(seriesRaw)-SG_w):
        sub_series = []
        for j in np.arange(-SG_w,SG_w+1,1):
            sub_series.append(seriesRaw[i+j])
        sub_series = np.array(sub_series)
        tmp.append(sum(sub_series*coeff))

    tmp_pre = []
    for j in range(SG_w):
        sub_series = []
        for i in range(SG_w-j):
            sub_series.append(series[j+SG_w-i])
        for i in range(SG_w*2+1-len(sub_series)):
            sub_series.append(series[i])
        sub_series = np.array(sub_series)
        tmp_pre.append(sum(sub_series*coeff)) 
        
    tmp_post = []
    for j in list(reversed(range(1,SG_w+1))):
        sub_series = []
        for i in range(SG_w-j+1):
            sub_series.append(series[-j-SG_w+i])
        sub_series = np.array(series[-j-SG_w:] + list(reversed(sub_series)))
        
        tmp_post.append(sum(sub_series*coeff))
    seriesSmooth = list(tmp_pre)+list(tmp)+list(tmp_post)    

    return seriesSmooth

def SavitzkyGolaySmooth_ZeroPadding(series,SG_w,SG_k):
    coeff = np.array(SGsmooth_coeff(SG_w,SG_k))
    
    seriesRaw = series.copy()
    tmp = []
    for i in range(SG_w,len(seriesRaw)-SG_w):
        sub_series = []
        for j in np.arange(-SG_w,SG_w+1,1):
            sub_series.append(seriesRaw[i+j])
        sub_series = np.array(sub_series)
        tmp.append(sum(sub_series*coeff))

    tmp_pre = []
    for j in range(SG_w):
        sub_series = []
        for i in range(SG_w-j):
            sub_series.append(0)
        for i in range(SG_w*2+1-len(sub_series)):
            sub_series.append(series[i])
        sub_series = np.array(sub_series)
        tmp_pre.append(sum(sub_series*coeff)) 
        
    tmp_post = []
    for j in list(reversed(range(1,SG_w+1))):
        sub_series = []
        for i in range(SG_w-j+1):
            sub_series.append(0)
        sub_series = np.array(series[-j-SG_w:] + list(reversed(sub_series)))
        
        tmp_post.append(sum(sub_series*coeff))
        
    seriesSmooth = list(tmp_pre)+list(tmp)+list(tmp_post)    

    return seriesSmooth

def SavitzkyGolaySmooth_MinimumPadding(series,SG_w,SG_k):
    coeff = np.array(SGsmooth_coeff(SG_w,SG_k))
    
    seriesRaw = series.copy()
    tmp = []
    for i in range(SG_w,len(seriesRaw)-SG_w):
        sub_series = []
        for j in np.arange(-SG_w,SG_w+1,1):
            sub_series.append(seriesRaw[i+j])
        sub_series = np.array(sub_series)
        tmp.append(sum(sub_series*coeff))

    tmp_pre = []
    for j in range(SG_w):
        sub_series = []
        for i in range(SG_w-j):
            sub_series.append(min(seriesRaw))
        for i in range(SG_w*2+1-len(sub_series)):
            sub_series.append(series[i])
        sub_series = np.array(sub_series)
        tmp_pre.append(sum(sub_series*coeff)) 
        
    tmp_post = []
    for j in list(reversed(range(1,SG_w+1))):
        sub_series = []
        for i in range(SG_w-j+1):
            sub_series.append(min(seriesRaw))
        sub_series = np.array(series[-j-SG_w:] + list(reversed(sub_series)))
        
        tmp_post.append(sum(sub_series*coeff))
        
    seriesSmooth = list(tmp_pre)+list(tmp)+list(tmp_post)    

    return seriesSmooth

def SavitzkyGolaySmooth_MarginExtension(series,SG_w,SG_k):
    coeff = np.array(SGsmooth_coeff(SG_w,SG_k))
    
    seriesRaw = series.copy()
    tmp = []
    for i in range(SG_w,len(seriesRaw)-SG_w):
        sub_series = []
        for j in np.arange(-SG_w,SG_w+1,1):
            sub_series.append(seriesRaw[i+j])
        sub_series = np.array(sub_series)
        tmp.append(sum(sub_series*coeff))

    tmp_pre = []
    for j in range(SG_w):
        sub_series = []
        for i in range(SG_w-j):
            sub_series.append(seriesRaw[1])
        for i in range(SG_w*2+1-len(sub_series)):
            sub_series.append(series[i])
        sub_series = np.array(sub_series)
        tmp_pre.append(sum(sub_series*coeff)) 
        
    tmp_post = []
    for j in list(reversed(range(1,SG_w+1))):
        sub_series = []
        for i in range(SG_w-j+1):
            sub_series.append(seriesRaw[-1])
        sub_series = np.array(series[-j-SG_w:] + list(reversed(sub_series)))
        
        tmp_post.append(sum(sub_series*coeff))
        
    seriesSmooth = list(tmp_pre)+list(tmp)+list(tmp_post)    

    return seriesSmooth
####################################################################################################
#Regression
def regression3to6_param(pol_level,x,y):
    if pol_level == 6:
        X = np.column_stack((x, np.array(pd.array(x)**2,dtype='int64'), np.array(pd.array(x)**3,dtype='int64'),np.array(pd.array(x)**4,dtype='int64'),np.array(pd.array(x)**5,dtype='int64'),np.array(pd.array(x)**6,dtype='int64')))
        X = sm.add_constant(X)
        Y = y
        model = sm.OLS(Y,X)
        results = model.fit()
        g,f,e,d,c,b,a = results.params[0],results.params[1],results.params[2],results.params[3],results.params[4],results.params[5],results.params[6]
        return [a,b,c,d,e,f,g]
    elif pol_level == 5:
        X = np.column_stack((x, np.array(pd.array(x)**2,dtype='int64'), np.array(pd.array(x)**3,dtype='int64'),np.array(pd.array(x)**4,dtype='int64'),np.array(pd.array(x)**5,dtype='int64')))
        X = sm.add_constant(X)
        Y = y
        model = sm.OLS(Y,X)
        results = model.fit()
        f,e,d,c,b,a = results.params[0],results.params[1],results.params[2],results.params[3],results.params[4],results.params[5]
        return [a,b,c,d,e,f]
    elif pol_level == 4:
        X = np.column_stack((x, np.array(pd.array(x)**2,dtype='int64'), np.array(pd.array(x)**3,dtype='int64'),np.array(pd.array(x)**4,dtype='int64')))
        X = sm.add_constant(X)
        Y = y
        model = sm.OLS(Y,X)
        results = model.fit()
        e,d,c,b,a = results.params[0],results.params[1],results.params[2],results.params[3],results.params[4]
        return [a,b,c,d,e]
    elif pol_level == 3:
        X = np.column_stack((x, np.array(pd.array(x)**2,dtype='int64'), np.array(pd.array(x)**3,dtype='int64')))
        X = sm.add_constant(X)
        Y = y
        model = sm.OLS(Y,X)
        results = model.fit()
        d,c,b,a = results.params[0],results.params[1],results.params[2],results.params[3]
        return [a,b,c,d]