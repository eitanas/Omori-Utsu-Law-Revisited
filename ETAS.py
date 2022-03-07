#!/usr/bin/env python
# coding: utf-8

# Yossi Ashkenazy
# Eitan Asher
# Yongwen Zhang
# The source is https://github.com/zhangyongwen77/Yongwen-ETAS2-model
# The relevant paper is: http://arxiv.org/abs/2003.12539
#%%
import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
import math

bvalue=1

params = {"mu":0.2, 
         "A":6.26,
         "c":0.007,
         "alpha1":1.5,
         "alpha2":1.5,
         "P":1.13,
         "Q":2,
         "D":0.03,
         "gamma":0.48,
         "coss":200,
         "M0":3}
'''
params = {"mu":0.2, 
         "A":3.26,
         "c":0.007,
         "alpha1":2,
         "alpha2":1.4,
         "P":1.13,
         "Q":2,
         "D":0.03,
         "gamma":0.48,
         "coss":200,
         "M0":3}
'''

def ETAS_T(data, params, n, bvalue, M0):

    time=1*data.time
    mag=1*data.mag
    if (len(time)==0):
        time=math.inf 
        mag=0
        
    mag=mag-M0
  
    ti=max(time)
    
    # find the earthquake rate (number of events per day)
    Rmax=conditional_intensity(mag,time,ti,params)
  
    while ti<n:
        if (Rmax>0):
            tau=np.random.exponential(scale=1/Rmax)
        else:
            tau=math.inf
        ti=ti+tau
        rate=conditional_intensity(mag,time, ti, params)
        if (np.random.uniform()<=(rate/Rmax)):
        # select a magnitude
            new_mag=np.random.exponential(scale=(1/(bvalue*np.log(10))))
            time=np.r_[time,ti]
            mag=np.r_[mag,new_mag]
      
            ti=max(time)
            Rmax=conditional_intensity(mag,time, ti, params)
        else:
            Rmax=rate
  
    a,b = time, mag+M0
    
    return(a,b)

def conditional_intensity(mag,time, T, params):
    mu=params["mu"]
    A=params["A"]
    alpha1=params["alpha1"]
    alpha2=params["alpha2"]
    CC=params["c"]
    P=params["P"] 
    coss=params["coss"]
            
    N_time=len(time)
    N_mag=len(mag)
    if (N_time>0):
        if (N_time>coss):
            ci1=A*sum(np.exp(alpha1*mag[N_mag-coss:N_mag])*(1+(T-time[N_mag-coss:N_mag])/CC)**(-P))
            ci2=A*sum(np.exp(alpha2*mag[0:N_mag-coss])*(1+(T-time[0:N_mag-coss])/CC)**(-P))
            ci=ci1+ci2
        else:
            ci=A*sum(np.exp(alpha1*mag)*(1+(T-time)/CC)**(-P))
    else:
        ci=0
        
    ci=mu+ci
    return ci

def aftershock(x0,y0,mag,time,lmax,xx,yy,bg,params):
    x=x0
    y=y0
  
    mu=params["mu"]
    A=params["A"]
    alpha1=params["alpha1"]
    alpha2=params["alpha2"]
    CC=params["c"]
    P=params["P"] 
    coss=params["coss"]
    D =params["D"]
    Q =params["Q"]
    gamma = params["gamma"]
    
    N_mag=len(mag)
    for t0 in range(lmax,N_mag):
        ci1=A*np.exp(alpha1*mag[t0-coss+1:t0])*(1+(time[t0]-time[t0-coss+1:t0])/CC)**(-P)
        ci2=A*np.exp(alpha2*mag[t0-lmax:t0-coss])*(1+(time[t0]-time[t0-lmax:t0-coss])/CC)**(-P)
        
        ci=np.r_[ci2,ci1,mu]
                
        #index=np.random.choice(np.r_[0:len(ci)], size=1, replace=True, p=ci)
        
        ci_df = pd.DataFrame(ci, columns=['w'])
        # Check: df['a'].sample(4, replace=True, weights=df['b'])
        res = ci_df.sample(1, replace=True, weights=ci_df['w'])
        tmp = res.iloc[0]
        index = tmp.name
        if (index<lmax):
            sigma=D*np.exp(gamma*mag[t0])
            index2=t0-lmax+index-1
            
            
            R=sigma*np.sqrt(np.random.uniform()**(1/(1-Q))-1)
            theta=np.random.uniform(0,2*np.pi,(1,))
            rx=R*np.cos(theta)+x[index2]
            ry=R*np.sin(theta)+y[index2]          
        else:
                                  
            res =ci_df.sample(1, replace=True, weights=bg)
            a = res.iloc[0]
            index3 = a.name
                
            rx=xx[index3]
            ry=yy[index3]
              
        x=np.r_[x,rx]
        y=np.r_[y,ry]
      
    return x,y
# X,Y = aftershock(list(df2['x']), list(df2['y']), list(sim_catalog['mag']),
#                list(sim_catalog['time']),
#                l_max, list( italy_bg_density['lon']), list(italy_bg_density['lat']), 
#                list(italy_bg_density['prb']),
#                params)

#%%
# Global catalog

# df = pd.read_csv('database.csv')   # database.csv is the global catalog (m>5)

#df.describe()
#df.head()
'''
times = df["Time"]
dates = df["Date"]
mags = df["Magnitude"].to_numpy()

#fullDates = datetime.combine()

# Create d0, the first catalog date:
format1 = '%m/%d/%Y %H:%M:%S'
format2 = '%m/%d/%Y %H:%M:%S.%f'

d = str(dates[0])
t = str(times[0])

dateStr = d + ' ' + t

d0 = datetime.strptime(dateStr,format1)
print("The first date of the current catalog is " + dateStr)
# creare dates list with days and hours together and store it in EVENT_DATES
EVENT_DATES =[]

for i in range(len(dates)):
    print(i)
    d = dates[i]
    t = times[i]
    dateStr = d + ' ' + t
    print(dateStr)
    try:
        dateTimeObj = datetime.strptime(dateStr, format1)  - d0
        dateTimeObj=dateTimeObj.total_seconds()/86400
    except Exception:  
        dateTimeObj = datetime.strptime(dateStr, format2) - d0
        dateTimeObj=dateTimeObj.total_seconds()/86400

       
    print(dateTimeObj)
    EVENT_DATES.append(dateTimeObj)

#%% Concatenate dates and mags into a df and run ETAS2 to build a catalog
data = np.array( [EVENT_DATES, mags]).T

df2=pd.DataFrame(data, columns=('time', 'mag'))
  
# i<-1
# sim_catlog<-ETAS_T(data=history_data, params=params, sequence.length=st, seed=i, bvalue=1, M0=M0)

sim_catalog = ETAS_T(df2,params, max(dates), seed, bvalue, M0)

''' ;
#%% Italy catalog:
    
df = pd.read_csv('italy_catalog.csv')

t = df['time']
mag = df['mag']
x0 = df['x']
y0 = df['y']

data = np.array( [t, mag, x0, y0]).T

df2=pd.DataFrame(data, columns=('time', 'mag', 'x', 'y'))
l_max=4000  # Number of events that we use from the original catalog to generate ETAS
df2 = df2[-l_max:] 
df2.reset_index(inplace = True, drop = True)
 
#%%
#  ETAS_T(data, params, n, seed, bvalue, M0):  
n_newDays = 365*28
d_final = max(df2.time) + n_newDays    

sim_t,sim_m = ETAS_T(df2,params, d_final, bvalue, params['M0'])
# Aftershocks

italy_bg_density = pd.read_csv('italy_bg_density.csv')
prb = italy_bg_density
# data.xy<-aftershock.xy(history_data$x,history_data$y,sim_catlog$magnitudes-M0,
#                           sim_catlog$times,lmax,
#                           italy_bg_density$lon,
#                           italy_bg_density$lat,
#                           italy_bg_density$prb,
#                           params, seed=i)

t_m = np.array([sim_t,sim_m]).T    
sim_catalog = pd.DataFrame(t_m, columns=('time', 'mag'))    

X,Y = aftershock(df2['x'], df2['y'], sim_catalog['mag'] - params['M0'],
                 sim_catalog['time'],
                 l_max, italy_bg_density['lon'], italy_bg_density['lat'], 
                 italy_bg_density['prb'],
                 params)

data = np.array([sim_t,sim_m,X,Y]).T
sim_catalog = pd.DataFrame(data, columns=('time', 'mag', 'X', 'Y'))    
sim_catalog = sim_catalog[l_max:]


#  def aftershock(x0,y0,mag,time,lmax,xx,yy,bg,params):
#%% plotting:

plt.subplot(2, 1, 1)

#rndIdx = np.random.randint(0, n-n_sim)
rndIdx = 0
#df2.iloc[rndIdx:rndIdx+n_sim, 1].T.plot()

plt.plot(df2['time'], df2['mag'])
#plt.plot(df2['mag'])

plt.title('Original')
plt.ylim(3,7)
plt.subplot(2, 1, 2)

plt.plot(list( sim_catalog['time'])[-l_max:] , list(sim_catalog['mag'])[-l_max:] )

#plt.plot( list(sim_catalog['mag'])[:l_max] )

plt.title('Simulated ETAS2')
plt.ylim(3,7)

plt.show()


#%% HIstogram:

plt.figure()

plt.hist(sim_catalog['mag'], alpha=0.5, label='Sim.', density=True, )
plt.hist(df2['mag'], alpha=0.5, label='Orig.', density=True)
plt.yscale('log')
plt.legend()


plt.figure()

plt.scatter(X,Y)
plt.scatter(df2['x'], df2['y'])  # Original locations

# We want:  time, mag, x, y


data = np.array([sim_catalog['time'],sim_catalog['mag'], 
                        sim_catalog['X'],sim_catalog['Y']]).T
ETAS_cat = pd.DataFrame(data, columns=('time', 'mag', 'X', 'Y'))

ETAS_cat.to_csv('ETAS_sim_cat.csv')












