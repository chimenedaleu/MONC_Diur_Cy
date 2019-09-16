from scipy import interpolate
import numpy as N 
import matplotlib.pyplot as plt
import os 
import netCDF4

import numpy
from scipy.stats import pearsonr
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt


L_vap                = 2.501e6     #J/kg
C_p                  = 1005.0      #J/kgK
R                    = 287.058 
kappa                = R/C_p    
gra                  = 9.8
L_sub=2.85e6 #2.834e6
nz=99
ny=512
nx=512
dx=200
dy=200
length_x=N.zeros(nx)
length_y=N.zeros(ny)

length_x[0]=00.
length_y[0]=00.
for j in N.arange(ny-1):
     length_x[j+1]=length_x[j]+dx
     length_y[j+1]=length_y[j]+dy

length_x=length_x/1000.
length_y=length_y/1000.

hein=N.zeros((nz))
hein = np.load('hein.npy')

precip_xy_ctrl = np.load('surf_precip_xy_time_CTRL.npy')
time_10min_ctrl = np.load('modets_time_10min_CTRL.npy')

precip_xy_ctrl_d2=N.zeros((n_30min,nx,ny))
precip_xy_ctrl_d3=N.zeros((n_30min,nx,ny))
precip_xy_ctrl_d4=N.zeros((n_30min,nx,ny))
precip_xy_ctrl_d5=N.zeros((n_30min,nx,ny))
precip_xy_ctrl_d6=N.zeros((n_30min,nx,ny))
precip_xy_ctrl_d7=N.zeros((n_30min,nx,ny))
precip_xy_ctrl_d8=N.zeros((n_30min,nx,ny))
precip_xy_ctrl_d9=N.zeros((n_30min,nx,ny))
precip_xy_ctrl_d10=N.zeros((n_30min,nx,ny))

hour=N.zeros((n_30min))

time_ctrl_d2[0:92]=time_10min_ctrl[90:182]-24.0
time_ctrl_d3[0:92]=time_10min_ctrl[182:274]-48.0
time_ctrl_d4[0:92]=time_10min_ctrl[273:365]-72.0
time_ctrl_d5[0:92]=time_10min_ctrl[363:455]-96.0
time_ctrl_d6[0:92]=time_10min_ctrl[455:547]-120.0
time_ctrl_d7[0:92]=time_10min_ctrl[547:639]-144.0
time_ctrl_d8[0:92]=time_10min_ctrl[637:729]-168.0
time_ctrl_d9[0:92]=time_10min_ctrl[728:820]-192.0
time_ctrl_d10[0:92]=time_10min_ctrl[819:911]-216.0



precip_xy_ctrl_d2[0:92,:,:]=3600.*precip_xy_ctrl[90:182,:,:]
precip_xy_ctrl_d3[0:92,:,:]=3600.*precip_xy_ctrl[182:274,:,:]
precip_xy_ctrl_d4[0:92,:,:]=3600.*precip_xy_ctrl[273:365,:,:]
precip_xy_ctrl_d5[0:92,:,:]=3600.*precip_xy_ctrl[363:455,:,:]
precip_xy_ctrl_d6[0:92,:,:]=3600.*precip_xy_ctrl[455:547,:,:]
precip_xy_ctrl_d7[0:92,:,:]=3600.*precip_xy_ctrl[547:639,:,:]
precip_xy_ctrl_d8[0:92,:,:]=3600.*precip_xy_ctrl[637:729,:,:]
precip_xy_ctrl_d9[0:92,:,:]=3600.*precip_xy_ctrl[728:820,:,:]
precip_xy_ctrl_d10[0:92,:,:]=3600.*precip_xy_ctrl[819:911,:,:]


n_30min=93
ntr=1
MRR_thr=N.zeros((ntr))
Surf_precip_threshold_Day_mean=N.zeros((ntr))
for i in N.arange(ntr):
     Surf_precip_threshold_Day_mean[i]=0.1
# 50% of the domian--mean daily-mean precipitation rate in the control simulation



nx4=26
ny4=26
length_x4=N.zeros(nx4)
length_y4=N.zeros(ny4)
dx4=4000
dy4=4000


precip_44_ctrl_d2=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d3=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d4=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d5=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d6=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d7=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d8=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d9=N.zeros((n_30min,nx4,ny4))
precip_44_ctrl_d10=N.zeros((n_30min,nx4,ny4))


precip_44_ctrl_d2_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d3_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d4_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d5_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d6_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d7_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d8_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d9_Day_mean=N.zeros((nx4,ny4))
precip_44_ctrl_d10_Day_mean=N.zeros((nx4,ny4))

surf_precip_44_contr_hist_Day_mean=N.zeros((81))
precip_44_ctrl_d6_hist_Day_mean=N.zeros((81))
precip_44_ctrl_d4_hist_Day_mean=N.zeros((81))
precip_44_ctrl_d5_hist_Day_mean=N.zeros((81))
precip_44_ctrl_d3_hist_Day_mean=N.zeros((81))

precip_44_ctrl_d2_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR_lag=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR_lag=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag1=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag2=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag3=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag4=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag5=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag6=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag7=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag8=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag9=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag10=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag11=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag12=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag13=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag14=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag15=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag16=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag17=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag18=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag19=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag20=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag21=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag22=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag23=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag24=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag25=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag26=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag27=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag28=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag29=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag30=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag31=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_lag32=N.zeros((ntr, n_30min))


precip_44_ctrl_d2_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag1=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag2=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag3=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag4=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag6=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag5=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag7=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag8=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag9=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag10=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag11=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag12=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag13=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag14=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag15=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag16=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag17=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag18=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag19=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag20=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag21=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag22=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag23=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag24=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag25=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag26=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag27=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag28=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag29=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag30=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag31=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PR2_lag32=N.zeros((ntr, n_30min))



precip_44_ctrl_d2_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag1=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag2=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag3=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag4=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag5=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag6=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag7=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag8=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag9=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag10=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag11=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag12=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag13=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag14=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag15=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag16=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag17=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag18=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag19=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag20=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag21=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag22=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag23=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag24=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag25=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag26=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag27=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag28=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag29=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag30=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag31=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag32=N.zeros((ntr, n_30min))




precip_44_ctrl_d2_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag33=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag34=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag35=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag36=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag37=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag38=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag39=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag40=N.zeros((ntr, n_30min))



precip_44_ctrl_d2_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag41=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag42=N.zeros((ntr, n_30min))




precip_44_ctrl_d2_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag43=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag44=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag45=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag46=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag47=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag48=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag49=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag50=N.zeros((ntr, n_30min))






precip_44_ctrl_d2_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag1_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag2_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag3_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag4_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag5_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag6_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag7_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag8_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag9_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag10_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag11_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag12_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag13_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag14_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag15_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag16_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag17_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag18_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag19_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag20_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag21_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag22_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag23_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag24_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag25_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag26_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag27_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag28_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag29_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag30_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag31_PR=N.zeros((ntr, n_30min))

precip_44_ctrl_d2_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d3_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d4_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d5_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d6_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d7_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d8_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d9_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))
precip_44_ctrl_d10_inst_PRR_PR2_lag32_PR=N.zeros((ntr, n_30min))

#coarse-grained averaging A=4*4 kilometre square
for j in N.arange(60):
     for ix in N.arange(nx4):
          for iy in N.arange(ny4):
               precip_44_ctrl_d2[j,ix,iy]=precip_xy_ctrl_d2[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d3[j,ix,iy]=precip_xy_ctrl_d3[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d4[j,ix,iy]=precip_xy_ctrl_d4[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d5[j,ix,iy]=precip_xy_ctrl_d5[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d6[j,ix,iy]=precip_xy_ctrl_d6[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d7[j,ix,iy]=precip_xy_ctrl_d7[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d8[j,ix,iy]=precip_xy_ctrl_d8[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d9[j,ix,iy]=precip_xy_ctrl_d9[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               precip_44_ctrl_d10[j,ix,iy]=precip_xy_ctrl_d10[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+19].mean()
               if (ix == 25) :
                    precip_44_ctrl_d2[j,ix,iy]=precip_xy_ctrl_d2[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d3[j,ix,iy]=precip_xy_ctrl_d3[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d4[j,ix,iy]=precip_xy_ctrl_d4[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d5[j,ix,iy]=precip_xy_ctrl_d5[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d6[j,ix,iy]=precip_xy_ctrl_d6[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d7[j,ix,iy]=precip_xy_ctrl_d7[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d8[j,ix,iy]=precip_xy_ctrl_d8[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d9[j,ix,iy]=precip_xy_ctrl_d9[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
                    precip_44_ctrl_d10[j,ix,iy]=precip_xy_ctrl_d10[j,(ix)*20:(ix)*20+11,(iy)*20:(iy)*20+19].mean()
               if (iy == 25) :
                    precip_44_ctrl_d2[j,ix,iy]=precip_xy_ctrl_d2[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d3[j,ix,iy]=precip_xy_ctrl_d3[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d4[j,ix,iy]=precip_xy_ctrl_d4[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d5[j,ix,iy]=precip_xy_ctrl_d5[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d6[j,ix,iy]=precip_xy_ctrl_d6[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d7[j,ix,iy]=precip_xy_ctrl_d7[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d8[j,ix,iy]=precip_xy_ctrl_d8[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d9[j,ix,iy]=precip_xy_ctrl_d9[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()
                    precip_44_ctrl_d10[j,ix,iy]=precip_xy_ctrl_d10[j,(ix)*20:(ix)*20+19,(iy)*20:(iy)*20+11].mean()



#probability of rain at time t given rain a previous time t-dt

for ix in N.arange(nx4):
     for iy in N.arange(ny4):
          for  j in N.arange(n_30min-1):
               for itr in N.arange(ntr):
                    if ( precip_44_ctrl_d2[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d2_inst_PR_lag[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d3[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d3_inst_PR_lag[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d4[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d4_inst_PR_lag[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d5[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d5_inst_PR_lag[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d6[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d6_inst_PR_lag[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d7[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d7_inst_PR_lag[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d8[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d8_inst_PR_lag[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d9[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d9_inst_PR_lag[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j]+1.
                    if ( precip_44_ctrl_d10[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                         precip_44_ctrl_d10_inst_PR_lag[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j]+1.


                    if ( j > 0 ):
                         if ( precip_44_ctrl_d2[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d2[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag1[itr,j]=precip_44_ctrl_d2_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d2[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag2[itr,j]=precip_44_ctrl_d2_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d2[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag3[itr,j]=precip_44_ctrl_d2_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d2[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag4[itr,j]=precip_44_ctrl_d2_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d2[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag5[itr,j]=precip_44_ctrl_d2_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d2[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag6[itr,j]=precip_44_ctrl_d2_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d2[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag7[itr,j]=precip_44_ctrl_d2_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d2[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag8[itr,j]=precip_44_ctrl_d2_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d2[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag9[itr,j]=precip_44_ctrl_d2_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d2[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag10[itr,j]=precip_44_ctrl_d2_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d2[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag11[itr,j]=precip_44_ctrl_d2_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d2[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag12[itr,j]=precip_44_ctrl_d2_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d2[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag13[itr,j]=precip_44_ctrl_d2_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d2[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag14[itr,j]=precip_44_ctrl_d2_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d2[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag15[itr,j]=precip_44_ctrl_d2_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d2[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag16[itr,j]=precip_44_ctrl_d2_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d2[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag17[itr,j]=precip_44_ctrl_d2_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d2[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag18[itr,j]=precip_44_ctrl_d2_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d2[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag19[itr,j]=precip_44_ctrl_d2_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d2[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag20[itr,j]=precip_44_ctrl_d2_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d2[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag21[itr,j]=precip_44_ctrl_d2_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d2[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag22[itr,j]=precip_44_ctrl_d2_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d2[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag23[itr,j]=precip_44_ctrl_d2_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d2[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag24[itr,j]=precip_44_ctrl_d2_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d2[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag25[itr,j]=precip_44_ctrl_d2_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d2[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag26[itr,j]=precip_44_ctrl_d2_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d2[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag27[itr,j]=precip_44_ctrl_d2_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d2[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag28[itr,j]=precip_44_ctrl_d2_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d2[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag29[itr,j]=precip_44_ctrl_d2_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d2[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag30[itr,j]=precip_44_ctrl_d2_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d2[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag31[itr,j]=precip_44_ctrl_d2_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d2[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d2_inst_PRR_lag32[itr,j]=precip_44_ctrl_d2_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d3[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d3[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag1[itr,j]=precip_44_ctrl_d3_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d3[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag2[itr,j]=precip_44_ctrl_d3_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d3[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag3[itr,j]=precip_44_ctrl_d3_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d3[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag4[itr,j]=precip_44_ctrl_d3_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d3[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag5[itr,j]=precip_44_ctrl_d3_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d3[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag6[itr,j]=precip_44_ctrl_d3_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d3[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag7[itr,j]=precip_44_ctrl_d3_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d3[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag8[itr,j]=precip_44_ctrl_d3_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d3[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag9[itr,j]=precip_44_ctrl_d3_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d3[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag10[itr,j]=precip_44_ctrl_d3_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d3[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag11[itr,j]=precip_44_ctrl_d3_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d3[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag12[itr,j]=precip_44_ctrl_d3_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d3[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag13[itr,j]=precip_44_ctrl_d3_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d3[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag14[itr,j]=precip_44_ctrl_d3_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d3[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag15[itr,j]=precip_44_ctrl_d3_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d3[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag16[itr,j]=precip_44_ctrl_d3_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d3[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag17[itr,j]=precip_44_ctrl_d3_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d3[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag18[itr,j]=precip_44_ctrl_d3_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d3[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag19[itr,j]=precip_44_ctrl_d3_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d3[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag20[itr,j]=precip_44_ctrl_d3_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d3[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag21[itr,j]=precip_44_ctrl_d3_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d3[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag22[itr,j]=precip_44_ctrl_d3_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d3[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag23[itr,j]=precip_44_ctrl_d3_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d3[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag24[itr,j]=precip_44_ctrl_d3_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d3[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag25[itr,j]=precip_44_ctrl_d3_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d3[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag26[itr,j]=precip_44_ctrl_d3_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d3[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag27[itr,j]=precip_44_ctrl_d3_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d3[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag28[itr,j]=precip_44_ctrl_d3_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d3[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag29[itr,j]=precip_44_ctrl_d3_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d3[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag30[itr,j]=precip_44_ctrl_d3_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d3[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag31[itr,j]=precip_44_ctrl_d3_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d3[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d3_inst_PRR_lag32[itr,j]=precip_44_ctrl_d3_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d4[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d4[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag1[itr,j]=precip_44_ctrl_d4_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d4[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag2[itr,j]=precip_44_ctrl_d4_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d4[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag3[itr,j]=precip_44_ctrl_d4_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d4[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag4[itr,j]=precip_44_ctrl_d4_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d4[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag5[itr,j]=precip_44_ctrl_d4_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d4[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag6[itr,j]=precip_44_ctrl_d4_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d4[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag7[itr,j]=precip_44_ctrl_d4_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d4[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag8[itr,j]=precip_44_ctrl_d4_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d4[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag9[itr,j]=precip_44_ctrl_d4_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d4[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag10[itr,j]=precip_44_ctrl_d4_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d4[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag11[itr,j]=precip_44_ctrl_d4_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d4[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag12[itr,j]=precip_44_ctrl_d4_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d4[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag13[itr,j]=precip_44_ctrl_d4_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d4[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag14[itr,j]=precip_44_ctrl_d4_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d4[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag15[itr,j]=precip_44_ctrl_d4_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d4[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag16[itr,j]=precip_44_ctrl_d4_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d4[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag17[itr,j]=precip_44_ctrl_d4_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d4[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag18[itr,j]=precip_44_ctrl_d4_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d4[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag19[itr,j]=precip_44_ctrl_d4_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d4[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag20[itr,j]=precip_44_ctrl_d4_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d4[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag21[itr,j]=precip_44_ctrl_d4_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d4[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag22[itr,j]=precip_44_ctrl_d4_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d4[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag23[itr,j]=precip_44_ctrl_d4_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d4[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag24[itr,j]=precip_44_ctrl_d4_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d4[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag25[itr,j]=precip_44_ctrl_d4_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d4[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag26[itr,j]=precip_44_ctrl_d4_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d4[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag27[itr,j]=precip_44_ctrl_d4_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d4[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag28[itr,j]=precip_44_ctrl_d4_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d4[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag29[itr,j]=precip_44_ctrl_d4_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d4[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag30[itr,j]=precip_44_ctrl_d4_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d4[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag31[itr,j]=precip_44_ctrl_d4_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d4[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d4_inst_PRR_lag32[itr,j]=precip_44_ctrl_d4_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d5[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d5[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag1[itr,j]=precip_44_ctrl_d5_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d5[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag2[itr,j]=precip_44_ctrl_d5_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d5[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag3[itr,j]=precip_44_ctrl_d5_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d5[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag4[itr,j]=precip_44_ctrl_d5_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d5[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag5[itr,j]=precip_44_ctrl_d5_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d5[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag6[itr,j]=precip_44_ctrl_d5_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d5[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag7[itr,j]=precip_44_ctrl_d5_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d5[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag8[itr,j]=precip_44_ctrl_d5_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d5[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag9[itr,j]=precip_44_ctrl_d5_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d5[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag10[itr,j]=precip_44_ctrl_d5_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d5[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag11[itr,j]=precip_44_ctrl_d5_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d5[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag12[itr,j]=precip_44_ctrl_d5_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d5[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag13[itr,j]=precip_44_ctrl_d5_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d5[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag14[itr,j]=precip_44_ctrl_d5_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d5[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag15[itr,j]=precip_44_ctrl_d5_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d5[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag16[itr,j]=precip_44_ctrl_d5_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d5[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag17[itr,j]=precip_44_ctrl_d5_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d5[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag18[itr,j]=precip_44_ctrl_d5_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d5[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag19[itr,j]=precip_44_ctrl_d5_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d5[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag20[itr,j]=precip_44_ctrl_d5_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d5[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag21[itr,j]=precip_44_ctrl_d5_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d5[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag22[itr,j]=precip_44_ctrl_d5_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d5[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag23[itr,j]=precip_44_ctrl_d5_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d5[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag24[itr,j]=precip_44_ctrl_d5_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d5[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag25[itr,j]=precip_44_ctrl_d5_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d5[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag26[itr,j]=precip_44_ctrl_d5_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d5[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag27[itr,j]=precip_44_ctrl_d5_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d5[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag28[itr,j]=precip_44_ctrl_d5_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d5[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag29[itr,j]=precip_44_ctrl_d5_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d5[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag30[itr,j]=precip_44_ctrl_d5_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d5[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag31[itr,j]=precip_44_ctrl_d5_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d5[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d5_inst_PRR_lag32[itr,j]=precip_44_ctrl_d5_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d6[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d6[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag1[itr,j]=precip_44_ctrl_d6_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d6[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag2[itr,j]=precip_44_ctrl_d6_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d6[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag3[itr,j]=precip_44_ctrl_d6_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d6[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag4[itr,j]=precip_44_ctrl_d6_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d6[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag5[itr,j]=precip_44_ctrl_d6_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d6[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag6[itr,j]=precip_44_ctrl_d6_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d6[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag7[itr,j]=precip_44_ctrl_d6_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d6[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag8[itr,j]=precip_44_ctrl_d6_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d6[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag9[itr,j]=precip_44_ctrl_d6_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d6[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag10[itr,j]=precip_44_ctrl_d6_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d6[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag11[itr,j]=precip_44_ctrl_d6_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d6[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag12[itr,j]=precip_44_ctrl_d6_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d6[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag13[itr,j]=precip_44_ctrl_d6_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d6[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag14[itr,j]=precip_44_ctrl_d6_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d6[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag15[itr,j]=precip_44_ctrl_d6_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d6[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag16[itr,j]=precip_44_ctrl_d6_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d6[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag17[itr,j]=precip_44_ctrl_d6_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d6[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag18[itr,j]=precip_44_ctrl_d6_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d6[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag19[itr,j]=precip_44_ctrl_d6_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d6[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag20[itr,j]=precip_44_ctrl_d6_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d6[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag21[itr,j]=precip_44_ctrl_d6_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d6[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag22[itr,j]=precip_44_ctrl_d6_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d6[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag23[itr,j]=precip_44_ctrl_d6_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d6[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag24[itr,j]=precip_44_ctrl_d6_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d6[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag25[itr,j]=precip_44_ctrl_d6_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d6[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag26[itr,j]=precip_44_ctrl_d6_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d6[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag27[itr,j]=precip_44_ctrl_d6_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d6[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag28[itr,j]=precip_44_ctrl_d6_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d6[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag29[itr,j]=precip_44_ctrl_d6_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d6[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag30[itr,j]=precip_44_ctrl_d6_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d6[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag31[itr,j]=precip_44_ctrl_d6_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d6[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d6_inst_PRR_lag32[itr,j]=precip_44_ctrl_d6_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d7[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d7[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag1[itr,j]=precip_44_ctrl_d7_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d7[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag2[itr,j]=precip_44_ctrl_d7_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d7[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag3[itr,j]=precip_44_ctrl_d7_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d7[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag4[itr,j]=precip_44_ctrl_d7_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d7[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag5[itr,j]=precip_44_ctrl_d7_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d7[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag6[itr,j]=precip_44_ctrl_d7_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d7[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag7[itr,j]=precip_44_ctrl_d7_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d7[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag8[itr,j]=precip_44_ctrl_d7_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d7[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag9[itr,j]=precip_44_ctrl_d7_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d7[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag10[itr,j]=precip_44_ctrl_d7_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d7[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag11[itr,j]=precip_44_ctrl_d7_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d7[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag12[itr,j]=precip_44_ctrl_d7_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d7[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag13[itr,j]=precip_44_ctrl_d7_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d7[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag14[itr,j]=precip_44_ctrl_d7_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d7[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag15[itr,j]=precip_44_ctrl_d7_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d7[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag16[itr,j]=precip_44_ctrl_d7_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d7[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag17[itr,j]=precip_44_ctrl_d7_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d7[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag18[itr,j]=precip_44_ctrl_d7_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d7[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag19[itr,j]=precip_44_ctrl_d7_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d7[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag20[itr,j]=precip_44_ctrl_d7_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d7[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag21[itr,j]=precip_44_ctrl_d7_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d7[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag22[itr,j]=precip_44_ctrl_d7_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d7[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag23[itr,j]=precip_44_ctrl_d7_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d7[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag24[itr,j]=precip_44_ctrl_d7_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d7[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag25[itr,j]=precip_44_ctrl_d7_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d7[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag26[itr,j]=precip_44_ctrl_d7_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d7[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag27[itr,j]=precip_44_ctrl_d7_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d7[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag28[itr,j]=precip_44_ctrl_d7_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d7[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag29[itr,j]=precip_44_ctrl_d7_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d7[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag30[itr,j]=precip_44_ctrl_d7_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d7[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag31[itr,j]=precip_44_ctrl_d7_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d7[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d7_inst_PRR_lag32[itr,j]=precip_44_ctrl_d7_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d8[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d8[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag1[itr,j]=precip_44_ctrl_d8_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d8[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag2[itr,j]=precip_44_ctrl_d8_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d8[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag3[itr,j]=precip_44_ctrl_d8_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d8[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag4[itr,j]=precip_44_ctrl_d8_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d8[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag5[itr,j]=precip_44_ctrl_d8_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d8[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag6[itr,j]=precip_44_ctrl_d8_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d8[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag7[itr,j]=precip_44_ctrl_d8_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d8[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag8[itr,j]=precip_44_ctrl_d8_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d8[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag9[itr,j]=precip_44_ctrl_d8_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d8[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag10[itr,j]=precip_44_ctrl_d8_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d8[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag11[itr,j]=precip_44_ctrl_d8_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d8[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag12[itr,j]=precip_44_ctrl_d8_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d8[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag13[itr,j]=precip_44_ctrl_d8_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d8[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag14[itr,j]=precip_44_ctrl_d8_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d8[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag15[itr,j]=precip_44_ctrl_d8_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d8[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag16[itr,j]=precip_44_ctrl_d8_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d8[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag17[itr,j]=precip_44_ctrl_d8_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d8[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag18[itr,j]=precip_44_ctrl_d8_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d8[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag19[itr,j]=precip_44_ctrl_d8_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d8[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag20[itr,j]=precip_44_ctrl_d8_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d8[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag21[itr,j]=precip_44_ctrl_d8_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d8[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag22[itr,j]=precip_44_ctrl_d8_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d8[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag23[itr,j]=precip_44_ctrl_d8_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d8[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag24[itr,j]=precip_44_ctrl_d8_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d8[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag25[itr,j]=precip_44_ctrl_d8_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d8[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag26[itr,j]=precip_44_ctrl_d8_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d8[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag27[itr,j]=precip_44_ctrl_d8_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d8[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag28[itr,j]=precip_44_ctrl_d8_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d8[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag29[itr,j]=precip_44_ctrl_d8_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d8[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag30[itr,j]=precip_44_ctrl_d8_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d8[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag31[itr,j]=precip_44_ctrl_d8_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d8[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d8_inst_PRR_lag32[itr,j]=precip_44_ctrl_d8_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d9[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d9[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag1[itr,j]=precip_44_ctrl_d9_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d9[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag2[itr,j]=precip_44_ctrl_d9_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d9[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag3[itr,j]=precip_44_ctrl_d9_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d9[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag4[itr,j]=precip_44_ctrl_d9_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d9[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag5[itr,j]=precip_44_ctrl_d9_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d9[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag6[itr,j]=precip_44_ctrl_d9_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d9[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag7[itr,j]=precip_44_ctrl_d9_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d9[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag8[itr,j]=precip_44_ctrl_d9_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d9[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag9[itr,j]=precip_44_ctrl_d9_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d9[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag10[itr,j]=precip_44_ctrl_d9_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d9[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag11[itr,j]=precip_44_ctrl_d9_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d9[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag12[itr,j]=precip_44_ctrl_d9_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d9[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag13[itr,j]=precip_44_ctrl_d9_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d9[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag14[itr,j]=precip_44_ctrl_d9_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d9[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag15[itr,j]=precip_44_ctrl_d9_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d9[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag16[itr,j]=precip_44_ctrl_d9_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d9[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag17[itr,j]=precip_44_ctrl_d9_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d9[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag18[itr,j]=precip_44_ctrl_d9_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d9[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag19[itr,j]=precip_44_ctrl_d9_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d9[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag20[itr,j]=precip_44_ctrl_d9_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d9[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag21[itr,j]=precip_44_ctrl_d9_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d9[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag22[itr,j]=precip_44_ctrl_d9_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d9[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag23[itr,j]=precip_44_ctrl_d9_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d9[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag24[itr,j]=precip_44_ctrl_d9_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d9[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag25[itr,j]=precip_44_ctrl_d9_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d9[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag26[itr,j]=precip_44_ctrl_d9_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d9[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag27[itr,j]=precip_44_ctrl_d9_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d9[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag28[itr,j]=precip_44_ctrl_d9_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d9[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag29[itr,j]=precip_44_ctrl_d9_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d9[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag30[itr,j]=precip_44_ctrl_d9_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d9[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag31[itr,j]=precip_44_ctrl_d9_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d9[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d9_inst_PRR_lag32[itr,j]=precip_44_ctrl_d9_inst_PRR_lag32[itr,j]+1.

                         if ( precip_44_ctrl_d10[j,ix,iy] >= Surf_precip_threshold_Day_mean[itr] ) :
                              if ( j >=1 and precip_44_ctrl_d10[j-1,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag1[itr,j]=precip_44_ctrl_d10_inst_PRR_lag1[itr,j]+1.
                              if ( j >=2 and precip_44_ctrl_d10[j-2,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag2[itr,j]=precip_44_ctrl_d10_inst_PRR_lag2[itr,j]+1.
                              if ( j >=3 and precip_44_ctrl_d10[j-3,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag3[itr,j]=precip_44_ctrl_d10_inst_PRR_lag3[itr,j]+1.
                              if ( j >=4 and precip_44_ctrl_d10[j-4,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag4[itr,j]=precip_44_ctrl_d10_inst_PRR_lag4[itr,j]+1.
                              if ( j >=5 and precip_44_ctrl_d10[j-5,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag5[itr,j]=precip_44_ctrl_d10_inst_PRR_lag5[itr,j]+1.
                              if ( j >=6 and precip_44_ctrl_d10[j-6,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag6[itr,j]=precip_44_ctrl_d10_inst_PRR_lag6[itr,j]+1.
                              if ( j >=7 and precip_44_ctrl_d10[j-7,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag7[itr,j]=precip_44_ctrl_d10_inst_PRR_lag7[itr,j]+1.
                              if ( j >=8 and precip_44_ctrl_d10[j-8,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag8[itr,j]=precip_44_ctrl_d10_inst_PRR_lag8[itr,j]+1.
                              if ( j >=9 and precip_44_ctrl_d10[j-9,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag9[itr,j]=precip_44_ctrl_d10_inst_PRR_lag9[itr,j]+1.
                              if ( j >=10 and precip_44_ctrl_d10[j-10,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag10[itr,j]=precip_44_ctrl_d10_inst_PRR_lag10[itr,j]+1.
                              if ( j >=11 and precip_44_ctrl_d10[j-11,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag11[itr,j]=precip_44_ctrl_d10_inst_PRR_lag11[itr,j]+1.
                              if ( j >=12 and precip_44_ctrl_d10[j-12,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag12[itr,j]=precip_44_ctrl_d10_inst_PRR_lag12[itr,j]+1.
                              if ( j >=13 and precip_44_ctrl_d10[j-13,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag13[itr,j]=precip_44_ctrl_d10_inst_PRR_lag13[itr,j]+1.
                              if ( j >=14 and precip_44_ctrl_d10[j-14,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag14[itr,j]=precip_44_ctrl_d10_inst_PRR_lag14[itr,j]+1.
                              if ( j >=15 and precip_44_ctrl_d10[j-15,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag15[itr,j]=precip_44_ctrl_d10_inst_PRR_lag15[itr,j]+1.
                              if ( j >=16 and precip_44_ctrl_d10[j-16,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag16[itr,j]=precip_44_ctrl_d10_inst_PRR_lag16[itr,j]+1.
                              if ( j >=17 and precip_44_ctrl_d10[j-17,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag17[itr,j]=precip_44_ctrl_d10_inst_PRR_lag17[itr,j]+1.
                              if ( j >=18 and precip_44_ctrl_d10[j-18,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag18[itr,j]=precip_44_ctrl_d10_inst_PRR_lag18[itr,j]+1.
                              if ( j >=19 and precip_44_ctrl_d10[j-19,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag19[itr,j]=precip_44_ctrl_d10_inst_PRR_lag19[itr,j]+1.
                              if ( j >=20 and precip_44_ctrl_d10[j-20,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag20[itr,j]=precip_44_ctrl_d10_inst_PRR_lag20[itr,j]+1.
                              if ( j >=21 and precip_44_ctrl_d10[j-21,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag21[itr,j]=precip_44_ctrl_d10_inst_PRR_lag21[itr,j]+1.
                              if ( j >=22 and precip_44_ctrl_d10[j-22,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag22[itr,j]=precip_44_ctrl_d10_inst_PRR_lag22[itr,j]+1.
                              if ( j >=23 and precip_44_ctrl_d10[j-23,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag23[itr,j]=precip_44_ctrl_d10_inst_PRR_lag23[itr,j]+1.
                              if ( j >=24 and precip_44_ctrl_d10[j-24,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag24[itr,j]=precip_44_ctrl_d10_inst_PRR_lag24[itr,j]+1.
                              if ( j >=25 and precip_44_ctrl_d10[j-25,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag25[itr,j]=precip_44_ctrl_d10_inst_PRR_lag25[itr,j]+1.
                              if ( j >=26 and precip_44_ctrl_d10[j-26,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag26[itr,j]=precip_44_ctrl_d10_inst_PRR_lag26[itr,j]+1.
                              if ( j >=27 and precip_44_ctrl_d10[j-27,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag27[itr,j]=precip_44_ctrl_d10_inst_PRR_lag27[itr,j]+1.
                              if ( j >=28 and precip_44_ctrl_d10[j-28,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag28[itr,j]=precip_44_ctrl_d10_inst_PRR_lag28[itr,j]+1.
                              if ( j >=29 and precip_44_ctrl_d10[j-29,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag29[itr,j]=precip_44_ctrl_d10_inst_PRR_lag29[itr,j]+1.
                              if ( j >=30 and precip_44_ctrl_d10[j-30,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag30[itr,j]=precip_44_ctrl_d10_inst_PRR_lag30[itr,j]+1.
                              if ( j >=31 and precip_44_ctrl_d10[j-31,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag31[itr,j]=precip_44_ctrl_d10_inst_PRR_lag31[itr,j]+1.
                              if ( j >=32 and precip_44_ctrl_d10[j-32,ix,iy] >= Surf_precip_threshold_Day_mean[itr]):
                                   precip_44_ctrl_d10_inst_PRR_lag32[itr,j]=precip_44_ctrl_d10_inst_PRR_lag32[itr,j]+1.




precip_44_ctrl_d2_inst_PR_lag[:,:]=precip_44_ctrl_d2_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag1[:,:]=precip_44_ctrl_d2_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag2[:,:]=precip_44_ctrl_d2_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag3[:,:]=precip_44_ctrl_d2_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag4[:,:]=precip_44_ctrl_d2_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag5[:,:]=precip_44_ctrl_d2_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag6[:,:]=precip_44_ctrl_d2_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag7[:,:]=precip_44_ctrl_d2_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag8[:,:]=precip_44_ctrl_d2_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag9[:,:]=precip_44_ctrl_d2_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag10[:,:]=precip_44_ctrl_d2_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag11[:,:]=precip_44_ctrl_d2_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag12[:,:]=precip_44_ctrl_d2_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag13[:,:]=precip_44_ctrl_d2_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag14[:,:]=precip_44_ctrl_d2_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag15[:,:]=precip_44_ctrl_d2_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag16[:,:]=precip_44_ctrl_d2_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag17[:,:]=precip_44_ctrl_d2_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag18[:,:]=precip_44_ctrl_d2_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag19[:,:]=precip_44_ctrl_d2_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag20[:,:]=precip_44_ctrl_d2_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag21[:,:]=precip_44_ctrl_d2_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag22[:,:]=precip_44_ctrl_d2_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag23[:,:]=precip_44_ctrl_d2_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag24[:,:]=precip_44_ctrl_d2_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag25[:,:]=precip_44_ctrl_d2_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag26[:,:]=precip_44_ctrl_d2_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag27[:,:]=precip_44_ctrl_d2_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag28[:,:]=precip_44_ctrl_d2_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag29[:,:]=precip_44_ctrl_d2_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag30[:,:]=precip_44_ctrl_d2_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag31[:,:]=precip_44_ctrl_d2_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d2_inst_PRR_lag32[:,:]=precip_44_ctrl_d2_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d3_inst_PR_lag[:,:]=precip_44_ctrl_d3_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag1[:,:]=precip_44_ctrl_d3_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag2[:,:]=precip_44_ctrl_d3_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag3[:,:]=precip_44_ctrl_d3_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag4[:,:]=precip_44_ctrl_d3_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag5[:,:]=precip_44_ctrl_d3_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag6[:,:]=precip_44_ctrl_d3_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag7[:,:]=precip_44_ctrl_d3_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag8[:,:]=precip_44_ctrl_d3_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag9[:,:]=precip_44_ctrl_d3_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag10[:,:]=precip_44_ctrl_d3_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag11[:,:]=precip_44_ctrl_d3_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag12[:,:]=precip_44_ctrl_d3_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag13[:,:]=precip_44_ctrl_d3_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag14[:,:]=precip_44_ctrl_d3_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag15[:,:]=precip_44_ctrl_d3_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag16[:,:]=precip_44_ctrl_d3_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag17[:,:]=precip_44_ctrl_d3_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag18[:,:]=precip_44_ctrl_d3_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag19[:,:]=precip_44_ctrl_d3_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag20[:,:]=precip_44_ctrl_d3_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag21[:,:]=precip_44_ctrl_d3_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag22[:,:]=precip_44_ctrl_d3_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag23[:,:]=precip_44_ctrl_d3_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag24[:,:]=precip_44_ctrl_d3_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag25[:,:]=precip_44_ctrl_d3_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag26[:,:]=precip_44_ctrl_d3_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag27[:,:]=precip_44_ctrl_d3_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag28[:,:]=precip_44_ctrl_d3_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag29[:,:]=precip_44_ctrl_d3_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag30[:,:]=precip_44_ctrl_d3_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag31[:,:]=precip_44_ctrl_d3_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d3_inst_PRR_lag32[:,:]=precip_44_ctrl_d3_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)


precip_44_ctrl_d4_inst_PR_lag[:,:]=precip_44_ctrl_d4_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag1[:,:]=precip_44_ctrl_d4_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRNR_lag1[:,:]=precip_44_ctrl_d4_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag2[:,:]=precip_44_ctrl_d4_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag3[:,:]=precip_44_ctrl_d4_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag4[:,:]=precip_44_ctrl_d4_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag5[:,:]=precip_44_ctrl_d4_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag6[:,:]=precip_44_ctrl_d4_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag7[:,:]=precip_44_ctrl_d4_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag8[:,:]=precip_44_ctrl_d4_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag9[:,:]=precip_44_ctrl_d4_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag10[:,:]=precip_44_ctrl_d4_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag11[:,:]=precip_44_ctrl_d4_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag12[:,:]=precip_44_ctrl_d4_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag13[:,:]=precip_44_ctrl_d4_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag14[:,:]=precip_44_ctrl_d4_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag15[:,:]=precip_44_ctrl_d4_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag16[:,:]=precip_44_ctrl_d4_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag17[:,:]=precip_44_ctrl_d4_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag18[:,:]=precip_44_ctrl_d4_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag19[:,:]=precip_44_ctrl_d4_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag20[:,:]=precip_44_ctrl_d4_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag21[:,:]=precip_44_ctrl_d4_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag22[:,:]=precip_44_ctrl_d4_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag23[:,:]=precip_44_ctrl_d4_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag24[:,:]=precip_44_ctrl_d4_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag25[:,:]=precip_44_ctrl_d4_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag26[:,:]=precip_44_ctrl_d4_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag27[:,:]=precip_44_ctrl_d4_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag28[:,:]=precip_44_ctrl_d4_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag29[:,:]=precip_44_ctrl_d4_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag30[:,:]=precip_44_ctrl_d4_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag31[:,:]=precip_44_ctrl_d4_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d4_inst_PRR_lag32[:,:]=precip_44_ctrl_d4_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d5_inst_PR_lag[:,:]=precip_44_ctrl_d5_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag1[:,:]=precip_44_ctrl_d5_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRNR_lag1[:,:]=precip_44_ctrl_d5_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag2[:,:]=precip_44_ctrl_d5_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag3[:,:]=precip_44_ctrl_d5_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag4[:,:]=precip_44_ctrl_d5_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag5[:,:]=precip_44_ctrl_d5_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag6[:,:]=precip_44_ctrl_d5_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag7[:,:]=precip_44_ctrl_d5_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag8[:,:]=precip_44_ctrl_d5_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag9[:,:]=precip_44_ctrl_d5_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag10[:,:]=precip_44_ctrl_d5_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag11[:,:]=precip_44_ctrl_d5_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag12[:,:]=precip_44_ctrl_d5_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag13[:,:]=precip_44_ctrl_d5_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag14[:,:]=precip_44_ctrl_d5_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag15[:,:]=precip_44_ctrl_d5_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag16[:,:]=precip_44_ctrl_d5_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag17[:,:]=precip_44_ctrl_d5_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag18[:,:]=precip_44_ctrl_d5_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag19[:,:]=precip_44_ctrl_d5_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag20[:,:]=precip_44_ctrl_d5_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag21[:,:]=precip_44_ctrl_d5_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag22[:,:]=precip_44_ctrl_d5_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag23[:,:]=precip_44_ctrl_d5_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag24[:,:]=precip_44_ctrl_d5_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag25[:,:]=precip_44_ctrl_d5_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag26[:,:]=precip_44_ctrl_d5_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag27[:,:]=precip_44_ctrl_d5_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag28[:,:]=precip_44_ctrl_d5_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag29[:,:]=precip_44_ctrl_d5_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag30[:,:]=precip_44_ctrl_d5_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag31[:,:]=precip_44_ctrl_d5_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d5_inst_PRR_lag32[:,:]=precip_44_ctrl_d5_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d6_inst_PR_lag[:,:]=precip_44_ctrl_d6_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag1[:,:]=precip_44_ctrl_d6_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRNR_lag1[:,:]=precip_44_ctrl_d6_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag2[:,:]=precip_44_ctrl_d6_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag3[:,:]=precip_44_ctrl_d6_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag4[:,:]=precip_44_ctrl_d6_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag5[:,:]=precip_44_ctrl_d6_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag6[:,:]=precip_44_ctrl_d6_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag7[:,:]=precip_44_ctrl_d6_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag8[:,:]=precip_44_ctrl_d6_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag9[:,:]=precip_44_ctrl_d6_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag10[:,:]=precip_44_ctrl_d6_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag11[:,:]=precip_44_ctrl_d6_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag12[:,:]=precip_44_ctrl_d6_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag13[:,:]=precip_44_ctrl_d6_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag14[:,:]=precip_44_ctrl_d6_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag15[:,:]=precip_44_ctrl_d6_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag16[:,:]=precip_44_ctrl_d6_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag17[:,:]=precip_44_ctrl_d6_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag18[:,:]=precip_44_ctrl_d6_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag19[:,:]=precip_44_ctrl_d6_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag20[:,:]=precip_44_ctrl_d6_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag21[:,:]=precip_44_ctrl_d6_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag22[:,:]=precip_44_ctrl_d6_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag23[:,:]=precip_44_ctrl_d6_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag24[:,:]=precip_44_ctrl_d6_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag25[:,:]=precip_44_ctrl_d6_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag26[:,:]=precip_44_ctrl_d6_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag27[:,:]=precip_44_ctrl_d6_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag28[:,:]=precip_44_ctrl_d6_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag29[:,:]=precip_44_ctrl_d6_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag30[:,:]=precip_44_ctrl_d6_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag31[:,:]=precip_44_ctrl_d6_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d6_inst_PRR_lag32[:,:]=precip_44_ctrl_d6_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d7_inst_PR_lag[:,:]=precip_44_ctrl_d7_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag1[:,:]=precip_44_ctrl_d7_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRNR_lag1[:,:]=precip_44_ctrl_d7_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag2[:,:]=precip_44_ctrl_d7_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag3[:,:]=precip_44_ctrl_d7_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag4[:,:]=precip_44_ctrl_d7_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag5[:,:]=precip_44_ctrl_d7_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag6[:,:]=precip_44_ctrl_d7_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag7[:,:]=precip_44_ctrl_d7_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag8[:,:]=precip_44_ctrl_d7_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag9[:,:]=precip_44_ctrl_d7_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag10[:,:]=precip_44_ctrl_d7_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag11[:,:]=precip_44_ctrl_d7_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag12[:,:]=precip_44_ctrl_d7_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag13[:,:]=precip_44_ctrl_d7_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag14[:,:]=precip_44_ctrl_d7_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag15[:,:]=precip_44_ctrl_d7_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag16[:,:]=precip_44_ctrl_d7_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag17[:,:]=precip_44_ctrl_d7_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag18[:,:]=precip_44_ctrl_d7_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag19[:,:]=precip_44_ctrl_d7_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag20[:,:]=precip_44_ctrl_d7_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag21[:,:]=precip_44_ctrl_d7_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag22[:,:]=precip_44_ctrl_d7_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag23[:,:]=precip_44_ctrl_d7_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag24[:,:]=precip_44_ctrl_d7_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag25[:,:]=precip_44_ctrl_d7_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag26[:,:]=precip_44_ctrl_d7_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag27[:,:]=precip_44_ctrl_d7_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag28[:,:]=precip_44_ctrl_d7_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag29[:,:]=precip_44_ctrl_d7_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag30[:,:]=precip_44_ctrl_d7_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag31[:,:]=precip_44_ctrl_d7_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d7_inst_PRR_lag32[:,:]=precip_44_ctrl_d7_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d8_inst_PR_lag[:,:]=precip_44_ctrl_d8_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag1[:,:]=precip_44_ctrl_d8_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRNR_lag1[:,:]=precip_44_ctrl_d8_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag2[:,:]=precip_44_ctrl_d8_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag3[:,:]=precip_44_ctrl_d8_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag4[:,:]=precip_44_ctrl_d8_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag5[:,:]=precip_44_ctrl_d8_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag6[:,:]=precip_44_ctrl_d8_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag7[:,:]=precip_44_ctrl_d8_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag8[:,:]=precip_44_ctrl_d8_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag9[:,:]=precip_44_ctrl_d8_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag10[:,:]=precip_44_ctrl_d8_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag11[:,:]=precip_44_ctrl_d8_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag12[:,:]=precip_44_ctrl_d8_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag13[:,:]=precip_44_ctrl_d8_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag14[:,:]=precip_44_ctrl_d8_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag15[:,:]=precip_44_ctrl_d8_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag16[:,:]=precip_44_ctrl_d8_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag17[:,:]=precip_44_ctrl_d8_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag18[:,:]=precip_44_ctrl_d8_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag19[:,:]=precip_44_ctrl_d8_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag20[:,:]=precip_44_ctrl_d8_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag21[:,:]=precip_44_ctrl_d8_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag22[:,:]=precip_44_ctrl_d8_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag23[:,:]=precip_44_ctrl_d8_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag24[:,:]=precip_44_ctrl_d8_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag25[:,:]=precip_44_ctrl_d8_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag26[:,:]=precip_44_ctrl_d8_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag27[:,:]=precip_44_ctrl_d8_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag28[:,:]=precip_44_ctrl_d8_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag29[:,:]=precip_44_ctrl_d8_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag30[:,:]=precip_44_ctrl_d8_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag31[:,:]=precip_44_ctrl_d8_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d8_inst_PRR_lag32[:,:]=precip_44_ctrl_d8_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d9_inst_PR_lag[:,:]=precip_44_ctrl_d9_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag1[:,:]=precip_44_ctrl_d9_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRNR_lag1[:,:]=precip_44_ctrl_d9_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag2[:,:]=precip_44_ctrl_d9_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag3[:,:]=precip_44_ctrl_d9_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag4[:,:]=precip_44_ctrl_d9_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag5[:,:]=precip_44_ctrl_d9_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag6[:,:]=precip_44_ctrl_d9_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag7[:,:]=precip_44_ctrl_d9_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag8[:,:]=precip_44_ctrl_d9_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag9[:,:]=precip_44_ctrl_d9_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag10[:,:]=precip_44_ctrl_d9_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag11[:,:]=precip_44_ctrl_d9_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag12[:,:]=precip_44_ctrl_d9_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag13[:,:]=precip_44_ctrl_d9_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag14[:,:]=precip_44_ctrl_d9_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag15[:,:]=precip_44_ctrl_d9_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag16[:,:]=precip_44_ctrl_d9_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag17[:,:]=precip_44_ctrl_d9_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag18[:,:]=precip_44_ctrl_d9_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag19[:,:]=precip_44_ctrl_d9_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag20[:,:]=precip_44_ctrl_d9_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag21[:,:]=precip_44_ctrl_d9_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag22[:,:]=precip_44_ctrl_d9_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag23[:,:]=precip_44_ctrl_d9_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag24[:,:]=precip_44_ctrl_d9_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag25[:,:]=precip_44_ctrl_d9_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag26[:,:]=precip_44_ctrl_d9_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag27[:,:]=precip_44_ctrl_d9_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag28[:,:]=precip_44_ctrl_d9_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag29[:,:]=precip_44_ctrl_d9_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag30[:,:]=precip_44_ctrl_d9_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag31[:,:]=precip_44_ctrl_d9_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d9_inst_PRR_lag32[:,:]=precip_44_ctrl_d9_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)

precip_44_ctrl_d10_inst_PR_lag[:,:]=precip_44_ctrl_d10_inst_PR_lag[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag1[:,:]=precip_44_ctrl_d10_inst_PRR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRNR_lag1[:,:]=precip_44_ctrl_d10_inst_PRNR_lag1[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag2[:,:]=precip_44_ctrl_d10_inst_PRR_lag2[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag3[:,:]=precip_44_ctrl_d10_inst_PRR_lag3[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag4[:,:]=precip_44_ctrl_d10_inst_PRR_lag4[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag5[:,:]=precip_44_ctrl_d10_inst_PRR_lag5[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag6[:,:]=precip_44_ctrl_d10_inst_PRR_lag6[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag7[:,:]=precip_44_ctrl_d10_inst_PRR_lag7[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag8[:,:]=precip_44_ctrl_d10_inst_PRR_lag8[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag9[:,:]=precip_44_ctrl_d10_inst_PRR_lag9[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag10[:,:]=precip_44_ctrl_d10_inst_PRR_lag10[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag11[:,:]=precip_44_ctrl_d10_inst_PRR_lag11[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag12[:,:]=precip_44_ctrl_d10_inst_PRR_lag12[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag13[:,:]=precip_44_ctrl_d10_inst_PRR_lag13[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag14[:,:]=precip_44_ctrl_d10_inst_PRR_lag14[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag15[:,:]=precip_44_ctrl_d10_inst_PRR_lag15[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag16[:,:]=precip_44_ctrl_d10_inst_PRR_lag16[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag17[:,:]=precip_44_ctrl_d10_inst_PRR_lag17[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag18[:,:]=precip_44_ctrl_d10_inst_PRR_lag18[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag19[:,:]=precip_44_ctrl_d10_inst_PRR_lag19[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag20[:,:]=precip_44_ctrl_d10_inst_PRR_lag20[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag21[:,:]=precip_44_ctrl_d10_inst_PRR_lag21[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag22[:,:]=precip_44_ctrl_d10_inst_PRR_lag22[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag23[:,:]=precip_44_ctrl_d10_inst_PRR_lag23[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag24[:,:]=precip_44_ctrl_d10_inst_PRR_lag24[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag25[:,:]=precip_44_ctrl_d10_inst_PRR_lag25[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag26[:,:]=precip_44_ctrl_d10_inst_PRR_lag26[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag27[:,:]=precip_44_ctrl_d10_inst_PRR_lag27[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag28[:,:]=precip_44_ctrl_d10_inst_PRR_lag28[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag29[:,:]=precip_44_ctrl_d10_inst_PRR_lag29[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag30[:,:]=precip_44_ctrl_d10_inst_PRR_lag30[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag31[:,:]=precip_44_ctrl_d10_inst_PRR_lag31[:,:]/(nx4*ny4*1.0)
precip_44_ctrl_d10_inst_PRR_lag32[:,:]=precip_44_ctrl_d10_inst_PRR_lag32[:,:]/(nx4*ny4*1.0)


for itr in N.arange(ntr):
     for j in N.arange(n_30min-1):


          if ( j >=1 ) :
               precip_44_ctrl_d2_inst_PR2_lag1[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-1]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag1[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-1]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag1[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-1]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag1[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-1]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag1[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-1]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag1[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-1]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag1[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-1]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag1[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-1]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag1[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-1]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d2_inst_PRR_lag1[itr,j]-precip_44_ctrl_d2_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d3_inst_PRR_lag1[itr,j]-precip_44_ctrl_d3_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d4_inst_PRR_lag1[itr,j]-precip_44_ctrl_d4_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d5_inst_PRR_lag1[itr,j]-precip_44_ctrl_d5_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d6_inst_PRR_lag1[itr,j]-precip_44_ctrl_d6_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d7_inst_PRR_lag1[itr,j]-precip_44_ctrl_d7_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d8_inst_PRR_lag1[itr,j]-precip_44_ctrl_d8_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d9_inst_PRR_lag1[itr,j]-precip_44_ctrl_d9_inst_PR2_lag1[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag1[itr,j]=precip_44_ctrl_d10_inst_PRR_lag1[itr,j]-precip_44_ctrl_d10_inst_PR2_lag1[itr,j]


          if ( j >=2 ) :
               precip_44_ctrl_d2_inst_PR2_lag2[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-2]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag2[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-2]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag2[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-2]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag2[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-2]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag2[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-2]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag2[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-2]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag2[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-2]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag2[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-2]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag2[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-2]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d2_inst_PRR_lag2[itr,j]-precip_44_ctrl_d2_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d3_inst_PRR_lag2[itr,j]-precip_44_ctrl_d3_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d4_inst_PRR_lag2[itr,j]-precip_44_ctrl_d4_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d5_inst_PRR_lag2[itr,j]-precip_44_ctrl_d5_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d6_inst_PRR_lag2[itr,j]-precip_44_ctrl_d6_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d7_inst_PRR_lag2[itr,j]-precip_44_ctrl_d7_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d8_inst_PRR_lag2[itr,j]-precip_44_ctrl_d8_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d9_inst_PRR_lag2[itr,j]-precip_44_ctrl_d9_inst_PR2_lag2[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag2[itr,j]=precip_44_ctrl_d10_inst_PRR_lag2[itr,j]-precip_44_ctrl_d10_inst_PR2_lag2[itr,j]


          if ( j >=3 ) :
               precip_44_ctrl_d2_inst_PR2_lag3[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-3]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag3[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-3]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag3[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-3]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag3[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-3]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag3[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-3]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag3[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-3]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag3[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-3]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag3[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-3]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag3[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-3]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d2_inst_PRR_lag3[itr,j]-precip_44_ctrl_d2_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d3_inst_PRR_lag3[itr,j]-precip_44_ctrl_d3_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d4_inst_PRR_lag3[itr,j]-precip_44_ctrl_d4_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d5_inst_PRR_lag3[itr,j]-precip_44_ctrl_d5_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d6_inst_PRR_lag3[itr,j]-precip_44_ctrl_d6_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d7_inst_PRR_lag3[itr,j]-precip_44_ctrl_d7_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d8_inst_PRR_lag3[itr,j]-precip_44_ctrl_d8_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d9_inst_PRR_lag3[itr,j]-precip_44_ctrl_d9_inst_PR2_lag3[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag3[itr,j]=precip_44_ctrl_d10_inst_PRR_lag3[itr,j]-precip_44_ctrl_d10_inst_PR2_lag3[itr,j]


          if ( j >=4 ) :
               precip_44_ctrl_d2_inst_PR2_lag4[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-4]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag4[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-4]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag4[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-4]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag4[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-4]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag4[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-4]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag4[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-4]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag4[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-4]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag4[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-4]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag4[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-4]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d2_inst_PRR_lag4[itr,j]-precip_44_ctrl_d2_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d3_inst_PRR_lag4[itr,j]-precip_44_ctrl_d3_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d4_inst_PRR_lag4[itr,j]-precip_44_ctrl_d4_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d5_inst_PRR_lag4[itr,j]-precip_44_ctrl_d5_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d6_inst_PRR_lag4[itr,j]-precip_44_ctrl_d6_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d7_inst_PRR_lag4[itr,j]-precip_44_ctrl_d7_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d8_inst_PRR_lag4[itr,j]-precip_44_ctrl_d8_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d9_inst_PRR_lag4[itr,j]-precip_44_ctrl_d9_inst_PR2_lag4[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag4[itr,j]=precip_44_ctrl_d10_inst_PRR_lag4[itr,j]-precip_44_ctrl_d10_inst_PR2_lag4[itr,j]


          if ( j >=5 ) :
               precip_44_ctrl_d2_inst_PR2_lag5[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-5]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag5[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-5]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag5[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-5]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag5[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-5]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag5[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-5]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag5[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-5]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag5[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-5]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag5[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-5]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag5[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-5]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d2_inst_PRR_lag5[itr,j]-precip_44_ctrl_d2_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d3_inst_PRR_lag5[itr,j]-precip_44_ctrl_d3_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d4_inst_PRR_lag5[itr,j]-precip_44_ctrl_d4_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d5_inst_PRR_lag5[itr,j]-precip_44_ctrl_d5_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d6_inst_PRR_lag5[itr,j]-precip_44_ctrl_d6_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d7_inst_PRR_lag5[itr,j]-precip_44_ctrl_d7_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d8_inst_PRR_lag5[itr,j]-precip_44_ctrl_d8_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d9_inst_PRR_lag5[itr,j]-precip_44_ctrl_d9_inst_PR2_lag5[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag5[itr,j]=precip_44_ctrl_d10_inst_PRR_lag5[itr,j]-precip_44_ctrl_d10_inst_PR2_lag5[itr,j]


          if ( j >=6 ) :
               precip_44_ctrl_d2_inst_PR2_lag6[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-6]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag6[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-6]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag6[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-6]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag6[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-6]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag6[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-6]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag6[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-6]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag6[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-6]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag6[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-6]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag6[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-6]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d2_inst_PRR_lag6[itr,j]-precip_44_ctrl_d2_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d3_inst_PRR_lag6[itr,j]-precip_44_ctrl_d3_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d4_inst_PRR_lag6[itr,j]-precip_44_ctrl_d4_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d5_inst_PRR_lag6[itr,j]-precip_44_ctrl_d5_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d6_inst_PRR_lag6[itr,j]-precip_44_ctrl_d6_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d7_inst_PRR_lag6[itr,j]-precip_44_ctrl_d7_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d8_inst_PRR_lag6[itr,j]-precip_44_ctrl_d8_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d9_inst_PRR_lag6[itr,j]-precip_44_ctrl_d9_inst_PR2_lag6[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag6[itr,j]=precip_44_ctrl_d10_inst_PRR_lag6[itr,j]-precip_44_ctrl_d10_inst_PR2_lag6[itr,j]


          if ( j >=7 ) :
               precip_44_ctrl_d2_inst_PR2_lag7[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-7]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag7[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-7]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag7[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-7]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag7[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-7]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag7[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-7]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag7[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-7]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag7[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-7]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag7[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-7]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag7[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-7]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d2_inst_PRR_lag7[itr,j]-precip_44_ctrl_d2_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d3_inst_PRR_lag7[itr,j]-precip_44_ctrl_d3_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d4_inst_PRR_lag7[itr,j]-precip_44_ctrl_d4_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d5_inst_PRR_lag7[itr,j]-precip_44_ctrl_d5_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d6_inst_PRR_lag7[itr,j]-precip_44_ctrl_d6_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d7_inst_PRR_lag7[itr,j]-precip_44_ctrl_d7_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d8_inst_PRR_lag7[itr,j]-precip_44_ctrl_d8_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d9_inst_PRR_lag7[itr,j]-precip_44_ctrl_d9_inst_PR2_lag7[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag7[itr,j]=precip_44_ctrl_d10_inst_PRR_lag7[itr,j]-precip_44_ctrl_d10_inst_PR2_lag7[itr,j]

          if ( j >=8 ) :
               precip_44_ctrl_d2_inst_PR2_lag8[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-8]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag8[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-8]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag8[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-8]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag8[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-8]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag8[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-8]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag8[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-8]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag8[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-8]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag8[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-8]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag8[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-8]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d2_inst_PRR_lag8[itr,j]-precip_44_ctrl_d2_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d3_inst_PRR_lag8[itr,j]-precip_44_ctrl_d3_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d4_inst_PRR_lag8[itr,j]-precip_44_ctrl_d4_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d5_inst_PRR_lag8[itr,j]-precip_44_ctrl_d5_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d6_inst_PRR_lag8[itr,j]-precip_44_ctrl_d6_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d7_inst_PRR_lag8[itr,j]-precip_44_ctrl_d7_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d8_inst_PRR_lag8[itr,j]-precip_44_ctrl_d8_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d9_inst_PRR_lag8[itr,j]-precip_44_ctrl_d9_inst_PR2_lag8[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag8[itr,j]=precip_44_ctrl_d10_inst_PRR_lag8[itr,j]-precip_44_ctrl_d10_inst_PR2_lag8[itr,j]

          if ( j >=9 ) :
               precip_44_ctrl_d2_inst_PR2_lag9[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-9]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag9[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-9]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag9[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-9]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag9[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-9]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag9[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-9]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag9[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-9]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag9[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-9]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag9[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-9]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag9[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-9]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d2_inst_PRR_lag9[itr,j]-precip_44_ctrl_d2_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d3_inst_PRR_lag9[itr,j]-precip_44_ctrl_d3_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d4_inst_PRR_lag9[itr,j]-precip_44_ctrl_d4_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d5_inst_PRR_lag9[itr,j]-precip_44_ctrl_d5_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d6_inst_PRR_lag9[itr,j]-precip_44_ctrl_d6_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d7_inst_PRR_lag9[itr,j]-precip_44_ctrl_d7_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d8_inst_PRR_lag9[itr,j]-precip_44_ctrl_d8_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d9_inst_PRR_lag9[itr,j]-precip_44_ctrl_d9_inst_PR2_lag9[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag9[itr,j]=precip_44_ctrl_d10_inst_PRR_lag9[itr,j]-precip_44_ctrl_d10_inst_PR2_lag9[itr,j]


          if ( j >=10 ) :
               precip_44_ctrl_d2_inst_PR2_lag10[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-10]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag10[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-10]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag10[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-10]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag10[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-10]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag10[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-10]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag10[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-10]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag10[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-10]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag10[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-10]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag10[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-10]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d2_inst_PRR_lag10[itr,j]-precip_44_ctrl_d2_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d3_inst_PRR_lag10[itr,j]-precip_44_ctrl_d3_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d4_inst_PRR_lag10[itr,j]-precip_44_ctrl_d4_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d5_inst_PRR_lag10[itr,j]-precip_44_ctrl_d5_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d6_inst_PRR_lag10[itr,j]-precip_44_ctrl_d6_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d7_inst_PRR_lag10[itr,j]-precip_44_ctrl_d7_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d8_inst_PRR_lag10[itr,j]-precip_44_ctrl_d8_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d9_inst_PRR_lag10[itr,j]-precip_44_ctrl_d9_inst_PR2_lag10[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag10[itr,j]=precip_44_ctrl_d10_inst_PRR_lag10[itr,j]-precip_44_ctrl_d10_inst_PR2_lag10[itr,j]


          if ( j >=11 ) :
               precip_44_ctrl_d2_inst_PR2_lag11[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-11]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag11[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-11]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag11[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-11]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag11[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-11]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag11[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-11]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag11[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-11]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag11[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-11]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag11[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-11]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag11[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-11]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d2_inst_PRR_lag11[itr,j]-precip_44_ctrl_d2_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d3_inst_PRR_lag11[itr,j]-precip_44_ctrl_d3_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d4_inst_PRR_lag11[itr,j]-precip_44_ctrl_d4_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d5_inst_PRR_lag11[itr,j]-precip_44_ctrl_d5_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d6_inst_PRR_lag11[itr,j]-precip_44_ctrl_d6_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d7_inst_PRR_lag11[itr,j]-precip_44_ctrl_d7_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d8_inst_PRR_lag11[itr,j]-precip_44_ctrl_d8_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d9_inst_PRR_lag11[itr,j]-precip_44_ctrl_d9_inst_PR2_lag11[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag11[itr,j]=precip_44_ctrl_d10_inst_PRR_lag11[itr,j]-precip_44_ctrl_d10_inst_PR2_lag11[itr,j]


          if ( j >=12 ) :
               precip_44_ctrl_d2_inst_PR2_lag12[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-12]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag12[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-12]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag12[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-12]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag12[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-12]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag12[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-12]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag12[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-12]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag12[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-12]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag12[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-12]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag12[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-12]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d2_inst_PRR_lag12[itr,j]-precip_44_ctrl_d2_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d3_inst_PRR_lag12[itr,j]-precip_44_ctrl_d3_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d4_inst_PRR_lag12[itr,j]-precip_44_ctrl_d4_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d5_inst_PRR_lag12[itr,j]-precip_44_ctrl_d5_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d6_inst_PRR_lag12[itr,j]-precip_44_ctrl_d6_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d7_inst_PRR_lag12[itr,j]-precip_44_ctrl_d7_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d8_inst_PRR_lag12[itr,j]-precip_44_ctrl_d8_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d9_inst_PRR_lag12[itr,j]-precip_44_ctrl_d9_inst_PR2_lag12[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag12[itr,j]=precip_44_ctrl_d10_inst_PRR_lag12[itr,j]-precip_44_ctrl_d10_inst_PR2_lag12[itr,j]


          if ( j >=13 ) :
               precip_44_ctrl_d2_inst_PR2_lag13[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-13]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag13[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-13]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag13[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-13]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag13[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-13]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag13[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-13]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag13[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-13]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag13[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-13]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag13[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-13]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag13[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-13]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d2_inst_PRR_lag13[itr,j]-precip_44_ctrl_d2_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d3_inst_PRR_lag13[itr,j]-precip_44_ctrl_d3_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d4_inst_PRR_lag13[itr,j]-precip_44_ctrl_d4_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d5_inst_PRR_lag13[itr,j]-precip_44_ctrl_d5_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d6_inst_PRR_lag13[itr,j]-precip_44_ctrl_d6_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d7_inst_PRR_lag13[itr,j]-precip_44_ctrl_d7_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d8_inst_PRR_lag13[itr,j]-precip_44_ctrl_d8_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d9_inst_PRR_lag13[itr,j]-precip_44_ctrl_d9_inst_PR2_lag13[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag13[itr,j]=precip_44_ctrl_d10_inst_PRR_lag13[itr,j]-precip_44_ctrl_d10_inst_PR2_lag13[itr,j]


          if ( j >=14 ) :
               precip_44_ctrl_d2_inst_PR2_lag14[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-14]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag14[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-14]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag14[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-14]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag14[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-14]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag14[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-14]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag14[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-14]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag14[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-14]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag14[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-14]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag14[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-14]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d2_inst_PRR_lag14[itr,j]-precip_44_ctrl_d2_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d3_inst_PRR_lag14[itr,j]-precip_44_ctrl_d3_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d4_inst_PRR_lag14[itr,j]-precip_44_ctrl_d4_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d5_inst_PRR_lag14[itr,j]-precip_44_ctrl_d5_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d6_inst_PRR_lag14[itr,j]-precip_44_ctrl_d6_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d7_inst_PRR_lag14[itr,j]-precip_44_ctrl_d7_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d8_inst_PRR_lag14[itr,j]-precip_44_ctrl_d8_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d9_inst_PRR_lag14[itr,j]-precip_44_ctrl_d9_inst_PR2_lag14[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag14[itr,j]=precip_44_ctrl_d10_inst_PRR_lag14[itr,j]-precip_44_ctrl_d10_inst_PR2_lag14[itr,j]

          if ( j >=15 ) :
               precip_44_ctrl_d2_inst_PR2_lag15[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-15]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag15[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-15]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag15[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-15]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag15[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-15]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag15[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-15]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag15[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-15]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag15[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-15]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag15[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-15]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag15[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-15]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d2_inst_PRR_lag15[itr,j]-precip_44_ctrl_d2_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d3_inst_PRR_lag15[itr,j]-precip_44_ctrl_d3_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d4_inst_PRR_lag15[itr,j]-precip_44_ctrl_d4_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d5_inst_PRR_lag15[itr,j]-precip_44_ctrl_d5_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d6_inst_PRR_lag15[itr,j]-precip_44_ctrl_d6_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d7_inst_PRR_lag15[itr,j]-precip_44_ctrl_d7_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d8_inst_PRR_lag15[itr,j]-precip_44_ctrl_d8_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d9_inst_PRR_lag15[itr,j]-precip_44_ctrl_d9_inst_PR2_lag15[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag15[itr,j]=precip_44_ctrl_d10_inst_PRR_lag15[itr,j]-precip_44_ctrl_d10_inst_PR2_lag15[itr,j]


          if ( j >=16 ) :
               precip_44_ctrl_d2_inst_PR2_lag16[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-16]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag16[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-16]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag16[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-16]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag16[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-16]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag16[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-16]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag16[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-16]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag16[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-16]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag16[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-16]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag16[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-16]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d2_inst_PRR_lag16[itr,j]-precip_44_ctrl_d2_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d3_inst_PRR_lag16[itr,j]-precip_44_ctrl_d3_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d4_inst_PRR_lag16[itr,j]-precip_44_ctrl_d4_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d5_inst_PRR_lag16[itr,j]-precip_44_ctrl_d5_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d6_inst_PRR_lag16[itr,j]-precip_44_ctrl_d6_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d7_inst_PRR_lag16[itr,j]-precip_44_ctrl_d7_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d8_inst_PRR_lag16[itr,j]-precip_44_ctrl_d8_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d9_inst_PRR_lag16[itr,j]-precip_44_ctrl_d9_inst_PR2_lag16[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag16[itr,j]=precip_44_ctrl_d10_inst_PRR_lag16[itr,j]-precip_44_ctrl_d10_inst_PR2_lag16[itr,j]

          if ( j >=17 ) :
               precip_44_ctrl_d2_inst_PR2_lag17[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-17]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag17[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-17]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag17[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-17]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag17[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-17]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag17[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-17]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag17[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-17]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag17[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-17]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag17[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-17]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag17[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-17]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d2_inst_PRR_lag17[itr,j]-precip_44_ctrl_d2_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d3_inst_PRR_lag17[itr,j]-precip_44_ctrl_d3_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d4_inst_PRR_lag17[itr,j]-precip_44_ctrl_d4_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d5_inst_PRR_lag17[itr,j]-precip_44_ctrl_d5_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d6_inst_PRR_lag17[itr,j]-precip_44_ctrl_d6_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d7_inst_PRR_lag17[itr,j]-precip_44_ctrl_d7_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d8_inst_PRR_lag17[itr,j]-precip_44_ctrl_d8_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d9_inst_PRR_lag17[itr,j]-precip_44_ctrl_d9_inst_PR2_lag17[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag17[itr,j]=precip_44_ctrl_d10_inst_PRR_lag17[itr,j]-precip_44_ctrl_d10_inst_PR2_lag17[itr,j]

   

          if ( j >=18 ) :
               precip_44_ctrl_d2_inst_PR2_lag18[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-18]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag18[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-18]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag18[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-18]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag18[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-18]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag18[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-18]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag18[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-18]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag18[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-18]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag18[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-18]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag18[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-18]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d2_inst_PRR_lag18[itr,j]-precip_44_ctrl_d2_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d3_inst_PRR_lag18[itr,j]-precip_44_ctrl_d3_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d4_inst_PRR_lag18[itr,j]-precip_44_ctrl_d4_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d5_inst_PRR_lag18[itr,j]-precip_44_ctrl_d5_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d6_inst_PRR_lag18[itr,j]-precip_44_ctrl_d6_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d7_inst_PRR_lag18[itr,j]-precip_44_ctrl_d7_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d8_inst_PRR_lag18[itr,j]-precip_44_ctrl_d8_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d9_inst_PRR_lag18[itr,j]-precip_44_ctrl_d9_inst_PR2_lag18[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag18[itr,j]=precip_44_ctrl_d10_inst_PRR_lag18[itr,j]-precip_44_ctrl_d10_inst_PR2_lag18[itr,j]

 

          if ( j >=19 ) :
               precip_44_ctrl_d2_inst_PR2_lag19[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-19]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag19[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-19]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag19[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-19]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag19[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-19]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag19[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-19]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag19[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-19]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag19[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-19]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag19[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-19]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag19[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-19]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d2_inst_PRR_lag19[itr,j]-precip_44_ctrl_d2_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d3_inst_PRR_lag19[itr,j]-precip_44_ctrl_d3_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d4_inst_PRR_lag19[itr,j]-precip_44_ctrl_d4_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d5_inst_PRR_lag19[itr,j]-precip_44_ctrl_d5_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d6_inst_PRR_lag19[itr,j]-precip_44_ctrl_d6_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d7_inst_PRR_lag19[itr,j]-precip_44_ctrl_d7_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d8_inst_PRR_lag19[itr,j]-precip_44_ctrl_d8_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d9_inst_PRR_lag19[itr,j]-precip_44_ctrl_d9_inst_PR2_lag19[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag19[itr,j]=precip_44_ctrl_d10_inst_PRR_lag19[itr,j]-precip_44_ctrl_d10_inst_PR2_lag19[itr,j]



          if ( j >=20 ) :
               precip_44_ctrl_d2_inst_PR2_lag20[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-20]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag20[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-20]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag20[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-20]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag20[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-20]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag20[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-20]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag20[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-20]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag20[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-20]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag20[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-20]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag20[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-20]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d2_inst_PRR_lag20[itr,j]-precip_44_ctrl_d2_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d3_inst_PRR_lag20[itr,j]-precip_44_ctrl_d3_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d4_inst_PRR_lag20[itr,j]-precip_44_ctrl_d4_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d5_inst_PRR_lag20[itr,j]-precip_44_ctrl_d5_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d6_inst_PRR_lag20[itr,j]-precip_44_ctrl_d6_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d7_inst_PRR_lag20[itr,j]-precip_44_ctrl_d7_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d8_inst_PRR_lag20[itr,j]-precip_44_ctrl_d8_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d9_inst_PRR_lag20[itr,j]-precip_44_ctrl_d9_inst_PR2_lag20[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag20[itr,j]=precip_44_ctrl_d10_inst_PRR_lag20[itr,j]-precip_44_ctrl_d10_inst_PR2_lag20[itr,j]

          if ( j >=21 ) :
               precip_44_ctrl_d2_inst_PR2_lag21[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-21]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag21[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-21]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag21[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-21]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag21[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-21]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag21[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-21]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag21[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-21]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag21[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-21]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag21[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-21]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag21[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-21]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d2_inst_PRR_lag21[itr,j]-precip_44_ctrl_d2_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d3_inst_PRR_lag21[itr,j]-precip_44_ctrl_d3_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d4_inst_PRR_lag21[itr,j]-precip_44_ctrl_d4_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d5_inst_PRR_lag21[itr,j]-precip_44_ctrl_d5_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d6_inst_PRR_lag21[itr,j]-precip_44_ctrl_d6_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d7_inst_PRR_lag21[itr,j]-precip_44_ctrl_d7_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d8_inst_PRR_lag21[itr,j]-precip_44_ctrl_d8_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d9_inst_PRR_lag21[itr,j]-precip_44_ctrl_d9_inst_PR2_lag21[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag21[itr,j]=precip_44_ctrl_d10_inst_PRR_lag21[itr,j]-precip_44_ctrl_d10_inst_PR2_lag21[itr,j]


          if ( j >=22 ) :
               precip_44_ctrl_d2_inst_PR2_lag22[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-22]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag22[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-22]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag22[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-22]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag22[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-22]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag22[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-22]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag22[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-22]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag22[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-22]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag22[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-22]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag22[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-22]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d2_inst_PRR_lag22[itr,j]-precip_44_ctrl_d2_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d3_inst_PRR_lag22[itr,j]-precip_44_ctrl_d3_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d4_inst_PRR_lag22[itr,j]-precip_44_ctrl_d4_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d5_inst_PRR_lag22[itr,j]-precip_44_ctrl_d5_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d6_inst_PRR_lag22[itr,j]-precip_44_ctrl_d6_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d7_inst_PRR_lag22[itr,j]-precip_44_ctrl_d7_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d8_inst_PRR_lag22[itr,j]-precip_44_ctrl_d8_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d9_inst_PRR_lag22[itr,j]-precip_44_ctrl_d9_inst_PR2_lag22[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag22[itr,j]=precip_44_ctrl_d10_inst_PRR_lag22[itr,j]-precip_44_ctrl_d10_inst_PR2_lag22[itr,j]


          if ( j >=23 ) :
               precip_44_ctrl_d2_inst_PR2_lag23[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-23]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag23[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-23]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag23[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-23]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag23[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-23]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag23[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-23]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag23[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-23]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag23[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-23]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag23[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-23]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag23[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-23]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d2_inst_PRR_lag23[itr,j]-precip_44_ctrl_d2_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d3_inst_PRR_lag23[itr,j]-precip_44_ctrl_d3_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d4_inst_PRR_lag23[itr,j]-precip_44_ctrl_d4_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d5_inst_PRR_lag23[itr,j]-precip_44_ctrl_d5_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d6_inst_PRR_lag23[itr,j]-precip_44_ctrl_d6_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d7_inst_PRR_lag23[itr,j]-precip_44_ctrl_d7_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d8_inst_PRR_lag23[itr,j]-precip_44_ctrl_d8_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d9_inst_PRR_lag23[itr,j]-precip_44_ctrl_d9_inst_PR2_lag23[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag23[itr,j]=precip_44_ctrl_d10_inst_PRR_lag23[itr,j]-precip_44_ctrl_d10_inst_PR2_lag23[itr,j]

  

          if ( j >=24 ) :
               precip_44_ctrl_d2_inst_PR2_lag24[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-24]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag24[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-24]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag24[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-24]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag24[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-24]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag24[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-24]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag24[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-24]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag24[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-24]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag24[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-24]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag24[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-24]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d2_inst_PRR_lag24[itr,j]-precip_44_ctrl_d2_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d3_inst_PRR_lag24[itr,j]-precip_44_ctrl_d3_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d4_inst_PRR_lag24[itr,j]-precip_44_ctrl_d4_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d5_inst_PRR_lag24[itr,j]-precip_44_ctrl_d5_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d6_inst_PRR_lag24[itr,j]-precip_44_ctrl_d6_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d7_inst_PRR_lag24[itr,j]-precip_44_ctrl_d7_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d8_inst_PRR_lag24[itr,j]-precip_44_ctrl_d8_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d9_inst_PRR_lag24[itr,j]-precip_44_ctrl_d9_inst_PR2_lag24[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag24[itr,j]=precip_44_ctrl_d10_inst_PRR_lag24[itr,j]-precip_44_ctrl_d10_inst_PR2_lag24[itr,j]

 

          if ( j >=25 ) :
               precip_44_ctrl_d2_inst_PR2_lag25[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-25]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag25[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-25]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag25[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-25]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag25[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-25]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag25[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-25]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag25[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-25]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag25[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-25]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag25[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-25]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag25[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-25]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d2_inst_PRR_lag25[itr,j]-precip_44_ctrl_d2_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d3_inst_PRR_lag25[itr,j]-precip_44_ctrl_d3_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d4_inst_PRR_lag25[itr,j]-precip_44_ctrl_d4_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d5_inst_PRR_lag25[itr,j]-precip_44_ctrl_d5_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d6_inst_PRR_lag25[itr,j]-precip_44_ctrl_d6_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d7_inst_PRR_lag25[itr,j]-precip_44_ctrl_d7_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d8_inst_PRR_lag25[itr,j]-precip_44_ctrl_d8_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d9_inst_PRR_lag25[itr,j]-precip_44_ctrl_d9_inst_PR2_lag25[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag25[itr,j]=precip_44_ctrl_d10_inst_PRR_lag25[itr,j]-precip_44_ctrl_d10_inst_PR2_lag25[itr,j]


          if ( j >=26 ) :
               precip_44_ctrl_d2_inst_PR2_lag26[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-26]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag26[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-26]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag26[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-26]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag26[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-26]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag26[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-26]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag26[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-26]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag26[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-26]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag26[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-26]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag26[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-26]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d2_inst_PRR_lag26[itr,j]-precip_44_ctrl_d2_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d3_inst_PRR_lag26[itr,j]-precip_44_ctrl_d3_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d4_inst_PRR_lag26[itr,j]-precip_44_ctrl_d4_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d5_inst_PRR_lag26[itr,j]-precip_44_ctrl_d5_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d6_inst_PRR_lag26[itr,j]-precip_44_ctrl_d6_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d7_inst_PRR_lag26[itr,j]-precip_44_ctrl_d7_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d8_inst_PRR_lag26[itr,j]-precip_44_ctrl_d8_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d9_inst_PRR_lag26[itr,j]-precip_44_ctrl_d9_inst_PR2_lag26[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag26[itr,j]=precip_44_ctrl_d10_inst_PRR_lag26[itr,j]-precip_44_ctrl_d10_inst_PR2_lag26[itr,j]


          if ( j >=27 ) :
               precip_44_ctrl_d2_inst_PR2_lag27[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-27]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag27[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-27]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag27[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-27]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag27[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-27]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag27[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-27]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag27[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-27]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag27[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-27]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag27[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-27]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag27[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-27]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d2_inst_PRR_lag27[itr,j]-precip_44_ctrl_d2_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d3_inst_PRR_lag27[itr,j]-precip_44_ctrl_d3_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d4_inst_PRR_lag27[itr,j]-precip_44_ctrl_d4_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d5_inst_PRR_lag27[itr,j]-precip_44_ctrl_d5_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d6_inst_PRR_lag27[itr,j]-precip_44_ctrl_d6_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d7_inst_PRR_lag27[itr,j]-precip_44_ctrl_d7_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d8_inst_PRR_lag27[itr,j]-precip_44_ctrl_d8_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d9_inst_PRR_lag27[itr,j]-precip_44_ctrl_d9_inst_PR2_lag27[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag27[itr,j]=precip_44_ctrl_d10_inst_PRR_lag27[itr,j]-precip_44_ctrl_d10_inst_PR2_lag27[itr,j]


          if ( j >=28 ) :
               precip_44_ctrl_d2_inst_PR2_lag28[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-28]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag28[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-28]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag28[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-28]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag28[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-28]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag28[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-28]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag28[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-28]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag28[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-28]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag28[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-28]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag28[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-28]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d2_inst_PRR_lag28[itr,j]-precip_44_ctrl_d2_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d3_inst_PRR_lag28[itr,j]-precip_44_ctrl_d3_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d4_inst_PRR_lag28[itr,j]-precip_44_ctrl_d4_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d5_inst_PRR_lag28[itr,j]-precip_44_ctrl_d5_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d6_inst_PRR_lag28[itr,j]-precip_44_ctrl_d6_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d7_inst_PRR_lag28[itr,j]-precip_44_ctrl_d7_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d8_inst_PRR_lag28[itr,j]-precip_44_ctrl_d8_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d9_inst_PRR_lag28[itr,j]-precip_44_ctrl_d9_inst_PR2_lag28[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag28[itr,j]=precip_44_ctrl_d10_inst_PRR_lag28[itr,j]-precip_44_ctrl_d10_inst_PR2_lag28[itr,j]


          if ( j >=29 ) :
               precip_44_ctrl_d2_inst_PR2_lag29[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-29]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag29[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-29]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag29[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-29]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag29[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-29]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag29[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-29]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag29[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-29]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag29[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-29]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag29[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-29]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag29[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-29]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d2_inst_PRR_lag29[itr,j]-precip_44_ctrl_d2_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d3_inst_PRR_lag29[itr,j]-precip_44_ctrl_d3_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d4_inst_PRR_lag29[itr,j]-precip_44_ctrl_d4_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d5_inst_PRR_lag29[itr,j]-precip_44_ctrl_d5_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d6_inst_PRR_lag29[itr,j]-precip_44_ctrl_d6_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d7_inst_PRR_lag29[itr,j]-precip_44_ctrl_d7_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d8_inst_PRR_lag29[itr,j]-precip_44_ctrl_d8_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d9_inst_PRR_lag29[itr,j]-precip_44_ctrl_d9_inst_PR2_lag29[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag29[itr,j]=precip_44_ctrl_d10_inst_PRR_lag29[itr,j]-precip_44_ctrl_d10_inst_PR2_lag29[itr,j]


          if ( j >=30 ) :
               precip_44_ctrl_d2_inst_PR2_lag30[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-30]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag30[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-30]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag30[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-30]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag30[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-30]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag30[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-30]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag30[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-30]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag30[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-30]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag30[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-30]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag30[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-30]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d2_inst_PRR_lag30[itr,j]-precip_44_ctrl_d2_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d3_inst_PRR_lag30[itr,j]-precip_44_ctrl_d3_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d4_inst_PRR_lag30[itr,j]-precip_44_ctrl_d4_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d5_inst_PRR_lag30[itr,j]-precip_44_ctrl_d5_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d6_inst_PRR_lag30[itr,j]-precip_44_ctrl_d6_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d7_inst_PRR_lag30[itr,j]-precip_44_ctrl_d7_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d8_inst_PRR_lag30[itr,j]-precip_44_ctrl_d8_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d9_inst_PRR_lag30[itr,j]-precip_44_ctrl_d9_inst_PR2_lag30[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag30[itr,j]=precip_44_ctrl_d10_inst_PRR_lag30[itr,j]-precip_44_ctrl_d10_inst_PR2_lag30[itr,j]



          if ( j >=31 ) :
               precip_44_ctrl_d2_inst_PR2_lag31[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-31]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag31[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-31]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag31[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-31]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag31[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-31]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag31[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-31]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag31[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-31]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag31[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-31]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag31[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-31]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag31[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-31]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d2_inst_PRR_lag31[itr,j]-precip_44_ctrl_d2_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d3_inst_PRR_lag31[itr,j]-precip_44_ctrl_d3_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d4_inst_PRR_lag31[itr,j]-precip_44_ctrl_d4_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d5_inst_PRR_lag31[itr,j]-precip_44_ctrl_d5_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d6_inst_PRR_lag31[itr,j]-precip_44_ctrl_d6_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d7_inst_PRR_lag31[itr,j]-precip_44_ctrl_d7_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d8_inst_PRR_lag31[itr,j]-precip_44_ctrl_d8_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d9_inst_PRR_lag31[itr,j]-precip_44_ctrl_d9_inst_PR2_lag31[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag31[itr,j]=precip_44_ctrl_d10_inst_PRR_lag31[itr,j]-precip_44_ctrl_d10_inst_PR2_lag31[itr,j]

          if ( j >=32 ) :
               precip_44_ctrl_d2_inst_PR2_lag32[itr,j]=precip_44_ctrl_d2_inst_PR_lag[itr,j-32]*precip_44_ctrl_d2_inst_PR_lag[itr,j]
               precip_44_ctrl_d3_inst_PR2_lag32[itr,j]=precip_44_ctrl_d3_inst_PR_lag[itr,j-32]*precip_44_ctrl_d3_inst_PR_lag[itr,j]
               precip_44_ctrl_d4_inst_PR2_lag32[itr,j]=precip_44_ctrl_d4_inst_PR_lag[itr,j-32]*precip_44_ctrl_d4_inst_PR_lag[itr,j]
               precip_44_ctrl_d5_inst_PR2_lag32[itr,j]=precip_44_ctrl_d5_inst_PR_lag[itr,j-32]*precip_44_ctrl_d5_inst_PR_lag[itr,j]
               precip_44_ctrl_d6_inst_PR2_lag32[itr,j]=precip_44_ctrl_d6_inst_PR_lag[itr,j-32]*precip_44_ctrl_d6_inst_PR_lag[itr,j]
               precip_44_ctrl_d7_inst_PR2_lag32[itr,j]=precip_44_ctrl_d7_inst_PR_lag[itr,j-32]*precip_44_ctrl_d7_inst_PR_lag[itr,j]
               precip_44_ctrl_d8_inst_PR2_lag32[itr,j]=precip_44_ctrl_d8_inst_PR_lag[itr,j-32]*precip_44_ctrl_d8_inst_PR_lag[itr,j]
               precip_44_ctrl_d9_inst_PR2_lag32[itr,j]=precip_44_ctrl_d9_inst_PR_lag[itr,j-32]*precip_44_ctrl_d9_inst_PR_lag[itr,j]
               precip_44_ctrl_d10_inst_PR2_lag32[itr,j]=precip_44_ctrl_d10_inst_PR_lag[itr,j-32]*precip_44_ctrl_d10_inst_PR_lag[itr,j]
          
               precip_44_ctrl_d2_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d2_inst_PRR_lag32[itr,j]-precip_44_ctrl_d2_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d3_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d3_inst_PRR_lag32[itr,j]-precip_44_ctrl_d3_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d4_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d4_inst_PRR_lag32[itr,j]-precip_44_ctrl_d4_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d5_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d5_inst_PRR_lag32[itr,j]-precip_44_ctrl_d5_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d6_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d6_inst_PRR_lag32[itr,j]-precip_44_ctrl_d6_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d7_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d7_inst_PRR_lag32[itr,j]-precip_44_ctrl_d7_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d8_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d8_inst_PRR_lag32[itr,j]-precip_44_ctrl_d8_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d9_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d9_inst_PRR_lag32[itr,j]-precip_44_ctrl_d9_inst_PR2_lag32[itr,j]
               precip_44_ctrl_d10_inst_PRR_PR2_lag32[itr,j]=precip_44_ctrl_d10_inst_PRR_lag32[itr,j]-precip_44_ctrl_d10_inst_PR2_lag32[itr,j]

proRR_ProR2t_lag_dt=N.zeros((nlag, n_30min))
proRR_ProR2t_lag_dt[0,:]=float('nan')
for j in N.arange(48):
     proRR_ProR2t_lag_dt[1,j]=proRR_proR2_lag1_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[2,j]=proRR_proR2_lag2_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[3,j]=proRR_proR2_lag3_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[4,j]=proRR_proR2_lag4_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[5,j]=proRR_proR2_lag5_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[6,j]=proRR_proR2_lag6_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[7,j]=proRR_proR2_lag7_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[8,j]=proRR_proR2_lag8_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[9,j]=proRR_proR2_lag9_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[10,j]=proRR_proR2_lag10_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[11,j]=proRR_proR2_lag11_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[12,j]=proRR_proR2_lag12_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[13,j]=proRR_proR2_lag13_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[14,j]=proRR_proR2_lag14_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[15,j]=proRR_proR2_lag15_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[16,j]=proRR_proR2_lag16_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[17,j]=proRR_proR2_lag17_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[18,j]=proRR_proR2_lag18_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[19,j]=proRR_proR2_lag19_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[20,j]=proRR_proR2_lag20_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[21,j]=proRR_proR2_lag21_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[22,j]=proRR_proR2_lag22_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[23,j]=proRR_proR2_lag23_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[24,j]=proRR_proR2_lag24_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[25,j]=proRR_proR2_lag25_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[26,j]=proRR_proR2_lag26_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[27,j]=proRR_proR2_lag27_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[28,j]=proRR_proR2_lag28_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[29,j]=proRR_proR2_lag29_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[30,j]=proRR_proR2_lag30_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[31,j]=proRR_proR2_lag31_ens_mean_inst_44_vs_AveRR_thr[0,j]
     proRR_ProR2t_lag_dt[32,j]=proRR_proR2_lag32_ens_mean_inst_44_vs_AveRR_thr[0,j]
     if ( j <=11):
          proRR_ProR2t_lag_dt[:,j]=float('nan')
     if ( j == 12):
          proRR_ProR2t_lag_dt[5:nlag,j]=float('nan')
     if ( j == 13):
          proRR_ProR2t_lag_dt[6:nlag,j]=float('nan')
     if ( j == 14):
          proRR_ProR2t_lag_dt[7:nlag,j]=float('nan')
     if ( j == 15):
          proRR_ProR2t_lag_dt[8:nlag,j]=float('nan')
     if ( j == 16):
          proRR_ProR2t_lag_dt[9:nlag,j]=float('nan')
     if ( j == 17):
          proRR_ProR2t_lag_dt[10:nlag,j]=float('nan')
     if ( j == 18):
          proRR_ProR2t_lag_dt[11:nlag,j]=float('nan')
     if ( j == 19):
          proRR_ProR2t_lag_dt[12:nlag,j]=float('nan')
     if ( j == 20):
          proRR_ProR2t_lag_dt[13:nlag,j]=float('nan')
     if ( j == 21):
          proRR_ProR2t_lag_dt[14:nlag,j]=float('nan')
     if ( j == 22):
          proRR_ProR2t_lag_dt[15:nlag,j]=float('nan')
     if ( j == 23):
          proRR_ProR2t_lag_dt[16:nlag,j]=float('nan')
     if ( j == 24):
          proRR_ProR2t_lag_dt[17:nlag,j]=float('nan')
     if ( j == 25):
          proRR_ProR2t_lag_dt[18:nlag,j]=float('nan')
     if ( j == 26):
          proRR_ProR2t_lag_dt[19:nlag,j]=float('nan')
     if ( j == 27):
          proRR_ProR2t_lag_dt[20:nlag,j]=float('nan')
     if ( j == 28):
          proRR_ProR2t_lag_dt[21:nlag,j]=float('nan')
     if ( j == 29):
          proRR_ProR2t_lag_dt[22:nlag,j]=float('nan')
     if ( j == 30):
          proRR_ProR2t_lag_dt[23:nlag,j]=float('nan')
     if ( j == 31):
          proRR_ProR2t_lag_dt[24:nlag,j]=float('nan')
     if ( j == 32):
          proRR_ProR2t_lag_dt[25:nlag,j]=float('nan')
     if ( j == 33):
          proRR_ProR2t_lag_dt[26:nlag,j]=float('nan')
     if ( j == 34):
          proRR_ProR2t_lag_dt[27:nlag,j]=float('nan')
     if ( j == 35):
          proRR_ProR2t_lag_dt[28:nlag,j]=float('nan')
     if ( j == 36):
          proRR_ProR2t_lag_dt[29:nlag,j]=float('nan')
     if ( j == 37):
          proRR_ProR2t_lag_dt[30:nlag,j]=float('nan')
     if ( j == 38):
          proRR_ProR2t_lag_dt[31:nlag,j]=float('nan')
     if ( j == 39):
          proRR_ProR2t_lag_dt[32:nlag,j]=float('nan')
     if ( j == 40):
          proRR_ProR2t_lag_dt[33:nlag,j]=float('nan')
     if ( j> 40):
          proRR_ProR2t_lag_dt[34:nlag,j]=float('nan')


np.save('proRR_ProR2t_lag_dt_CTRL',proRR_ProR2t_lag_dt)

