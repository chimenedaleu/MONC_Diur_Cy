# Import libraries
from mpl_toolkits.basemap import Basemap, cm
from scipy import interpolate
import numpy as N 
import matplotlib.pyplot as plt
import os 
import netCDF4
import matplotlib.patches as patches
from PIL import Image
import numpy
import matplotlib
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from matplotlib.ticker import FormatStrFormatter
#import matplotlib.pyplot as plt

nx=512
ny=512
nz=99
dx=200
dy=200
length_x=N.zeros(nx)
length_y=N.zeros(ny)
#th_3d_xy_contr_d2_Day_mean=N.zeros((nx,ny,nz))


L_vap                = 2.501e6     #J/kg
C_p                  = 1005.0      #J/kgK
R                    = 287.058 
kappa                = R/C_p    
gra                  = 9.8
L_sub=2.85e6 #2.834e6
pre_o=100000.0
# Setup file paths 


inpath = os.path.join( 'projects','nexcs-n06','cdaleu', 'MONC_DiurCy', 'vn0.8_DuirCir','diagnostic_files')
outpath = os.path.join( 'projects','nexcs-n06','cdaleu', 'PYTHON_NEW')



nrun=1
neq=5
nz=99




ncond=16
ndiag=10
n1_10min=90
n2_10min=92
n3_10min=90
n4_10min=91
n5_10min=92
n6_10min=92
n7_10min=93
n8_10min=92
n9_10min=92
n10_10min=93
n_10min=n1_10min+n2_10min+n3_10min+n4_10min+n5_10min+n6_10min+n7_10min+n8_10min+n9_10min+n10_10min+10
n1_30min=46
n2_30min=47
n3_30min=47
n4_30min=46
n5_30min=47
n6_30min=47
n7_30min=47
n8_30min=47
n9_30min=47
n10_30min=47
n_30min=n1_30min+n2_30min+n3_30min+n4_30min+n5_30min+n6_30min+n7_30min+n8_30min+n9_30min+n10_30min+10


pref=N.zeros((nz))
AC_MF_time=N.zeros((n_10min,nz))
AC_frac_time=N.zeros((n_10min,nz))
ACuzt_frac_time=N.zeros((nz,n_10min))
ACzt_MF_time=N.zeros((nz,n_10min))
ACuzt_MF_time=N.zeros((nz,n_10min))
ACzt_frac_time=N.zeros((nz,n_10min))

ACu_MF_time=N.zeros((n_10min,nz))
ACu_frac_time=N.zeros((n_10min,nz))
ACuS_MF_time=N.zeros((n_10min,nz))
ACuS_frac_time=N.zeros((n_10min,nz))

th_mean=N.zeros((n_30min,nz))
qv_mean=N.zeros((n_30min,nz))
cloud_top_max=N.zeros((n_10min))
cloud_base_max=N.zeros((n_10min))
cloudBCu_top_max=N.zeros((n_10min))
cloudBCu_base_max=N.zeros((n_10min))

cloud_MF_time_cb5     =N.zeros((n_10min))
cloud_frac_time_cb5=N.zeros((n_10min))
cloudu_MF_time_cb5     =N.zeros((n_10min))
cloudu_frac_time_cb5=N.zeros((n_10min))
clouduS_MF_time_cb5     =N.zeros((n_10min))
clouduS_frac_time_cb5=N.zeros((n_10min))
cloudBCu_MF_time_cb5     =N.zeros((n_10min))
cloudBCu_frac_time_cb5=N.zeros((n_10min))

cloud_MF_time_cb4     =N.zeros((n_10min))
cloud_frac_time_cb4=N.zeros((n_10min))
cloud_MF_time_cb3     =N.zeros((n_10min))
cloud_frac_time_cb3=N.zeros((n_10min))
cloud_MF_time_cb2     =N.zeros((n_10min))
cloud_frac_time_cb2=N.zeros((n_10min))
cloud_MF_time_cb1     =N.zeros((n_10min))
cloud_frac_time_cb1=N.zeros((n_10min))
cloud_MF_time_cb     =N.zeros((n_10min))
cloud_frac_time_cb=N.zeros((n_10min))

BCu_w_time=N.zeros((n_10min,nz))
BCu_MF_time=N.zeros((n_10min,nz))
BCuzt_frac_time=N.zeros((nz,n_10min))
BCuzt_MF_time=N.zeros((nz,n_10min))
BCu_frac_time=N.zeros((n_10min,nz))

surf_precip_xy_time=N.zeros((n_10min,nx,ny))
modets_time_1h=N.zeros((n_30min))
surf_precip_time=N.zeros((n_10min))
den_n_mean=N.zeros((nz))
hei=N.zeros((nz))
hein=N.zeros((nz))
den_mean=N.zeros((nz))
th_ref_prof=N.zeros((nz))
pre_ref_prof=N.zeros((nz))


shf_time=N.zeros((n_10min))
surf_flx=N.zeros((n_10min))
lhf_time=N.zeros((n_10min))
modets_time=N.zeros((n_10min))

th_std_time=N.zeros((n_30min,nz))
thzt_std_time=N.zeros((nz,n_30min))
qv_std_time=N.zeros((n_30min,nz))
qvzt_std_time=N.zeros((nz,n_30min))


thzy14_snap=N.zeros((nz,ny))
thzy24_snap=N.zeros((nz,ny))
qvzy14_snap=N.zeros((nz,ny))
qvzy24_snap=N.zeros((nz,ny))

dy=200
length_y[1]=200.
for j in N.arange(ny-1):
     length_y[j+1]=length_y[j]+dy
print length_y[0:10]




cwv_cond=0
for a in N.arange(nrun):
  #   a=a+1
     if (a+1 < 11):
          infile = '/'+inpath+'/newRC_CTRL/DC_newRC_CTRL_'+str((a+1)*864)+'00.0.nc'
          print infile
          print '.nc', a+1
          data = netCDF4.Dataset(infile)



# Height, pressure, th etc
     z  = data.variables['z'][:]
     zn = data.variables['zn'][:]
     den = data.variables['rho'][:,:] 
     den_n = data.variables['rhon'][:,:]
     th_ref = data.variables['thref'][:,:]
     pre_ref = data.variables['prefn'][:,:]
  
# Extract mean variables
     umean = data.variables['u_wind_mean'][:,:]
     vmean = data.variables['v_wind_mean'][:,:]
     thmean = data.variables['theta_mean'][:,:]

     vapmean = data.variables['vapour_mmr_mean'][:,:]
     CondDiags=data.variables['CondDiags_mean'][:,:]
     th = data.variables['th'][:,:]
     qv = data.variables['q_vapour'][:,:]

     surf_precip = data.variables['surface_precip'][:,:]
#Instantaneous variables
     vwp = data.variables['VWP_mean'][:]  #mm

     shf = data.variables['senhf_mean'][:]  #watt/m2
     lhf = data.variables['lathf_mean'][:] #watt/m2
     mod_ts= data.variables['time_series_50_900.0'][:] #
     mod_ts_1h= data.variables['time_series_50_1800.0'][:] #

     den_n_mean[:]=den_n[1,:]
     den_mean[:]=den[1,:]
     hei[:]=z[1,:]
     hein[:]=zn[:]
     th_ref_prof[:]=th_ref[10,:]
     pref[:]=pre_ref[10,:]


     if (a+1 == 1):
          for i in N.arange(n1_10min):

               shf_time[i]=shf[i]
               lhf_time[i]=lhf[i]
               modets_time[i]= mod_ts[i]/3600.
               surf_precip_xy_time[i,:,:]=surf_precip[i,:,:]
               surf_precip_time[i]=surf_precip[i,:,:].mean()

               BCu_frac_time[i,:]     = CondDiags[i,0,2,:]
               BCuzt_frac_time[:,i]     = BCu_frac_time[i,:] 

               AC_frac_time[i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[i,:]     = CondDiags[i,1,4,:]
               ACu_frac_time[i,:]     = CondDiags[i,0,5,:]
               ACuzt_frac_time[:,i]     = ACu_frac_time[i,:] 
               ACzt_frac_time[:,i]     = AC_frac_time[i,:] 

               ACu_MF_time[i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[i,:]     = CondDiags[i,1,6,:]

               for k in N.arange(nz):
                    BCu_MF_time[i,k]= BCu_MF_time[i,k]*den_n_mean[k]
                    BCuzt_MF_time[k,i]     = BCu_MF_time[i,k] 


                    AC_MF_time[i,k]= AC_MF_time[i,k]*den_n_mean[k]
                    ACzt_MF_time[k,i]     = AC_MF_time[i,k] 
                    ACu_MF_time[i,k]= ACu_MF_time[i,k]*den_n_mean[k]
                    ACuzt_MF_time[k,i]     = ACu_MF_time[i,k] 
                    ACuS_MF_time[i,k]= ACuS_MF_time[i,k]*den_n_mean[k]
 

                    if (AC_frac_time[i,k] > 0.0):
                         cloud_top_max[i]=hein[k]
                    if (BCu_frac_time[i,k] > 0.0):
                         cloudBCu_top_max[i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[i,k] > 0.0):
                         cloud_base_max[i]=hein[k]
                    if (BCu_frac_time[i,k] > 0.0):
                         cloudBCu_base_max[i]=hein[k]

          for i in N.arange(n1_30min):
               modets_time_1h[i]= mod_ts_1h[i]/3600.

               for k in N.arange(nz):
                    th_mean[i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[i,k]=1000.*qv[i,:,:,k].mean()

                    th_std_time[i,k]=np.std(th[i,:,:,k])
                    thzt_std_time[k,i]=th_std_time[i,k]
                    qv_std_time[i,k]=np.std(1000.*qv[i,:,:,k])
                    qvzt_std_time[k,i]=qv_std_time[i,k]

                    thzy14_snap[k,j]=th[27,256,j,k]-th[27,:,:,k].mean() 
                    thzy24_snap[k,j]=th[45,256,j,k]-th[45,:,:,k].mean() 
                    qvzy14_snap[k,j]=1000.*(qv[27,256,j,k]-qv[27,:,:,k].mean())  
                    qvzy24_snap[k,j]=1000.*(qv[45,256,j,k]-qv[45,:,:,k].mean()) 


     if (a+1 == 2):
          for i in N.arange(n2_10min):
               surf_precip_xy_time[n1_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n1_10min+i]=surf_precip[i,:,:].mean()

               shf_time[n1_10min+i]=shf[i]
               lhf_time[n1_10min+i]=lhf[i]
               modets_time[n1_10min+i]= mod_ts[i]/3600.

               BCu_frac_time[n1_10min+i,:]     = CondDiags[i,0,2,:]
               BCuzt_frac_time[:,n1_10min+i]     = BCu_frac_time[n1_10min+i,:] 

               BCu_MF_time[n1_10min+i,:]     = CondDiags[i,1,2,:]

               AC_frac_time[n1_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n1_10min+i,:]     = CondDiags[i,1,4,:]
               ACzt_frac_time[:,n1_10min+i]     = AC_frac_time[n1_10min+i,:] 
               ACu_frac_time[n1_10min+i,:]     = CondDiags[i,0,5,:]
               ACuzt_frac_time[:,n1_10min+i]     = ACu_frac_time[n1_10min+i,:] 
               ACu_MF_time[n1_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n1_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n1_10min+i,:]     = CondDiags[i,1,6,:]

               for k in N.arange(nz):
                    BCu_MF_time[n1_10min+i,k]= BCu_MF_time[n1_10min+i,k]*den_n_mean[k]
                    BCuzt_MF_time[k,n1_10min+i]     = BCu_MF_time[n1_10min+i,k] 
                   
                    AC_MF_time[n1_10min+i,k]= AC_MF_time[n1_10min+i,k]*den_n_mean[k]
                    ACzt_MF_time[k,n1_10min+i]     = AC_MF_time[n1_10min+i,k] 
                    ACu_MF_time[n1_10min+i,k]= ACu_MF_time[n1_10min+i,k]*den_n_mean[k]
                    ACuzt_MF_time[k,n1_10min+i]     = ACu_MF_time[n1_10min+i,k] 
                    ACuS_MF_time[n1_10min+i,k]= ACuS_MF_time[n1_10min+i,k]*den_n_mean[k]

                    if (AC_frac_time[n1_10min+i,k] > 0.0):
                         cloud_top_max[n1_10min+i]=hein[k]
                    if (BCu_frac_time[n1_10min+i,k] > 0.0):
                         cloudBCu_top_max[n1_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n1_10min+i,k] > 0.0):
                         cloud_base_max[n1_10min+i]=hein[k]
                    if (BCu_frac_time[n1_10min+i,k] > 0.0):
                         cloudBCu_base_max[n1_10min+i]=hein[k]


          for i in N.arange(n2_30min):
               modets_time_1h[n1_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n1_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n1_30min+i,k]=1000.*qv[i,:,:,k].mean()

     if (a+1 == 3):
          n_t2_10min=n1_10min+n2_10min
          n_t2_30min= n1_30min+n2_30min
          for i in N.arange(n3_10min):
               surf_precip_xy_time[n_t2_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n_t2_10min+i]=surf_precip[i,:,:].mean()

               shf_time[n_t2_10min+i]=shf[i]
               lhf_time[n_t2_10min+i]=lhf[i]
               modets_time[n_t2_10min+i]= mod_ts[i]/3600.

               BCu_frac_time[n_t2_10min+i,:]     = CondDiags[i,0,2,:]
               BCu_MF_time[n_t2_10min+i,:]     = CondDiags[i,1,2,:]
 
               AC_frac_time[n_t2_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n_t2_10min+i,:]     = CondDiags[i,1,4,:]
               ACzt_frac_time[:,n_t2_10min+i]     = AC_frac_time[n_t2_10min+i,:] 
               ACu_frac_time[n_t2_10min+i,:]     = CondDiags[i,0,5,:]
               ACuzt_frac_time[:,n_t2_10min+i]     = ACu_frac_time[n_t2_10min+i,:] 
               BCuzt_frac_time[:,n_t2_10min+i]     = BCu_frac_time[n_t2_10min+i,:] 
               ACu_MF_time[n_t2_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n_t2_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n_t2_10min+i,:]     = CondDiags[i,1,6,:]

               for k in N.arange(nz):
                    BCu_MF_time[n_t2_10min+i,k]= BCu_MF_time[n_t2_10min+i,k]*den_n_mean[k]
                    BCuzt_MF_time[k,n_t2_10min+i]     = BCu_MF_time[n_t2_10min+i,k] 

                    AC_MF_time[n_t2_10min+i,k]= AC_MF_time[n_t2_10min+i,k]*den_n_mean[k]
                    ACzt_MF_time[k,n_t2_10min+i]     = AC_MF_time[n_t2_10min+i,k] 
                    ACu_MF_time[n_t2_10min+i,k]= ACu_MF_time[n_t2_10min+i,k]*den_n_mean[k]
                    ACuzt_MF_time[k,n_t2_10min+i]     = ACu_MF_time[n_t2_10min+i,k] 
                    ACuS_MF_time[n_t2_10min+i,k]= ACuS_MF_time[n_t2_10min+i,k]*den_n_mean[k]
 
                    if (AC_frac_time[n_t2_10min+i,k] > 0.0):
                         cloud_top_max[n_t2_10min+i]=hein[k]
                    if (BCu_frac_time[n_t2_10min+i,k] > 0.0):
                         cloudBCu_top_max[n_t2_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n_t2_10min+i,k] > 0.0):
                         cloud_base_max[n_t2_10min+i]=hein[k]
                    if (BCu_frac_time[n_t2_10min+i,k] > 0.0):
                         cloudBCu_base_max[n_t2_10min+i]=hein[k]


          for i in N.arange(n3_30min):
               modets_time_1h[n_t2_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n_t2_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n_t2_30min+i,k]=1000.*qv[i,:,:,k].mean()

     if (a+1 == 4):
          n_t3_10min=n1_10min+n2_10min+n3_10min
          n_t3_30min=n1_30min+n2_30min+n3_30min 
          for i in N.arange(n4_10min):
               surf_precip_xy_time[n_t3_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n_t3_10min+i]=surf_precip[i,:,:].mean()

               shf_time[n_t3_10min+i]=shf[i]
               lhf_time[n_t3_10min+i]=lhf[i]
               modets_time[n_t3_10min+i]= mod_ts[i]/3600.

               BCu_frac_time[n_t3_10min+i,:]     = CondDiags[i,0,2,:]
               AC_frac_time[n_t3_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n_t3_10min+i,:]     = CondDiags[i,1,4,:]
               ACu_frac_time[n_t3_10min+i,:]     = CondDiags[i,0,5,:]
               ACuzt_frac_time[:,n_t3_10min+i]     = ACu_frac_time[n_t3_10min+i,:] 
               BCuzt_frac_time[:,n_t3_10min+i]     = BCu_frac_time[n_t3_10min+i,:] 
               ACu_MF_time[n_t3_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n_t3_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n_t3_10min+i,:]     = CondDiags[i,1,6,:]


               for k in N.arange(nz):
                    BCu_MF_time[n_t3_10min+i,k]= BCu_MF_time[n_t3_10min+i,k]*den_n_mean[k]
                    BCuzt_MF_time[k,n_t3_10min+i]     = BCu_MF_time[n_t3_10min+i,k] 

                    AC_MF_time[n_t3_10min+i,k]= AC_MF_time[n_t3_10min+i,k]*den_n_mean[k]
                    ACu_MF_time[n_t3_10min+i,k]= ACu_MF_time[n_t3_10min+i,k]*den_n_mean[k]
                    ACuzt_MF_time[k,n_t3_10min+i]     = ACu_MF_time[n_t3_10min+i,k]
                    ACuS_MF_time[n_t3_10min+i,k]= ACuS_MF_time[n_t3_10min+i,k]*den_n_mean[k]


                    if (AC_frac_time[n_t3_10min+i,k] > 0.0):
                         cloud_top_max[n_t3_10min+i]=hein[k]
                    if (BCu_frac_time[n_t3_10min+i,k] > 0.0):
                         cloudBCu_top_max[n_t3_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n_t3_10min+i,k] > 0.0):
                         cloud_base_max[n_t3_10min+i]=hein[k]
                    if (BCu_frac_time[n_t3_10min+i,k] > 0.0):
                         cloudBCu_base_max[n_t3_10min+i]=hein[k]


          for i in N.arange(n4_30min):
               modets_time_1h[n_t3_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n_t3_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n_t3_30min+i,k]=1000.*qv[i,:,:,k].mean()


     if (a+1 == 5):
          n_t4_10min=n1_10min+n2_10min+n3_10min+n4_10min
          n_t4_30min=n1_30min+n2_30min+n3_30min+n4_30min    

          for i in N.arange(n5_10min):
               surf_precip_xy_time[n_t4_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n_t4_10min+i]=surf_precip[i,:,:].mean()

               shf_time[n_t4_10min+i]=shf[i]
               lhf_time[n_t4_10min+i]=lhf[i]
               modets_time[n_t4_10min+i]= mod_ts[i]/3600.

               BCu_frac_time[n_t4_10min+i,:]     = CondDiags[i,0,2,:]
               BCu_MF_time[n_t4_10min+i,:]     = CondDiags[i,1,2,:]


               AC_frac_time[n_t4_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n_t4_10min+i,:]     = CondDiags[i,1,4,:]
               ACu_frac_time[n_t4_10min+i,:]     = CondDiags[i,0,5,:]
               ACu_MF_time[n_t4_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n_t4_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n_t4_10min+i,:]     = CondDiags[i,1,6,:]

               for k in N.arange(nz):
                    BCu_MF_time[n_t4_10min+i,k]= BCu_MF_time[n_t4_10min+i,k]*den_n_mean[k]

                    AC_MF_time[n_t4_10min+i,k]= AC_MF_time[n_t4_10min+i,k]*den_n_mean[k]
                    ACu_MF_time[n_t4_10min+i,k]= ACu_MF_time[n_t4_10min+i,k]*den_n_mean[k]
                    ACuS_MF_time[n_t4_10min+i,k]= ACuS_MF_time[n_t4_10min+i,k]*den_n_mean[k]

                    if (AC_frac_time[n_t4_10min+i,k] > 0.0):
                         cloud_top_max[n_t4_10min+i]=hein[k]
                    if (BCu_frac_time[n_t4_10min+i,k] > 0.0):
                         cloudBCu_top_max[n_t4_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n_t4_10min+i,k] > 0.0):
                         cloud_base_max[n_t4_10min+i]=hein[k]
                    if (BCu_frac_time[n_t4_10min+i,k] > 0.0):
                         cloudBCu_base_max[n_t4_10min+i]=hein[k]


          for i in N.arange(n5_30min):
               modets_time_1h[n_t4_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n_t4_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n_t4_30min+i,k]=1000.*qv[i,:,:,k].mean()

     if (a+1 == 6):
          n_t5_10min=n1_10min+n2_10min+n3_10min+n4_10min+n5_10min  
          n_t5_30min=n1_30min+n2_30min+n3_30min+n4_30min+n5_30min 
          for i in N.arange(n6_10min):
               surf_precip_xy_time[n_t5_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n_t5_10min+i]=surf_precip[i,:,:].mean()
 
               shf_time[n_t5_10min+i]=shf[i]
               lhf_time[n_t5_10min+i]=lhf[i]
               modets_time[n_t5_10min+i]= mod_ts[i]/3600.

               BCu_frac_time[n_t5_10min+i,:]     = CondDiags[i,0,2,:]
               BCu_MF_time[n_t5_10min+i,:]     = CondDiags[i,1,2,:]

               AC_frac_time[n_t5_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n_t5_10min+i,:]     = CondDiags[i,1,4,:]
               ACu_frac_time[n_t5_10min+i,:]     = CondDiags[i,0,5,:]
               ACu_MF_time[n_t5_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n_t5_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n_t5_10min+i,:]     = CondDiags[i,1,6,:]
               for k in N.arange(nz):
                    BCu_MF_time[n_t5_10min+i,k]= BCu_MF_time[n_t5_10min+i,k]*den_n_mean[k]
                    AC_MF_time[n_t5_10min+i,k]= AC_MF_time[n_t5_10min+i,k]*den_n_mean[k]
                    ACu_MF_time[n_t5_10min+i,k]= ACu_MF_time[n_t5_10min+i,k]*den_n_mean[k]
                    ACuS_MF_time[n_t5_10min+i,k]= ACuS_MF_time[n_t5_10min+i,k]*den_n_mean[k]

                    if (AC_frac_time[n_t5_10min+i,k] > 0.0):
                         cloud_top_max[n_t5_10min+i]=hein[k]
                    if (BCu_frac_time[n_t5_10min+i,k] > 0.0):
                         cloudBCu_top_max[n_t5_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n_t5_10min+i,k] > 0.0):
                         cloud_base_max[n_t5_10min+i]=hein[k]
                    if (BCu_frac_time[n_t5_10min+i,k] > 0.0):
                         cloudBCu_base_max[n_t5_10min+i]=hein[k]

          for i in N.arange(n6_30min):
               modets_time_1h[n_t5_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n_t5_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n_t5_30min+i,k]=1000.*qv[i,:,:,k].mean()


     if (a+1 == 7):
          n_t6_10min=n1_10min+n2_10min+n3_10min+n4_10min+n5_10min+n6_10min
          n_t6_30min=n1_30min+n2_30min+n3_30min+n4_30min+n5_30min+n6_30min
          for i in N.arange(n7_10min):
               surf_precip_xy_time[n_t6_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n_t6_10min+i]=surf_precip[i,:,:].mean()
               shf_time[n_t6_10min+i]=shf[i]
               lhf_time[n_t6_10min+i]=lhf[i]
               modets_time[n_t6_10min+i]= mod_ts[i]/3600.

               BCu_frac_time[n_t6_10min+i,:]     = CondDiags[i,0,2,:]
               BCu_MF_time[n_t6_10min+i,:]     = CondDiags[i,1,2,:]

               AC_frac_time[n_t6_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n_t6_10min+i,:]     = CondDiags[i,1,4,:]
               ACu_frac_time[n_t6_10min+i,:]     = CondDiags[i,0,5,:]
               ACu_MF_time[n_t6_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n_t6_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n_t6_10min+i,:]     = CondDiags[i,1,6,:]
               for k in N.arange(nz):
                    BCu_MF_time[n_t6_10min+i,k]= BCu_MF_time[n_t6_10min+i,k]*den_n_mean[k]
                    AC_MF_time[n_t6_10min+i,k]= AC_MF_time[n_t6_10min+i,k]*den_n_mean[k]
                    ACu_MF_time[n_t6_10min+i,k]= ACu_MF_time[n_t6_10min+i,k]*den_n_mean[k]

                    if (AC_frac_time[n_t6_10min+i,k] > 0.0):
                         cloud_top_max[n_t6_10min+i]=hein[k]
                    if (BCu_frac_time[n_t6_10min+i,k] > 0.0):
                         cloudBCu_top_max[n_t6_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n_t6_10min+i,k] > 0.0):
                         cloud_base_max[n_t6_10min+i]=hein[k]
                    if (BCu_frac_time[n_t6_10min+i,k] > 0.0):
                         cloudBCu_base_max[n_t6_10min+i]=hein[k]

          for i in N.arange(n7_30min):
               modets_time_1h[n_t6_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n_t6_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n_t6_30min+i,k]=1000.*qv[i,:,:,k].mean()


     if (a+1 == 8):
          n_t7_10min=n1_10min+n2_10min+n3_10min+n4_10min+n5_10min+n6_10min+n7_10min
          n_t7_30min=n1_30min+n2_30min+n3_30min+n4_30min+n5_30min+n6_30min+n7_30min
          for i in N.arange(n8_10min):
               surf_precip_xy_time[n_t7_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n_t7_10min+i]=surf_precip[i,:,:].mean()

               shf_time[n_t7_10min+i]=shf[i]
               lhf_time[n_t7_10min+i]=lhf[i]
               modets_time[n_t7_10min+i]= mod_ts[i]/3600.

               BCu_frac_time[n_t7_10min+i,:]     = CondDiags[i,0,2,:]
               BCu_MF_time[n_t7_10min+i,:]     = CondDiags[i,1,2,:]

               AC_frac_time[n_t7_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n_t7_10min+i,:]     = CondDiags[i,1,4,:]
               ACu_frac_time[n_t7_10min+i,:]     = CondDiags[i,0,5,:]
               ACu_MF_time[n_t7_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n_t7_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n_t7_10min+i,:]     = CondDiags[i,1,6,:]

               for k in N.arange(nz):
                    BCu_MF_time[n_t7_10min+i,k]= BCu_MF_time[n_t7_10min+i,k]*den_n_mean[k]
                    AC_MF_time[n_t7_10min+i,k]= AC_MF_time[n_t7_10min+i,k]*den_n_mean[k]
                    ACu_MF_time[n_t7_10min+i,k]= ACu_MF_time[n_t7_10min+i,k]*den_n_mean[k]
                    ACuS_MF_time[n_t7_10min+i,k]= ACuS_MF_time[n_t7_10min+i,k]*den_n_mean[k]
                    if (AC_frac_time[n_t7_10min+i,k] > 0.0):
                         cloud_top_max[n_t7_10min+i]=hein[k]
                    if (BCu_frac_time[n_t7_10min+i,k] > 0.0):
                         cloudBCu_top_max[n_t7_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n_t7_10min+i,k] > 0.0):
                         cloud_base_max[n_t7_10min+i]=hein[k]
                    if (BCu_frac_time[n_t7_10min+i,k] > 0.0):
                         cloudBCu_base_max[n_t7_10min+i]=hein[k]

          for i in N.arange(n8_30min):
               modets_time_1h[n_t7_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n_t7_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n_t7_30min+i,k]=1000.*qv[i,:,:,k].mean()

     if (a+1 == 9):
          n_t8_10min=n1_10min+n2_10min+n3_10min+n4_10min+n5_10min+n6_10min+n7_10min+n8_10min
          n_t8_30min=n1_30min+n2_30min+n3_30min+n4_30min+n5_30min+n6_30min+n7_30min+n8_30min
          for i in N.arange(n9_10min):
               surf_precip_xy_time[n_t8_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n_t8_10min+i]=surf_precip[i,:,:].mean()

               shf_time[n_t8_10min+i]=shf[i]
               lhf_time[n_t8_10min+i]=lhf[i]
               modets_time[n_t8_10min+i]= mod_ts[i]/3600.


               BCu_frac_time[n_t8_10min+i,:]     = CondDiags[i,0,2,:]
               BCu_MF_time[n_t8_10min+i,:]     = CondDiags[i,1,2,:]
               AC_frac_time[n_t8_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n_t8_10min+i,:]     = CondDiags[i,1,4,:]
               ACu_frac_time[n_t8_10min+i,:]     = CondDiags[i,0,5,:]
               ACu_MF_time[n_t8_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n_t8_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n_t8_10min+i,:]     = CondDiags[i,1,6,:]

               for k in N.arange(nz):
                    BCu_MF_time[n_t8_10min+i,k]= BCu_MF_time[n_t8_10min+i,k]*den_n_mean[k]
                    AC_MF_time[n_t8_10min+i,k]= AC_MF_time[n_t8_10min+i,k]*den_n_mean[k]
                    ACu_MF_time[n_t8_10min+i,k]= ACu_MF_time[n_t8_10min+i,k]*den_n_mean[k]
                    ACuS_MF_time[n_t8_10min+i,k]= ACuS_MF_time[n_t8_10min+i,k]*den_n_mean[k]

                    if (AC_frac_time[n_t8_10min+i,k] > 0.0):
                         cloud_top_max[n_t8_10min+i]=hein[k]
                    if (BCu_frac_time[n_t8_10min+i,k] > 0.0):
                         cloudBCu_top_max[n_t8_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n_t8_10min+i,k] > 0.0):
                         cloud_base_max[n_t8_10min+i]=hein[k]
                    if (BCu_frac_time[n_t8_10min+i,k] > 0.0):
                         cloudBCu_base_max[n_t8_10min+i]=hein[k]

          for i in N.arange(n9_30min):
               modets_time_1h[n_t8_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n_t8_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n_t8_30min+i,k]=1000.*qv[i,:,:,k].mean()
 
     if (a+1 == 10):
          n_t9_10min=n1_10min+n2_10min+n3_10min+n4_10min+n5_10min+n6_10min+n7_10min+n8_10min+n9_10min
          n_t9_30min=n1_30min+n2_30min+n3_30min+n4_30min+n5_30min+n6_30min+n7_30min+n8_30min+n9_30min
          for i in N.arange(n10_10min):
               surf_precip_xy_time[n_t9_10min+i,:,:]=surf_precip[i,:,:]
               surf_precip_time[n_t9_10min+i]=surf_precip[i,:,:].mean()

               shf_time[n_t9_10min+i]=shf[i]
               lhf_time[n_t9_10min+i]=lhf[i]
               modets_time[n_t9_10min+i]= mod_ts[i]/3600.

               BCu_frac_time[n_t9_10min+i,:]     = CondDiags[i,0,2,:]
               BCu_MF_time[n_t9_10min+i,:]     = CondDiags[i,1,2,:]

               AC_frac_time[n_t9_10min+i,:]     = CondDiags[i,0,4,:]
               AC_MF_time[n_t9_10min+i,:]     = CondDiags[i,1,4,:]
               ACu_frac_time[n_t9_10min+i,:]     = CondDiags[i,0,5,:]
               ACu_MF_time[n_t9_10min+i,:]     = CondDiags[i,1,5,:]
               ACuS_frac_time[n_t9_10min+i,:]     = CondDiags[i,0,6,:]
               ACuS_MF_time[n_t9_10min+i,:]     = CondDiags[i,1,6,:]

               for k in N.arange(nz):
                    BCu_MF_time[n_t9_10min+i,k]= BCu_MF_time[n_t9_10min+i,k]*den_n_mean[k]
                    AC_MF_time[n_t9_10min+i,k]= AC_MF_time[n_t9_10min+i,k]*den_n_mean[k]
                    ACu_MF_time[n_t9_10min+i,k]= ACu_MF_time[n_t9_10min+i,k]*den_n_mean[k]
                    ACuS_MF_time[n_t9_10min+i,k]= ACuS_MF_time[n_t9_10min+i,k]*den_n_mean[k]

                    if (AC_frac_time[n_t9_10min+i,k] > 0.0):
                         cloud_top_max[n_t9_10min+i]=hein[k]
                    if (BCu_frac_time[n_t9_10min+i,k] > 0.0):
                         cloudBCu_top_max[n_t9_10min+i]=hein[k]
                    k=nz-k-1
                    if (AC_frac_time[n_t9_10min+i,k] > 0.0):
                         cloud_base_max[n_t9_10min+i]=hein[k]
                    if (BCu_frac_time[n_t9_10min+i,k] > 0.0):
                         cloudBCu_base_max[n_t9_10min+i]=hein[k]


          for i in N.arange(n10_30min):
               modets_time_1h[n_t9_30min+i]= mod_ts_1h[i]/3600.
               for k in N.arange(nz):
                    th_mean[n_t9_30min+i,k]=th[i,:,:,k].mean()+th_ref_prof[k]
                    qv_mean[n_t9_30min+i,k]=1000.*qv[i,:,:,k].mean()


surf_flx[:]=shf_time[:]+lhf_time[:]

