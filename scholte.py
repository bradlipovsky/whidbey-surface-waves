from numpy.random import random_sample
from time import perf_counter
from tqdm import tqdm
import pickle
import pandas as pd
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress

def pysurf96(hr,vp,vs,rho,model_file='model-from-python.d'):

    lines = ['MODEL',
             'TEST MODEL',
            'ISOTROPIC',
            'KGS',
            'FLAT EARTH',
            '1-D',
            'CONSTANT VELOCITY',
            'LINE08',
            'LINE09',
            'LINE10',
            'LINE11',
            'HR    VP    VS    RHO   QP  QS  ETAP ETAS FREFP FREFS']
    for h,p,s,r in zip(hr,vp,vs,rho):
#         lines.append(f'{h:.3f} {p:.3f} {s:.3f} {r:.3f}  0.0 0.0 0.0  0.0  1.0   1.0')
        lines.append(f'{h} {p} {s} {r}  0.0 0.0 0.0  0.0  1.0   1.0')
    
    with open(model_file, 'w') as f:
        for line in lines:
            f.write(line)
            f.write('\n')
    
    result = subprocess.run(['./model.sh'], stdout=subprocess.PIPE)
    
    data=pd.read_fwf('SDISPR.TXT',header=1,colspecs=[(22,34),(38,52)])

    '''
    It's a bit messy, but the following parses the 0th and 1st scholte mode from the output file.
    '''
    N0=12
    N1=6

    f0 = np.array(data['FREQ'][0:N0].astype(float))
    c0 = np.array(data['PHASE VEL'][0:N0].astype(float)*1e3)
    k0 = f0/c0

    f1 = np.array(data['FREQ'][N0+2:N0+2+N1].astype(float))
    c1 = np.array(data['PHASE VEL'][N0+2:N0+2+N1].astype(float)*1e3)
    k1 = f1/c1

    f = [f0,f1]
    k = [k0,k1]
    return f,k

def wrap_four_layer(p,h,k_obs,
                        vp = [1.5,1.6,1.6,1.6,1.6], 
                        rho= [1,1.5,1.5,1.5,1.5]  ):

    this_vs = [0,p[0],p[1],p[2],p[3]]
    f,k = pysurf96(h,vp,this_vs,rho)
    misfit = np.abs(k[3:]-k_obs).sum()

    return misfit

def do_mcmc_4_param():
    t0 = perf_counter()

    step_size = 2.5e-4 # 0.25 m/s
    # k_sensitivity = 0.000418 * len(kobs)
    k_sensitivity = 7.5e-5
    x0 = [0.038,0.043,0.080,0.60]
    h0 = [0.1,0.005,0.010,0.020,0.0]
    misfit0 = wrap_four_layer(x0,h0,kobs)
    N = 100000

    x_list = []
    for i in tqdm(range(N)):
        x = x0 + step_size*(2*random_sample(len(x0))-1)
        misfit = wrap_four_layer(x,h0,kobs)
        L =  np.exp(-(misfit-misfit0)**2 / k_sensitivity**2)
        thr = random_sample()
        if L >= thr:
            x_list.append(x)
            x0 = x
            m0 = m
    x_list = np.array(x_list) * 1e3
    print(f'Runtime was {perf_counter()-t0}')
    
    with open('100k_4param_MCMC.pickle', 'wb') as handle:
        pickle.dump(x_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
def do_mcmc_8param(kobs, 
                   N = 200000, 
                   filename='100k_4param_MCMC.pickle',
                   fixhalfspace=True,
                   restartfile=None):
    
    t0 = perf_counter()

    step_size = 1.0e-4 # 0.1 m/s
    step_size_h = 0.01e-3 # 0.01 m
    # k_sensitivity = 0.000418 * len(kobs)
    k_sensitivity = 5e-5
    x0 = [0.038,0.043,0.080,0.60]
    h0 = [0.1,0.005,0.010,0.020,0.0]
    m0 = wrap_four_layer(x0,h0,kobs)

    x_list = []
    for i in tqdm(range(N)):
        x = x0 + step_size*(2*random_sample(len(x0))-1)
        if fixhalfspace:
            x[3]=x0[3]
        h = h0 + step_size_h*(2*random_sample(len(h0))-1)
        h[-1]=0
        m = wrap_four_layer(x,h,kobs)
        L =  np.exp(-(m-m0)**2 / k_sensitivity**2)
        thr = random_sample()
        if L >= thr:
            x_list.append( np.hstack((h0,x0)) )
            x0 = x
            h0 = h
            m0 = m
    x_list = np.array(x_list) * 1e3
    print(f'Runtime was {perf_counter()-t0}')
    
    if restartfile is not None:
        x_list_old = pickle.load(open(restartfile, 'rb'))
        x_list = np.vstack((x_list_old,x_list))
    
    with open(filename, 'wb') as handle:
        pickle.dump(x_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
def exp_vs(z,z0,exp,coeff):
    return (coeff*(z - z0))**exp


def plt_dispersion(f,k,fobs,kobs,z,vs,save=False):

    fig,ax=plt.subplots(1,2,figsize=(16,9))
    fig.patch.set_facecolor('w')
    plt.subplot(121)
    plt.step(vs,z,label='Vs (km/s)')
    # plt.step(p,z,label='Vp (km/s)')
    # plt.step(rho,z,'--',label='Density (kg/m**3)')
    plt.ylim([max(z),0])
    plt.grid()
    # plt.legend(fontsize=18)
    plt.ylabel('Depth (m)',fontsize=22)
    plt.xlabel('Shear Wave Speed (m/s)',fontsize=22)
    plt.text(0.1,50,'water',fontsize=18)
    plt.text(100,120,'fluidized \nsediments',fontsize=18)
    plt.text(400,175,'   compact \n sediments',fontsize=18)

    plt.xlim([0, 650])

    plt.subplot(122)
    ind = 0
    for kk,ff in zip(k,f):
        plt.plot(kk,ff,'-o',label=f'Model Mode {ind}',markersize=12)
        ind = ind+1
    ind = 0
    for kk,ff in zip(kobs,fobs):
        plt.plot( kk , ff, 'o' , label=f'Obs. Mode {ind}',markersize=10)
        ind = ind+1
    plt.ylabel('Frequency (Hz)',fontsize=22)
    plt.xlabel(r'Ordinary Wavenumber, $1/\lambda$ (1/m)',fontsize=22)

    
#     fit = linregress(kobs,fobs)
#     plt.plot(kobs, fit.slope*kobs + fit.intercept,'--',label=f'Constant vg ({fit.slope:.0f} m/s)')
#     print(fit.slope)
#     print(fit.intercept)
    plt.legend(fontsize=18)
    plt.grid()
    if save:
        plt.savefig('model.eps')
    else:
        plt.show()