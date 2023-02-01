from numpy.random import random_sample
from time import perf_counter
from tqdm import tqdm
import pickle
import pandas as pd
import subprocess
import numpy as np


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
    f = data['FREQ']
    k = data['FREQ']/data['PHASE VEL']/1e3
    return f,k

def wrap_four_layer(p,h,k_obs,vp = [1.5,1.6,1.6,1.6,1.6], rho= [1,1.5,1.5,1.5,1.5]):
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
    
def do_mcmc_8param(kobs, N = 200000, filename='100k_4param_MCMC.pickle',fixhalfspace=True):
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
    
    with open(filename, 'wb') as handle:
        pickle.dump(x_list, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
def exp_vs(z,z0,exp,coeff):
    return (coeff*(z - z0))**exp