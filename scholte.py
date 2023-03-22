from numpy.random import random_sample
from time import perf_counter
from tqdm import tqdm
import pickle
import pandas as pd
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
from scipy.interpolate import interp1d

def pysurf96(hr,vp,vs,rho,model_file='model-from-python.d',N0 = 18, N1 = 15):

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

    f = []
    f_this_mode = []
    k = []
    k_this_mode = []
    for this_c,this_f in zip(data['PHASE VEL'],data['FREQ']):

        if this_f == 'EIGH WAVE':
            continue
        elif this_f == 'FREQ':
            f.append( np.array(f_this_mode) )
            f_this_mode = []
            k.append( np.array(k_this_mode) )
            k_this_mode = []
        else:
            f_this_mode.append(float(this_f))
            k_this_mode.append(float(this_f)/float(this_c)/1000)
    f.append( np.array(f_this_mode) )
    k.append( np.array(k_this_mode) )
    return f,k

def function_wrapper(vs,h,k_obs, f_obs,
                        vp = [1.5,1.6,1.6,1.6,1.6], 
                        rho= [1,1.5,1.5,1.5,1.5]  ):

    f,k = pysurf96(h,vp,vs,rho)
    
    misfit = 0
    n = 0
    for fm_mode,km_mode,fo_mode,ko_mode in zip(f,k,f_obs,k_obs):
        # Loop over the different Scholte wave modes
        if len(km_mode) < 2:
            # surf96 failed? Don't use this mode.
            continue
        else:
            n = n + len(km_mode)
        f_interp = interp1d(km_mode,fm_mode)
        # index to account for the fact that surf96 may not have found
        # the dispersion curve where we have observations -- only analyze
        # where the model and observations overlap.
        index = (min(km_mode) < ko_mode) & (ko_mode < max(km_mode))
        f_model = f_interp(ko_mode[index])
        misfit = misfit + np.abs(f_model - fo_mode[index]).sum()
    misfit = misfit / n # divide to arrive at the mean absolute error

    return misfit
    
def do_mcmc( kobs, fobs, h, x0,
                    N = 200000, 
                    filename='test.pickle',
                    fixhalfspace=True,
                    fixdepths=True,
                    fix_water_layer=True,
                    restartfile=None,
                    step_size = 1.0e-4, # 0.1 m/s
                    step_size_h = 0.01e-3, # 0.01 m
                    f_sensitivity = 1e-3,
                    verbose=False
                  ):
    
    t0 = perf_counter()
    
    misfit0 = function_wrapper(x0,h,kobs,fobs)
    L0 =  -misfit0**2 / f_sensitivity**2
    
    if verbose: print(f'Probability of initial model = {np.exp(L0)}')
    if verbose: print(f'misfit0 = {misfit0}')
        
    x_list = []
    for i in tqdm(range(N)):
        x = x0 + np.random.default_rng().normal(0, step_size, len(x0))
        if fix_water_layer: x[0]=0 # zero shear wave speed in a water layer
            
        misfit = function_wrapper(x,h,kobs,fobs)
        Lp =  -misfit**2 / f_sensitivity**2
        alpha = np.exp(Lp)/np.exp(L0)
#         alpha = np.exp(Lp-L0)
        u = random_sample()
        if verbose: print(f'alpha={alpha},  u={u}')
        if  alpha > u: # accept the proposition
            
            x_list.append( x )
            x0 = x
            misfit0 = misfit
            
            
    x_list = np.array(x_list) * 1e3
    print(f'Runtime was {perf_counter()-t0}')
    print(f'Acceptances = {len(x_list)}, Acceptance Ratio = {len(x_list)/N}')
    
    if restartfile is not None:
        x_list_old = pickle.load(open(restartfile, 'rb'))
        x_list = np.vstack((x_list_old,x_list))
    
    with open(filename, 'wb') as handle:
        pickle.dump(x_list, handle, protocol=pickle.HIGHEST_PROTOCOL)

        

        
def load_dispersion_observations():
    kobs0_manual = np.array([0.04, 0.035, 0.030, 0.025, 0.020, 0.0150, 0.0113, 0.0059, 0.0029, 0.00042])
    fobs0_manual = np.array([1.43, 1.35, 1.250, 1.140, 1.02, 0.88, 0.77, 0.585, 0.46, 0.378])

    kobs1_manual = np.array([0.001, 0.002, 0.003, 0.004, 0.005, 0.007, 0.009, 0.011, 0.013])
    fobs1_manual = np.array([0.756,   0.80, 0.835, 0.865, 0.90, 0.95, 1.02, 1.09, 1.15])

    kobs2_manual = np.array([0.002,0.003,0.004,0.006,0.008,0.010,0.012,0.014])
    fobs2_manual = np.array([1.10,1.21,1.26,1.35,1.45,1.55,1.65,1.75])

    kobs = [kobs0_manual,kobs1_manual,kobs2_manual]
    fobs = [fobs0_manual,fobs1_manual,fobs2_manual]
    
    return kobs,fobs


