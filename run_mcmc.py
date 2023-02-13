from scholte import *
kobs = np.array([0.0343, 0.0305, 0.0246, 0.020, 0.0150, 0.0113, 0.0059, 0.0029, 0.00042])
for i in range(100):
    print(f'Restart {i+1}')
    do_mcmc_8param(kobs, N = 20000, filename='test.pickle',restartfile='test.pickle')
