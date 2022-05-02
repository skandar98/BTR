from cmath import tau
import autograd.numpy as np
import numpy as npy
import random as random
import scipy.io as mat
from autograd import make_jvp
from minimization import Minimizer
from numpy.linalg import eig
from sklearn.decomposition import PCA
from mayavi import mlab
from autograd import jacobian
from utils import *

path = 'FPAnalysis/simulation_results'
inputs = mat.loadmat(path + '/inputs_small_W_1rep.mat')
id = 'small_w_rep1_jac'

experiments = 3 #number of eperimnets simulated
n_trials = 1 #number of sampled trials (simulations). Max=200
n_tpoints = 1 #number of time points sampled from the same simualtion. Max=100


for experiment in range(experiments):
    # sampling states and activations
    print('Experiment:'+str(experiment+1)+'\n Sampling....')
    V_ff_sample, W_sample, V_rec_sample = sample_states(inputs, n_trials, n_tpoints, 1, 0.5)
    #print(npy.shape(V_ff), npy.shape(W), npy.shape(V_rec))
    # Find fixed points through minimizing the velocity function
    
    fps_good = npy.empty((512,1))
    jac_good  = npy.empty((512,512,1))
    good = 0
    for i in range(n_trials*n_tpoints):
        fp = npy.empty((512,1))

        print('trial:', str(int(n_tpoints/(i+1))), 'time point:', str((int(n_trials/(i+1)))))
        minimizer = Minimizer(alr_decayr=0.5, epsilon=0.1, init_agnc=1)
        fun_ds = build_ds(inputs, V_ff_sample[:,i], W_sample[:,:,i])
        fun = build_fun(inputs, V_ff_sample[:,i], W_sample[:,:,i])
        qfun = build_qfun(inputs, V_ff_sample[:,i], W_sample[:,:,i])
        fp[:,0] = minimizer.adam_optimization(fun_ds, V_rec_sample[:, i])
        q = qfun(fp[:,0])
        jac = make_jvp(fun)(fp[:,0])
        #jac =  jacobian(fun)(fp)
        if q < 1:
            if good == 0:
                fps_good[:,0] = fp[:,0]
                jac_good = jac(np.eye(512)[:,0])
            else:
                fps_good = npy.append(fps_good, fp, 1)
                #jac_good = npy.append(jac_good, jac, 2)
            good = good+1
    #fps_good = npy.delete(fps_good, 0, 0)
    #jac_good = npy.delete(jac_good, [0, 0], 0)
    if good == 0:
        print('No fixed points were found...')
    else:
        ux, idx = npy.unique(fps_good.round(decimals=5), axis=1, return_index=True)
        # Select and save unique fixed points
        unique_fps = fps_good[:,idx]
        unique_V_ff= V_ff_sample[:,idx]
        unique_W= W_sample[:,:,idx]
        unique_V_rec= V_rec_sample[:,idx]
        #unique_jac = jac
        #compute jacobian at fixed point
        #jac_fun= make_jvp(fun)(unique_fps)
        #value_of_fun, jac = jac_fun(unique_fps)
    

        mdict = {'fps': unique_fps, 'V_ff': unique_V_ff, 'V_rec': unique_V_rec, 'W': unique_W}
        mat.savemat(path+'/'+str(id)+'_exp'+str(experiment)+'.mat', mdict)


