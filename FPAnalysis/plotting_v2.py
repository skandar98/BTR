from cmath import tau
from turtle import width
import autograd.numpy as np
import numpy as npy
import scipy.io as mat
from numpy.linalg import eig
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import mayavi.mlab as mlab
from sklearn.preprocessing import scale
from utils import *

path = 'FPAnalysis/simulation_results'
experiments = 3
id = 11

Data = mat.loadmat(path + '/'+str(id)+'_exp0.mat')
n_fps = len(Data['V_rec'][1,:])
X_rec = npy.zeros((n_fps, 3, 3))
X_fps = npy.zeros((n_fps, 3, 3))

for exp in range(experiments):
     if exp==0:
          color = 'r'
     elif exp==1:
         color = 'g'
         Data = mat.loadmat(path + '/'+str(id)+'_exp'+str(exp)+'.mat')
     elif exp==2:
         color = 'm'
         Data = mat.loadmat(path + '/'+str(id)+'_exp'+str(exp)+'.mat')
     

     fps = npy.transpose(Data['fps'])
     V_ff = npy.transpose(Data['V_ff'])
     V_rec = npy.transpose(Data['V_rec'])
     W = Data['W']
     #jac = Data['jac']

     pca = PCA(n_components=3)
     pca.fit(fps)

     #for i in range(len(fps[0])):
          #eigvals, eigvec = eig(jac[:,i])
          #idx = npy.argwhere(npy.real(eigvals) > 0)
          #countgreaterzero = npy.sum(eigvals > 0)
          #if countgreaterzero == 0:
          #     print('stable fixed point was found.')
            
          #elif countgreaterzero > 0:
          #     print('saddle point was found.')
            
            
     X_rec[:,:,exp] = pca.transform(V_rec)
     X_fps[:,:,exp]  = pca.transform(fps)

     #mlab.plot3d(X_rec[:,0,0], X_rec[:,1,0], X_rec[:,2,0], color=(0.3, 1, 1), tube_radius=0.1)
     #mlab.points3d(X_fps[:,0,0], X_fps[:,1,0], X_fps[:,2,0], color=(0.3, 1, 1), scale_factor=.5)

    #plt.plot(X_rec[0], X_rec[1], color = color)
    #plt.scatter(X_fps[0], X_fps[1], color = color)
#s = npy.linspace(0.2, 0.25, num=n_fps)
mlab.plot3d(X_rec[:,0,0], X_rec[:,1,0], X_rec[:,2,0], color=(0.3, 1, 1), tube_radius=0.01)
mlab.points3d(X_fps[:,0,0], X_fps[:,1,0], X_fps[:,2,0], color=(1, 0.3, 1), scale_factor=0.25)
#mlab.plot3d(X_rec[:,0,1], X_rec[:,1,1], X_rec[:,2,1], color=(1, 0.3, 1), tube_radius=0.01)   
#mlab.points3d(X_fps[:,0,1], X_fps[:,1,1], X_fps[:,2,1], color=(1, 0.3, 1), scale_factor=.25) 
#mlab.plot3d(X_rec[:,0,2], X_rec[:,1,2], X_rec[:,2,2], color=(1, 1, 0.3), tube_radius=0.025)
#mlab.points3d(X_fps[:,0,2], X_fps[:,1,2], X_fps[:,2,2], color=(1, 1, 0.3), scale_factor=.1)
mlab.show()

