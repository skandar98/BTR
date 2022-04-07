from cmath import tau
import autograd.numpy as np
import numdifftools as nd
from numpy import dtype
import scipy.io as mat
from autograd import make_jvp
from minimization import Minimizer
from numpy.linalg import eig
from sklearn.decomposition import PCA

path = 'FPAnalysis/simulation_results'
inputs = mat.loadmat(path + '/inputs1.mat')

session = 8

def nonlin(x):
    return np.maximum(x, np.zeros(np.size(x)))


def build_fun(inputs, session):
    tau = inputs['tau']
    alpha = inputs['alpha']
    W = inputs['W'][:,:,session-1]
    V_ff = inputs['V_ff'][:,session-1]


    def fun(x):
        return (- x + V_ff + W * alpha * nonlin(x) / tau)
    return fun

def build_qfun(fun):
    
    def qfun(x):
        return 1 / 512 * np.sum(fun(x) ** 2)

    return qfun

def sample_activations(V_rec):
    samples = V_rec
    return samples 

def add_gaussian_noise(self, data, noise_scale=0.0):
        """ Adds IID Gaussian noise to Numpy data.
        Args:
            data: Numpy array.
            noise_scale: (Optional) non-negative scalar indicating the
            standard deviation of the Gaussian noise samples to be generated.
            Default: 0.0.
        Returns:
            Numpy array with shape matching that of data.
        Raises:
            ValueError if noise_scale is negative.
        """
        # Add IID Gaussian noise
        if noise_scale == 0.0:
            return data # no noise to add
        if noise_scale > 0.0:
            return data + noise_scale * self.rng.randn(*data.shape)
        elif noise_scale < 0.0:
            raise ValueError('noise_scale must be non-negative,'
                             ' but was %f' % noise_scale)


V_ff = inputs['V_ff'][:,session-1]
V_rec = inputs['V_rec'][:,session-1]

fun = build_fun(inputs,session)
qfun = build_qfun(fun)



minimizer = Minimizer()
fp = minimizer.adam_optimization(qfun, V_rec)

#compute jacobian at fixed point
jac_fun= make_jvp(fun)(fp)
value_of_fun, jac = jac_fun(fp)



eigenvalues,v=eig(jac)
print (eigenvalues)

X = np.stack((V_ff, V_rec), axis =1)
pca = PCA(n_components=2)
pca.fit(X)
print(pca.explained_variance_ratio_)
print(pca.singular_values_)

mdict = {'fp': fp, 'jacobian': jac}
mat.savemat(path+'/begin_fp.mat',mdict)
