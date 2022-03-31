from cmath import tau
import autograd.numpy as np
import numdifftools as nd
import scipy.io as mat
from minimization import Minimizer


path = 'BTR/FPAnalysis/simulation_results'
inputs = mat.loadmat(path + '/inputs1.mat')

session = 4

def nonlin(x):
    return np.maximum(x, np.zeros(np.size(x)))

def build_qfun(inputs, session):
    
    tau = inputs['tau']
    alpha = inputs['alpha']
    W = inputs['W'][:,:,session-1]
    V_ff = inputs['V_ff'][:,session-1]
    rate = inputs['rate'][:,session-1]

    def qfun(x):
        return 1 / 512 * np.sum((- x + V_ff + W * alpha * nonlin(x) / tau) ** 2)

    return qfun

V_rec = inputs['V_rec'][:,session-1]

qfun = build_qfun(inputs, session)


minimizer = Minimizer()
fp = minimizer.adam_optimization(qfun, V_rec)

print(fp)
