from codes.utils.sampling import sampling
from scipy.optimize import minimize, differential_evolution
from codes.utils.keyrate import keyrate
import numpy as np

class MyBounds:
    def __init__(self, xmax=[1-1e-14, 1-1e-14, 1-1e-14, 1-1e-14, 1-1e-14], xmin=[1e-14, 1e-14, 1e-14, 1e-14, 1e-14]):
        self.xmax = np.array(xmax)
        self.xmin = np.array(xmin)

    def __call__(self, x_new):
        x = x_new
        tmax = bool(np.all(x <= self.xmax))
        tmin = bool(np.all(x >= self.xmin))
        return tmax and tmin
    
def optimize(kwargs, l=0,x=[],maxit = 100):
    if x != []:
        pass
    else:
        x = sampling()

    def fun(*x):
        x = x[0]
        return -keyrate(x, l, kwargs)
    bnds = ((1e-14, 10-1e-14), (1e-5, 10-1e-14), (1e-14, 1-1e-14),
            (1e-14, 1-1e-14), (1e-14, 1-1e-14))
    opt = {'xatol': 1e-10, 'fatol': 1e-9,
           'adaptive': True, 'maxiter': 1000, 'maxfev': 5000}
    OptRes = minimize(fun, x, method='Nelder-Mead',
                      bounds=bnds, tol=None, callback=None, options=opt)
    # attributes of OptRes: class OptimizeResult, see: https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.OptimizeResult.html#scipy.optimize.OptimizeResult
    if OptRes.success:
        return OptRes.x, OptRes.success
    elif 'x' in OptRes:
        return OptRes.x, True
    else:
        return [], False


def global_optimize(kwargs,l=0, x=[], global_method='basinghopping',maxit = 100):
    def fun(*x):
        x = x[0]
        r = keyrate(x, l, kwargs)
        # print(r)
        return -r
    bnds = ((1e-14, 1-1e-14), (1e-14, 1-1e-14), (1e-14, 1-1e-14),
            (1e-14, 1-1e-14), (1e-14, 1-1e-14))
    mybounds = MyBounds()
    if x != []:
        pass
    else:
        x = sampling()

    bnds = ((1e-14, 1-1e-14), (1e-14, 1-1e-14), (1e-14, 1-1e-14),
            (1e-14, 1-1e-14), (1e-14, 1-1e-14))
    opt = {'xatol': 1e-10, 'fatol': 1e-9,
           'adaptive': True, 'maxiter': maxit, 'maxfev': maxit}
    minimizer_kwargs = {"method": 'Nelder-Mead', 'bounds': bnds,
                        'tol': None, 'callback': None, 'options': opt}

    def print_fun(x, f, accepted):
        print('at maximum %f accepted %d' % (-f, int(accepted)), x)
        ''''
    OptRes = basinhopping(fun,x,minimizer_kwargs=minimizer_kwargs,niter = 50,callback =None, accept_test=mybounds)
    '''
    '''
    OptRes1 = brute(fun,bnds,finish=optimize.fmin)
    if self._keyrate(OptRes1['x'])>self._keyrate(OptRes['x']):
        OptRes = OptRes1
    '''
    OptRes = differential_evolution(fun, bnds)
    '''
    if self._keyrate(OptRes1['x'])>self._keyrate(OptRes['x']):
        OptRes = OptRes1
        '''

    if OptRes['success']:
        return OptRes.x, OptRes.success
    elif 'x' in OptRes:
        return OptRes.x, True
    else:
        return [], False
