from codes.utils.optimize import global_optimize, optimize
from codes.utils.sampling import sampling
from codes.utils.keyrate import keyrate
def search_for_initial_point(kwargs):
    best_kr = 0
    best_x = []
    for i in range(1):
        x = sampling()
        x,b = global_optimize(kwargs,l=0,x=x,maxit=20)
        key_rate = keyrate(x,0,kwargs)
        if(key_rate>best_kr):
            best_kr = key_rate
            best_x = x
    return best_x
    