from root import ROOT_PATH

from sql_exe import create_tables,write_records,connect
from codes.utils.optimize import global_optimize, optimize


def optimizer(kwargs):
    cur = connect(kwargs['conf'],kwargs)
    testitems = kwargs["testitems"]
    create_tables(cur,kwargs)
    for testitem in testitems:
        if testitem == 'DP':
            for logNtot in kwargs['Ntots']:
                Ntot = 10 ^ logNtot
                kwargs['Ntot'] = Ntot


def dot_optimizer(x,lgts,kwargs):
    mode = kwargs['mode']
    for l in lgts:
        if x !=[]:
            xini = x
        if mode in {'DP','DP_Cao','CP_Ma','CPFK'}:
            x,b = global_optimize(kwargs,xini,l=l)
        elif mode in {}:
            x,b = optimize(kwargs,xini,l=l)
        else:
            x = [4.75111801e-01, 7.81356923e-06, 4.08562954e-01, 2.68708743e-01, 6.67068160e-01]
            b = True
        if  x==[]:
            print('optimize error')
            continue
        kwargs['Lgt'] = l
        x = list(x)
        if mode=='DP':
            logNtot = kwargs['logNtot']
        keyrate = keyrate(x,l,kwargs)
        if keyrate>0:
            write_records('DP10^%d' % logNtot,x,l,keyrate)
        else:
            keyrate = None
        print(mode,l,keyrate,'\n',x)
        if keyrate<0:
            break