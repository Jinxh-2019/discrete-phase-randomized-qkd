# from codes.utils.respoonse_rate import Qcorr
from omegaconf import OmegaConf
import sys
# from codes.utils.progress_bar import progress_bar
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from codes.runs.simulation import simulator
from codes.runs.graph import output_graphs
from codes.utils.keyrate import keyrate
from codes.utils.keyrate_decoy_Wang import keyrate as keyrate_Wang
# from root import ROOT_PATH
if __name__ == '__main__':
    conf = 'config_main'
    args = OmegaConf.load('./conf/%s.yaml'%conf)
    args['conf'] = conf
    args['mode'] = 'decoy'
    args['lemma'] = 'lemmaA1'
    args['Ntot'] = 1e16
    # output_graphs(args)
    # simulator(True, args)
    x = [0.660487544453102, 1e-05, 0.32833826346592443, 0.8941410211552083, 0.6761783356789675]

    # keyrate(x,0,args)
    # x =  [4.076347489223586e-06, 2.9444750131335518e-05, 0.3079141225269898, 0.998700379564296, 0.0008470767214430186]
    x = [0.6,0.2,1/3,0.5,0.5]
    # keyrate_Wang(x,175,args,ifdefault = True)
    X = range(180,400)
    p = [max([keyrate_Wang(x,l,args,ifdefault = False),0]) for l in X]
    args['mode'] = 'DP'
    # args['Ntot'] = 1e16
    q = [max([keyrate(x,l,args),0]) for l in X]
    
    import matplotlib.pyplot as plt
    plt.semilogy(X,p,label = 'decoy')
    plt.semilogy(X,q,label = 'DP')
    plt.legend()
    plt.show()
    # print(p)