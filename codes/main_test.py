# from codes.utils.respoonse_rate import Qcorr
from omegaconf import OmegaConf
import sys
# from codes.utils.progress_bar import progress_bar
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from codes.runs.simulation import simulator
from codes.runs.graph import output_graphs
from codes.utils.keyrate import keyrate
from codes.utils.keyrate_decoy_Wang import keyrate as keyrate_decoy
from codes.utils.keyrate_DPFK import keyrate as keyrate_DPFK
from codes.utils.keyrate_CP_Ma import keyrate as keyrate_CP_Ma
from codes.utils.keyrate_CPFK_mine import keyrate as keyrate_CPFK_mine
from codes.utils.keyrate_DP_Cao import keyrate as keyrate_DP_Cao

# from root import ROOT_PATH
if __name__ == '__main__':
    conf = 'config_main'
    args = OmegaConf.load('./conf/%s.yaml'%conf)
    args['conf'] = conf
    args['mode'] = 'CP'
    args['lemma'] = 'lemmaA1'
    args['N'] = 16
    args['Ntot'] = 1e21
    # output_graphs(args)
    # simulator(True, args)
    x = [0.660487544453102, 1e-05, 0.32833826346592443, 0.8941410211552083, 0.6761783356789675]

    # keyrate(x,0,args)
    # x =  [4.076347489223586e-06, 2.9444750131335518e-05, 0.3079141225269898, 0.998700379564296, 0.0008470767214430186]
    x = [0.6,0.2,0.99,0.99,0.5]
    # keyrate_Wang(x,175,args,ifdefault = True)
    X = range(268,300)
    # p = [max([keyrate_Wang(x,l,args,ifdefault = False),0]) for l in X]
    # args['mode'] = 'DP'
    # args['Ntot'] = 1e16
    import matplotlib.pyplot as plt
    # q = [max([keyrate_DP_Cao(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'DP_Cao')
    # args['Ntot'] = 1e20
    # q = [max([keyrate_CPFK_mine(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'CPFK_mine, Ntot=1e20')
    # args['Ntot'] = 1e21
    # q = [max([keyrate_CPFK_mine(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'CPFK_mine, Ntot=1e21')
    # args['Ntot'] = 1e22
    # q = [max([keyrate_CPFK_mine(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'CPFK_mine, Ntot=1e22')
    # args['Ntot'] = 1e23
    # q = [max([keyrate_CPFK_mine(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'CPFK_mine, Ntot=1e23')
    # args['Ntot'] = 1e24
    # q = [max([keyrate_CPFK_mine(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'CPFK_mine, Ntot=1e24')
    args['Ntot'] = 1e25
    q = [max([keyrate_CPFK_mine(x,l,args),0]) for l in X]
    plt.semilogy(X,q,label = 'CPFK_mine, Ntot=1e25, cpstates=16')
    args['cpstates'] = 20
    q = [max([keyrate_CPFK_mine(x,l,args),0]) for l in X]
    plt.semilogy(X,q,label = 'CPFK_mine, Ntot=1e25, cpstates=20')
    print(q)
    q = [max([keyrate_CP_Ma(x,l,args),0]) for l in X]
    plt.semilogy(X,q,label = 'CP_Ma')
    
    # q = [max([keyrate_decoy(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'decoy')
    # args['N'] = 12
    # q = [max([keyrate_DPFK(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'DPFK N=12')
    # args['N'] = 10
    # q = [max([keyrate_DPFK(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'DPFK N=10')
    # args['N'] = 8
    # q = [max([keyrate_DPFK(x,l,args),0]) for l in X]
    # plt.semilogy(X,q,label = 'DPFK N=8')
    
    
    # plt.semilogy(X,p,label = 'decoy')
    # plt.semilogy(X,q)
    plt.legend()
    plt.show()
    # print(p)