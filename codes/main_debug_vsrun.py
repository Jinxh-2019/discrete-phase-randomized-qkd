# from codes.utils.respoonse_rate import Qcorr
from omegaconf import OmegaConf
import sys
# from codes.utils.progress_bar import progress_bar
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from codes.runs.simulation import simulator
from codes.runs.graph import output_graphs
# from root import ROOT_PATH
if __name__ == '__main__':
    conf = 'config_main'
    args = OmegaConf.load('./conf/%s.yaml'%conf)
    args['conf'] = conf
    # output_graphs(args)
    
    simulator(True, args)

