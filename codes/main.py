# from codes.utils.respoonse_rate import Qcorr
from omegaconf import OmegaConf
import sys
from os import path
sys.path.append(path.dirname(path.dirname(path.abspath(__file__))))
from codes.runs.simulation import simulator
from codes.runs.graph import output_graphs
import argparse
# from root import ROOT_PATH
if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('--conf', type=str, default='config1')
    arg_parser.add_argument(
        '--draw', help='draw based on exist data, the database name is required', default='')
    arg_parser.add_argument(
        '--run', type=bool, default=False, help='run with optimizers or fixed parameters')
    args = arg_parser.parse_args()

    if args.draw:
        output_graphs(args)
    else:
        conf = args.conf
        args = OmegaConf.load('./conf'+conf+'.yaml')
        args['conf'] = conf
        simulator(args.run, args)
