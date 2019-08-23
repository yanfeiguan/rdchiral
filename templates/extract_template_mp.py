#!/usr/bin/python
from __future__ import print_function, absolute_import

from multiprocessing import Process, Manager
#######################################################################
# Molecular conformer generator in progress
# genConf.py -isdf file_input.sdf -osdf file_output.sdf
# -n number_of_conformers (optional, if not specified is based
# on the nomber of rotable bonds) -rtpre rms_threshold_pre_opt(optional)
# -rtpost rms_threshold_post_opt(optional) -e energy_window (optional, Kcal/mol)
# -t number_of_threads (if not specify 1)
#######################################################################

## known issues / to-do list
## logging doesn't work properly
## need to print out rotatable bond atoms for reference
## comparison of printing vs. old full monte - more options
## comparison of printing vs. GoodVibes
## tabulation of results in addition to sdf output file
## pythonify where possible
## currently RMS but torsion_list would be better?

from rdkit import Chem
from rdkit.Chem import AllChem
from concurrent import futures
import argparse, logging, os, sys, time, copy
import gzip
import pickle

from template_extractor import extract_from_reaction

# PHYSICAL CONSTANTS
GAS_CONSTANT, PLANCK_CONSTANT, BOLTZMANN_CONSTANT, SPEED_OF_LIGHT, AVOGADRO_CONSTANT, AMU_to_KG, atmos = 8.3144621, 6.62606957e-34, 1.3806488e-23, 2.99792458e10, 6.0221415e23, 1.66053886E-27, 101.325
# UNIT CONVERSION
j_to_au = 4.184 * 627.509541 * 1000.0

# version number
__version__ = "1.0.1"

# Formatted output to command line and log fill

#wrap the genConf in process so that the genConf can be stopped
class extractor:
    def __init__(self, entry, timeout):
        reactants = '.'.join([x['smiles'] for x in entry['reactants']])
        products = '.'.join([x['smiles'] for x in entry['products']])
        db_id = entry['db_id']

        self.reaction = {'reactants': reactants, 'products': products,
                         'db_id': db_id}


        self.name = m.GetProp('_Name')
        self.timeout = timeout

    def __call__(self):
        self.return_dict = Manager().dict()
        self.process = Process(target=extract_from_reaction, args=(self.reaction, self.return_dict))

        self.process.start()
        self.process.join(args.timeout)
        if 'template' in self.return_dict:
            entry.update(self.return_dict['template'])
            return entry
        else:
            self.terminate()
            return None

    def terminate(self):
        self.process.terminate()

# conformational search / handles parallel threads if more than one structure is defined
def extract_templates(df, out, args):
    post_df = []
    with futures.ProcessPoolExecutor(max_workers=args.threads) as executor:
        tasks = [extractor(next(df), args.timeout) for m in range(args.threads)]
        running_pool = {task.reaction.db_id: executor.submit(task) for task in tasks}

        while True:
            if len(running_pool) == 0: break

            if len(post_df) % 10000 == 0 and post_df:
                with open(out, 'wb') as handle:
                    pickle.dump(post_df, handle)

            new_tasks = []
            for mol_id in list(running_pool):
                future = running_pool[mol_id]
                if future.done():
                    entry_post = future.result(timeout=0)
                    if entry_post:
                        post_df.append(entry_post)
                    else:
                        #warning failed
                        pass

                    #add new task
                    del(running_pool[mol_id])

                    try:
                        task = extractor(next(df), args.timeout)
                    except StopIteration:
                        #reach end of the supp
                        pass
                    else:
                        running_pool[task.reaction.db_id] = executor.submit(task)

    return post_df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Molecular conformer generator')
    parser.add_argument('-ipickle', required=True,
                        help='pickle input file')
    parser.add_argument('-timeout', required=False, default=30,
                        help='timewall for each job')
    parser.add_argument('-threads', required=False, default=48,
                        help='number of threads required')
    args = parser.parse_args()

    name = args.icsv.split('.')[0]
    out = name+'_extract'+'.pikcle'

    with open(args.ipickle, 'rb') as handle:
        df = pickle.load(handle)

    df = (x for x in df)
    post_df = extract_templates(df, out, args)

    with open (out, 'wb') as handle:
        pickle.dump(post_df, handle)
