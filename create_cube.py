#!/usr/bin/env python3

import argparse
from tempfile import NamedTemporaryFile
import subprocess
import re
import sys

# edit this template for fine tuning
template = '''\
1 # select type:
1 # MO
3 # select operator:
{spin} # alpha (0) or beta (1)
5 # select output format:
7 # cube
4 # select resolution
{grid} # ngrid
2 # select orbital number
{index} # orbital index
10 # create file
11 # quit program'''

def call_orca_plot(orb_idx, spin, gbw, grid):
    'call orca_plot for orbitals orb_idx with operator spin on file gbw with grid'
    for i in orb_idx:
        with NamedTemporaryFile(mode='w+') as f: # only write temporary file
            content = template.format(spin=spin, grid=grid, index=i)
            f.write(content)
            
            f.seek(0) # reset file pointer for reading
            subprocess.run(['orca_plot', gbw, '-i'], stdin=f, stdout=subprocess.DEVNULL)
            print('INFO: Plotted orbital {} with operator {} from {} using grid {}'.format(i, spin, gbw, grid))

def expand_range(inp_str):
    inp_split = inp_str.split("-")
    int_list = list(range(int(inp_split[0]), int(inp_split[1]) + 1))
    str_list = [str(i) for i in int_list]
    return str_list

def convert_list(inp_str):
    inp_split = inp_str.split(",")
    res_list = []

    for i in inp_split:
        if "-" in i:
            res_list.extend(expand_range(i))
        else:
            res_list.append(i)
    
    return res_list
        


def run():
    parser = argparse.ArgumentParser(
        description='create cube from gbw')
    parser.add_argument('orca_gbw', metavar='ORCA_GBW', 
        help='ORCA gbw file')
    parser.add_argument('-a', '--alpha', metavar='i_a',  nargs='?', default="",
        help='indices of alpha orbitals')
    parser.add_argument('-b', '--beta', metavar='i_b',  nargs='?', default="",
        help='indices of beta orbitals')
    parser.add_argument('-g', '--grid', metavar='N', default=60,
        help='number of grid points')
    args = parser.parse_args()

    if not (args.alpha or args.beta):
        sys.exit("Please give orbital indices for processing!")
    else:
        match_alpha = re.findall(r"[^\d,-]", args.alpha)
        match_beta = re.findall(r"[^\d,-]", args.beta)

        if match_alpha:
            sys.exit(f"Unsupported characters in alpha orbital list: {' '.join(match_alpha)}")
        if match_beta:
            sys.exit(f"Unsupported characters in beta orbital list: {' '.join(match_beta)}")

    grid = args.grid
    gbw = args.orca_gbw

    if args.alpha:
        la = convert_list(args.alpha)
        call_orca_plot(orb_idx=la, spin=0, gbw=gbw, grid=grid)

    if args.beta:
        lb = convert_list(args.beta)
        call_orca_plot(orb_idx=lb, spin=1, gbw=gbw, grid=grid)
        
if __name__ == '__main__':
    run()

