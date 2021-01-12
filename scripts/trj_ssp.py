#!/opt/schrodinger/suites2019-3/run
# -*- coding: utf-8 -*-

import argparse
import sys

import matplotlib.pyplot as plt
import numpy as np

from schrodinger.application.desmond.packages import analysis, topo, traj_util

def main():
    parser = argparse.ArgumentParser(
        description='Calculate Secondary Structure along trajectory.')
    parser.add_argument('cms', help="Input cms file.")
    parser.add_argument(
        '-ASL', help="Atom selection for secondary structure calculation, default is 'protein'", default='protein')
    parser.add_argument('-o', help="Output filename.")
    parser.add_argument('-p', help="Plot results", action='store_true')
    args = parser.parse_args()

    msys, cms, trj = traj_util.read_cms_and_traj(args.cms)
    aids = cms.select_atom(args.ASL)
    ss_analyzer = analysis.SecondaryStructure(msys, cms, aids)
    rkey, ss = analysis.analyze(trj, ss_analyzer)
    ss = np.array(ss)

    if args.p:
        fig, ax = plt.subplots(2, 1, sharex=True)
        ax[0].plot([fr.time for fr in trj], [
                   100 * sse[sse > 0].sum() / ss.shape[1] for sse in ss])
        ax[0].set_ylabel('SSE %')
        ax[1].pcolormesh([fr.time for fr in trj], [int(r[6:]) for r in rkey], ss.T)
        ax[1].set_xlabel('time (ps)')
        ax[1].set_ylabel('residue nr.')
        plt.savefig(args.o + '.png')

    fh = open(args.o + '.dat', 'w') if args.o else sys.stdout
    for fr, sse in zip(trj, ss):
        fh.write(f"{fr.time} {sse[sse > 0].sum() / len(sse)}\n")
    # TODO: make some effort to store the existence matrix in a more convenient format!
    for fr, sse in zip(trj, ss):
        fh.write(f"{fr.time}, {', '.join(str(x) for x in sse)}\n")
    fh.close()

if __name__ == "__main__":
    main()
