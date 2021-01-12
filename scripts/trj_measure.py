from schrodinger.application.desmond.packages import topo, traj, traj_util, analysis
import argparse
import numpy as np
import typing

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument
    parser.add_argument("cms", help="Input cms file")
    parser.add_argument("trj", help="Input trajectory dir")
    parser.add_argument("out", help="output basename")
    parser.add_argument( "-e", action="append",
            help="""Expression in the form of ASL_1[:ASL_N]:MODE, where:
ASL is the ':'-separated list of asl expressions specifying the atoms to measure, if the asl expands to more than one atom,
the CenterOfMass of the selection is used; distance -> 2 tokens, angle -> 3 tokens, dihedral -> 4 or 6 tokens, xyz -> N tokens
MODE = distance, angle, dihedral or xyz (only output positons); 
Multiple expressions can be provided and will be executed sequentially from left to right.
e.g: -e "a.n 10:a.n 20-25:distance" will compute the distance between atoms 10 and the center of mass of atoms 20-25.""")
    parser.add_argument( "-pp",
        help='''Write the body of custom python function to be called upon the Nth expression results (usually those where MODE=xyz)
        in the form "N@body". The function scope exposes the variables: :pos: (N_tokens x N_frames x 3) numpy array; :fr: the frame object''', default=None)
    args = parser.parse_args()

    msys, cms = topo.read_cms(args.cms)
    trj: typing.List[traj.Frame] = traj.read_traj(args.trj)

    analyzers = []
    for e in args.e:
        e = e.split(':')
        mode = e.pop(-1)
        tokens_aids = [cms.select_atom(asl) for asl in e]
        token_gids = [topo.aids2gids(cms, aids) for aids in tokens_aids]
        tokens = [analysis.Com(msys, cms, gids=gids) for gids in token_gids]
        if mode == 'distance':
            analyzers.append(analysis.Distance(msys, cms, *tokens))
        elif mode == 'angle':
            analyzers.append(analysis.Angle(msys, cms, *tokens))
        elif mode == 'dihedral':
            if len(tokens) == 4:
                tokens = tokens[:3] + tokens[1:]
            analyzers.append(analysis.PlanarAngle(msys, cms, *tokens))
        elif mode == 'xyz':
            analyzers.append(*tokens)

    res = analysis.analyze(trj, *analyzers)

    for pp in args.pp:
        pass


if __name__ == "__main__":
    print('Not implemented yet. Goodbye!')
    # main()
