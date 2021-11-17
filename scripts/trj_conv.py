from schrodinger.application.desmond.packages import topo, traj, traj_util
from schrodinger.structure import PDBWriter
import argparse
from itertools import chain

def main():
    parser = argparse.ArgumentParser(description='Convert Desmond trajectory and topology from cms/dtr to pdb/xtc')
    parser.add_argument('cms', help='topology')
    parser.add_argument('-t', help='trajectory')
    parser.add_argument('-o', help='output base name', required=False)
    parser.add_argument('-s', help='slice trajectory', metavar='START:END:STEP')
    args = parser.parse_args()

    basename = args.o if args.o else args.cms.split['-out.cms'][0]

    if args.t:
        msys, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.t)
    else:
        msys, cms, trj = traj_util.read_cms_and_traj(args.cms)

    slicer = (
        slice(*[int(i) if i else None for i in args.s.split(":")])
        if args.s
        else slice(None, None)
    )
    trj = trj[slicer]

    # pseudo_aids = chain.from_iterable(cms.pseudoatoms.values())
    # pseudo_gids = topo.aids2gids(cms, pseudo_aids)

    new_trj = [fr.reduce(cms.allaid_gids) for fr in trj]

    traj.write_traj(new_trj, basename + '.xtc', format=traj.Fmt.XTC)
    writer = PDBWriter(basename + '.pdb', reorder_by_sequence=True)
    writer.write(cms)

if __name__ == "__main__":
    main()
