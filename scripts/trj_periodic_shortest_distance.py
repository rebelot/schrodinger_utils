from schrodinger.application.desmond.packages import traj_util, topo, traj
from schrodinger.structutils.measure import get_shortest_distance
from itertools import product
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Compute the shortest distance between two sets of atoms across orthorombic PBC"
    )
    parser.add_argument("cms", help="cms input file")
    parser.add_argument("-t", help="trajectory dir")
    parser.add_argument(
        "-g1", help="Define first group of atoms", metavar="ASL", required=True
    )
    parser.add_argument(
        "-g2", help="Define second gruoup of atoms", metavar="ASL", required=True
    )
    parser.add_argument("-periodic", help="Enable PBCs", action="store_true")
    parser.add_argument("-o", help="Save results in file", metavar="FILE")
    parser.add_argument("-p", help="Create plot", action="store_true")
    parser.add_argument("-s", help='Skip trajectory, format is "START:END:STEP"')
    parser.add_argument(
        "-noself",
        help="Do not compute distances within the same cell. This is useful to calculate the distance between a molecule and its periodic images",
        action="store_true",
    )
    args = parser.parse_args()

    slicer = (
        slice(*[int(x) if x else None for x in args.s.split(":")])
        if args.s
        else slice(None, None)
    )

    if args.t:
        msys, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.t)
    else:
        msys, cms, trj = traj_util.read_cms_and_traj(args.cms)

    trj = trj[slicer]

    g1_aids = cms.select_atom(args.g1)
    g1_gids = topo.aids2gids(cms, g1_aids)
    g1_st = cms.extract(g1_aids)
    g1_atoms = list(g1_st.atom)

    g2_aids = cms.select_atom(args.g2)
    g2_gids = topo.aids2gids(cms, g2_aids)
    g2_st = cms.extract(g2_aids)
    g2_atoms = list(g2_st.atom)

    _i = (0, 1, -1)  # --> (0, 0, 0), (0, 0, 1) ...
    CELLS = np.fromiter(product(_i, _i, _i), dtype=int)
    CELLS = CELLS if not args.noself else CELLS[1:]

    def get_periodic_images(P, box):
        for cell in CELLS * box:
            yield P + cell

    res = []
    for fr in trj:
        g1_st.setXYZ(fr.pos(g1_gids))

        if args.periodic:
            dist_imgs = []
            for g2_image in get_periodic_images(fr.pos(g2_gids), fr.box.diagonal()):
                g2_st.setXYZ(g2_image)
                dist_imgs.append(get_shortest_distance(g1_st, st2=g2_st))
            res.append(min(dist_imgs))

        else:
            g2_st.setXYZ(fr.pos(g2_gids))
            res.append(get_shortest_distance(g1_st, st2=g2_st))

    res = np.array(res)

    out = sys.stdout if not args.o else open(args.o + ".dat", "w")
    for (dist, a1, a2), fr in zip(res, trj):
        a1, a2 = int(a1) - 1, int(a2) - 1
        out.write(
            f"{fr.time} {dist} {g1_atoms[a1].pdbres}{g1_atoms[a1].resnum} ({g1_atoms[a1].index}) {g2_atoms[a2].pdbres}{g2_atoms[a2].resnum} ({g2_atoms[a2].index})\n"
        )
    out.close()

    if args.p:
        o = args.o if args.o else "trj_shortes_periodic_distance"
        plt.plot([fr.time / 1000 for fr in trj], res[:, 0])
        plt.xlabel("time (ns)")
        plt.ylabel("distance (Å)")
        plt.savefig(o + ".png")


if __name__ == "__main__":
    main()
