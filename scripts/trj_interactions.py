import argparse
import sys
from itertools import chain

import matplotlib.pyplot as plt
import numpy as np
from schrodinger.application.desmond.packages import analysis, topo, traj, traj_util
from schrodinger.structutils import analyze
from tabulate import tabulate

TITLE = {
    "HB": "H-Bonds",
    "SB": "Salt Bridges",
    "PiPi": "Pi-Pi interactions",
    "CatPi": "Pi-Cation interactions",
    "HPho": "Hydrophobic interactions",
    "all": "All interactions",
}


class InteractionOutput:
    """
    IO = InteractionOutput(data, cms)
    bonddict = IO.bonddict(aids1)
    em, keys = IO.existence_map(bonddict)
    nbonds = IO.count_bonds()
    """

    def __init__(self, data, btype, cms):
        self.data = self._preproc(data)
        self.cms = cms
        self.btype = btype
        self.atoms = list(cms.atom)

    def _preproc(self, data):
        pdata = []
        for frame in data:
            pbond = []
            for bond in frame:
                if isinstance(bond, analysis.CatPiInteraction):
                    x = bond.ring[0]
                    y = bond.cations[0]
                elif isinstance(bond, analysis.PiPiInteraction):
                    x = bond.ring1[0]
                    y = bond.ring2[0]
                else:
                    x, y = bond
                pbond.append([x, y])
            pdata.append(pbond)
        return pdata

    def aid2key(self, aid, start, stop):
        """
        key = [atom.chain, atom.pdbres, atom.resnum, atom.element, atom.index]
        :aid: 1-based atom indices
        :return: key[start:stop]
        """
        atom = self.atoms[aid - 1]
        key = [atom.chain, atom.pdbres, atom.resnum, atom.element, atom.index]
        return tuple(key[start:stop])

    def bonddict(self, aids1, i=0, j=None, k=0, l=None):
        """
        :start: and :stop: will slice the following atom `key` to unique-fy bonds
        key = [atom.chain, atom.pdbres, atom.resnum, atom.element, atom.index]
        """
        bonddict = {}
        for f, frame in enumerate(self.data):
            if frame:
                for bond in frame:
                    a1, a2 = bond if bond[0] in aids1 else reversed(bond)
                    key1 = self.aid2key(a1, i, j)
                    key2 = self.aid2key(a2, k, l)
                    key = (key1, key2)
                    bonddict.setdefault(key, [])
                    bonddict[key].append(f)
        if self.btype == "all":
           return { k: sorted(list(set(v))) for k, v in bonddict } 
        return bonddict

    def existence_map(self, bonddict):
        keys = list(bonddict.keys())
        em = np.zeros((len(keys), len(self.data)))

        for i, frames in enumerate(bonddict.values()):
            em[i, frames] = 1

        sorted_indices = np.argsort(em.mean(axis=1))

        return em[sorted_indices, :], [keys[i] for i in sorted_indices]

    def count_bonds(self):
        return [len(b) for b in self.data]

def plot_em(ax, em, keys, btype, times):
    ax.imshow(em, aspect="auto")
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([bkey2str(key) for key in keys], fontdict={"fontsize": 5})
    nlab = 5
    step = int(len(times) / (nlab - 1))
    pos = np.arange(0, len(times), step)
    ax.set_xticks(pos)
    ax.set_xticklabels(np.round(times[::step], decimals=2))
    ax.set_xlabel("time (ns)")
    ax.set_title(TITLE[btype])
    # ax.set_ylabel('bond')


def plot_allbonds(ax, allbonds, btype, times):
    ax.set_title(TITLE[btype])
    ax.plot(times, allbonds)
    ax.set_xlabel("time (ns)")
    ax.set_ylabel("# bonds")


def write_em(fname: str, em, keys, times):
    with open(fname, "w") as fh:
        tstr = ",".join(str(t) for t in times)
        fh.write("bond," + tstr + "\n")
        for bond, key in zip(em, keys):
            frames = ",".join(str(v) for v in bond)
            fh.write(f"{bkey2str(key)}, {frames}\n")


def write_allbonds(fname, all_bonds, times):
    with open(fname, "w") as fh:
        for t, bond_num in zip(times, all_bonds):
            fh.write(f"{t} {bond_num}\n")


def write_summary(fname, em, keys, btype, args):
    tdata = [
        [
            " ".join([str(i) for i in key[0]]),
            " ".join([str(i) for i in key[1]]),
            f"{row.mean():.2%}",
        ]
        for key, row in zip(keys, em)
    ]
    table = tabulate(reversed(tdata), headers=(args.g1, args.g2, "%"))
    with open(fname, "w") as fh:
        fh.write(f"{args.out} {TITLE[btype]} summary.\n")
        fh.write(table)
        fh.write("\n\n")


def bkey2str(key):
    return (
        f'{" ".join([str(i) for i in key[0]])} - {" ".join([str(i) for i in key[1]])}'
    )


def main():
    parser = argparse.ArgumentParser(
        description="""Calculate all nonbonded interactions within or
                between the specified atoms.  If group2 is given, then this
                script will return interacions between the two groups of
                atoms.  Else, nonbonded interactions within a single group of
                atoms will be returned."""
    )
    parser.add_argument("cms", help="Input cms file")
    parser.add_argument(
        "out",
        help="output base filename (descriptive extensions are added automatically)",
    )
    parser.add_argument("-t", help="Trajectory dir")
    parser.add_argument("-g1", help="define group 1", metavar="ASL", default="all")
    parser.add_argument("-g2", help="define group 2", metavar="ASL", default=None)
    parser.add_argument(
        "-p", "--plot", help="plot bond nr over time", action="store_true"
    )
    parser.add_argument(
        "-rkey",
        help="comma-separated indexes to slice the key that is used to identify unique reisudes:\
        key = [atom.chain, atom.pdbres, atom.resnum, atom.element, atom.index] \
        they may be specified indepentently for both groups in the order in which they appear:\
        e.g.: -rkey 0,3 -rkey 0,0 (0 means all)",
        action="append",
    )
    parser.add_argument(
        "-r",
        help="Filter results at residue level (alias for '-rkey 0,3[,0,3])",
        action="store_true",
    )
    parser.add_argument(
        "-b",
        help="comma-separated list of interaction types to compute:\
            default is hb,sb,pipi,catpi [optional: hpho]",
        default="hb,sb,pipi,catpi",
    )
    parser.add_argument(
        "-th",
        help="only report the existence of bonds with an existence % above N",
        metavar="N",
        type=float,
    )
    # parser.add_argument("-em", help="output indexed existence map", action="store_true")
    parser.add_argument(
        "-s", help="slice trajectory", metavar="START:END:STEP", default="::"
    )

    args = parser.parse_args()

    if args.t:
        msys, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.t)
    else:
        msys, cms, trj = traj_util.read_cms_and_traj(args.cms)

    slicer = slice(*[int(v) if v else None for v in args.s.split(":")])
    trj = trj[slicer]
    times = np.array(list(fr.time / 1000.0 for fr in trj))

    asl1 = args.g1
    asl2 = args.g2 or asl1
    aids1 = cms.select_atom(asl1)
    aids2 = cms.select_atom(asl2)

    analyzers = []
    btype = []
    for a in args.b.split(","):
        if a == "hb":
            analyzers.append(
                analysis.HydrogenBondFinder(msys, cms, aids1=aids1, aids2=aids2)
            )
            btype.append("HB")
        elif a == "sb":
            analyzers.append(
                analysis.SaltBridgeFinder(msys, cms, aids1=aids1, aids2=aids2)
            )
            btype.append("SB")
        elif a == "pipi":
            analyzers.append(analysis.PiPiFinder(msys, cms, aids1=aids1, aids2=aids2))
            btype.append("PiPi")
        elif a == "catpi":
            analyzers.append(analysis.CatPiFinder(msys, cms, aids1=aids1, aids2=aids2))
            btype.append("CatPi")
        elif a == "hpho":
            analyzers.append(
                analysis.HydrophobicInter(msys, cms, prot_asl=asl1, lig_asl=asl2)
            )
            btype.append("HPho")

    out = analysis.analyze(trj, *analyzers)

    out.append(chain(*out))
    btype.append("all")

    if args.rkey:
        rkeys = [
            int(key) if int(key) != 0 else None
            for keys in args.rkey
            for key in keys.split(",")
        ]
        assert len(rkeys) == 2 or len(rkeys) == 4
        rkeys[2:4] = rkeys[:3] if (len(rkeys) == 2 and not args.g2) else rkeys[2:4]
    else:
        rkeys = [0, None, 0, None]

    if args.r:
        rkeys = [0, 3, 0, 3]

    for data, btype in zip(out, btype):
        InterOut = InteractionOutput(data, btype, cms)
        bonddict = InterOut.bonddict(aids1, *rkeys)
        em, keys = InterOut.existence_map(bonddict)
        allbonds = InterOut.count_bonds()

        if args.th:
            existence = em.mean(axis=1)
            indexes = np.nonzero(existence * 100 >= args.th)[0]
            em, keys = em[indexes, :], [keys[i] for i in indexes]

        f_allbonds: str = args.out + f"_{btype}-ab.csv"
        f_em: str = args.out + f"_{btype}-em.csv"
        f_summary: str = args.out + f"_{btype}.txt"

        write_summary(f_summary, em, keys, btype, args)
        write_em(f_em, em, keys, times)
        write_allbonds(f_allbonds, allbonds, times)

        if keys and args.plot:

            fig, ax = plt.subplots(1, figsize=(10, 0.1 * len(keys) + 1))
            plot_em(ax, em, keys, btype, times)
            plt.tight_layout()
            fig.savefig(args.out + f"_{btype}-em.png")
            fig, ax = plt.subplots(1)
            plot_allbonds(ax, allbonds, btype, times)
            fig.savefig(args.out + f"_{btype}-ab.png")




if __name__ == "__main__":
    main()
