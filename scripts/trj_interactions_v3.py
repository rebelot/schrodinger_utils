from schrodinger.application.desmond.packages import topo, traj_util, traj, analysis
# analysis._PiPiInteraction
# analysis._PiCatInteraction
# from schrodinger.structutils import analyze
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse
from dataclasses import dataclass


@dataclass
class Residue:
    number: int
    name: str
    chain: str


def aid2info(aid, atoms):
    atom = atoms[aid]
    return Residue(atom.resnum, atom.pdbres.strip(), atom.chain)


def order_bonds(kind, out, aids1, atoms):
    result = []
    for frame in out:
        bonds = []
        for bond in frame:
            if kind in ['HB',   'SB']:
                x, y = bond
            elif kind == 'PP':
                x = bond.ring1[0]
                y = bond.ring2[0]
            else:  # kind == 'CP'
                x = bond.ring[0]
                y = bond.cations[0]
            if x in aids1:
                bonds.append((aid2info(x, atoms), aid2info(y, atoms)))
            else:
                bonds.append((aid2info(y, atoms), aid2info(x, atoms)))
        result.append(bonds)
    return result


def count_bonds(out):
    return [len(frame) for frame in out]


def reduce_residues(out):
    resdict = {}
    for i, frame in enumerate(out):
        if frame:
            for bond in frame:
                r1, r2 = bond
                key = f'{r1.name}{r1.number}_{r1.chain}-{r2.name}{r2.number}_{r2.chain}'
                resdict.setdefault(key, [])
                resdict[key].append(i)

    em = np.zeros((len(resdict.keys()), len(out)))

    for i, frames in enumerate(resdict.values()):
        em[i, frames] = 1

    sorted_indices = np.argsort(em.mean(axis=1))

    reslist = list(resdict.keys())
    # noinspection PyTypeChecker
    return em[sorted_indices, :], [reslist[i] for i in sorted_indices]


def plot_em(ax, em, keys, times):
    ax.imshow(em, aspect='auto')
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([bond for bond in keys], fontdict={'fontsize': 5})
    nlab = 5
    step = int(len(times) / (nlab - 1))
    pos = np.arange(0, len(times), step)
    ax.set_xticks(pos)
    ax.set_xticklabels(np.round(times[::step], decimals=2))
    ax.set_xlabel('time (ns)')
    # ax.set_ylabel('bond')


def plot_allbonds(ax, allbonds, times):
    ax.plot(times, allbonds)
    ax.set_xlabel('time (ns)')
    ax.set_ylabel('# bonds')


def bond2str(bond_key):
    return f'{" ".join(str(v) for v in bond_key[0])} - {" ".join(str(v) for v in bond_key[1])}'


def write_output(fh, em, keys):
    pass


def write_em(fh, em, keys, times):
    tstr = ','.join(str(t) for t in times)
    fh.write('bond,'+tstr+'\n')
    for bond, key in zip(em, keys):
        frames = ','.join(str(v) for v in bond)
        fh.write(f'{bond2str(key)}, {frames}\n')


def write_allbonds(fh, all_bonds, times):
    for t, bond_num in zip(times, all_bonds):
        fh.write(f'{t} {bond_num}\n')


def main():
    parser = argparse.ArgumentParser(
        description="""Calculate all nonbonded interactions within or
                between the specified atoms.  If group2 is given, then this
                script will return interacions between the two groups of
                atoms.  Else, nonbonded interactions within a single group of
                atoms will be returned.""")
    parser.add_argument("cms", help="Input cms file")
    parser.add_argument('-t', help="Trajectory dir")
    parser.add_argument("-g1", help="define group 1",
                        metavar="ASL", default='all')
    parser.add_argument("-g2", help="define group 2",
                        metavar="ASL", default=None)
    parser.add_argument("-p", "--plot", help="plot bond nr over time",
                        action="store_true")
    parser.add_argument(
        "-o", help="output base filename (descriptive extensions are added automatically)", metavar="FILE")
    # parser.add_argument("-em", help="output indexed existence map", action="store_true")
    parser.add_argument('-s', help='slice trajectory',
                        metavar='START:END:STEP', default='::')

    args = parser.parse_args()

    if args.t:
        msys, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.t)
    else:
        msys, cms, trj = traj_util.read_cms_and_traj(args.cms)

    slicer = slice(*[int(v) if v else None for v in args.s.split(':')])
    trj = trj[slicer]
    times = np.array(list(int(round(fr.time/1000.0)) for fr in trj))

    atoms = list(cms.atom)

    asl1 = args.g1
    asl2 = args.g2 or asl1
    aids1 = cms.select_atom(asl1)
    aids2 = cms.select_atom(asl2)
    # TODO: implement sidechain/backbone distinction
    # aids1_b = cms.select_atom(asl1 + ' and backbone')
    # aids2_b = cms.select_atom(asl2 + ' and backbone')
    # aids1_b = cms.select_atom(asl1 + ' and sidechain')
    # aids2_b = cms.select_atom(asl2 + ' and sidechain')

    analyzers = [
        analysis.HydrogenBondFinder(msys, cms, aids1=aids1, aids2=aids2),
        analysis.SaltBridgeFinder(msys, cms, aids1=aids1, aids2=aids2),
        analysis.PiPiFinder(msys, cms, aids1=aids1, aids2=aids2),
        analysis.CatPiFinder(msys, cms, aids1=aids1, aids2=aids2)]
    # analysis.WaterBridges(msys, cms, prot_asl, lig_asl),
    # analysis.HydrophobicInter(msys, cms, prot_asl, lig_asl),
    # analysis.MetalInter(msys, cms, prot_als, lig_asl, metal_asl)

    out = analysis.analyze(trj, *analyzers)

    for o, btype in zip(out, ('HB', 'SB', 'PP', 'CP')):
        out_ordered = order_bonds(btype, o, aids1, atoms)
        em, keys = reduce_residues(out_ordered)

        all_bonds = count_bonds(o)

        # fh_summary = open(args.o + '.txt', 'w') if args.o else sys.stdout
        fh_allbonds = open(args.o + f'_{btype}' + '-cum.dat', 'w') if args.o else sys.stdout
        fh_em = open(args.o + f'_{btype}' + '-em.dat', 'w') if args.o else sys.stdout

        # write_summary(fh_summary, out_sorted)
        write_em(fh_em, em, keys, times)
        fh_em.close()

        write_allbonds(fh_allbonds, all_bonds, times)
        fh_allbonds.close()

        if keys and args.plot:

            fig, ax = plt.subplots(1, figsize=(10, 0.1 * len(keys) + 1))
            plot_em(ax, em, keys, times)

            if btype == 'HB':
                plt.title('H-Bonds')
            elif btype == 'SB':
                plt.title('Salt Bridges')
            elif btype == 'PP':
                plt.title('Pi-Pi Interactions')
            elif btype == 'CP':
                plt.title('Pi-Cation Interactions')
            plt.tight_layout()
            fig.savefig(args.o + f'-em_{btype}.png')
            fig, ax = plt.subplots(1)
            plot_allbonds(ax, all_bonds, times)
            fig.savefig(args.o + f'-ab_{btype}.png')


if __name__ == "__main__":
    main()
