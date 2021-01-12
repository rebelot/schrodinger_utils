from schrodinger.application.desmond.packages import topo, traj_util, traj, analysis
from schrodinger.structutils import analyze
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse


def aid2info(aid, ATOMS):
    atom = ATOMS[aid]
    return atom.pdbres, atom.resnum, atom.chain, atom.element, atom.index


def sortbonds(out, aids1, ATOMS):
    return [[(aid2info(x, ATOMS), aid2info(y, ATOMS)) if x in aids1 else (aid2info(y, ATOMS), aid2info(x, ATOMS)) for x, y in Z] for Z in out]


def count_bonds(out):
    return [len(frame) for frame in out]


def reduceres(out):
    resdict = {}
    for i, frame in enumerate(out):
        for a1, a2 in frame:
            key = (a1[:-2], a2[:-2])
            resdict.setdefault(key, [])
            resdict[key].append(i)

    em = np.zeros((len(resdict.keys()), len(out)))

    for i, frames in enumerate(resdict.values()):
        em[i, frames] = 1

    sorted_indices = np.argsort(em.mean(axis=1))

    reslist = list(resdict.keys())
    return em[sorted_indices, :], [reslist[i] for i in sorted_indices]


def plot_em(ax, em, keys, times):
    ax.imshow(em, aspect='auto')
    ax.set_yticks(range(len(keys)))
    ax.set_yticklabels([bond2str(bond) for bond in keys], fontdict={'fontsize': 5})
    nlab = 5
    step = int(len(times) / (nlab - 1))
    pos = np.arange(0, len(times), step)
    ax.set_xticks(pos)
    ax.set_xticklabels(np.round(times[::step], decimals=2))
    ax.set_xlabel('time (ns)')
    ax.set_ylabel('bond')


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
        fh.write(f'{t},{bond_num}\n')

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
    times = np.array(list(fr.time for fr in trj))/1000

    ATOMS = list(cms.atom)

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
        analysis.SaltBridgeFinder(msys, cms, aids1=aids1, aids2=aids2)]
    # analysis.PiPiFinder(msys, cms, aids1=aids1, aids2=aids2),
    # analysis.CatPiFinder(msys, cms, aids1=aids1, aids2=aids2),
    # analysis.WaterBridges(msys, cms, prot_asl, lig_asl),
    # analysis.HydrophobicInter(msys, cms, prot_asl, lig_asl),
    # analysis.MetalInter(msys, cms, prot_als, lig_asl, metal_asl)

    out = analysis.analyze(trj, *analyzers)

    for o, btype in zip(out, ('HB', 'SB')):
        out_sorted = sortbonds(o, aids1, ATOMS)
        em, keys = reduceres(out_sorted)
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
            fig, ax = plt.subplots(1)
            plot_em(ax, em, keys, times)
            fig.savefig(args.o + f'-em_{btype}.png')

            fig, ax = plt.subplots(1)
            plot_allbonds(ax, all_bonds, times)
            fig.savefig(args.o + f'-ab_{btype}.png')


if __name__ == "__main__":
    main()
