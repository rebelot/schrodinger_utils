#!/opt/schrodinger/suites2019-3/run

from schrodinger.application.desmond.packages import topo, traj_util, traj
from schrodinger.structutils import interactions
from schrodinger.protein.analysis import analyze
import numpy as np
import matplotlib.pyplot as plt
import sys
import argparse


def cml():
    parser = argparse.ArgumentParser(
               description="""Calculate all nonbonded interactions within or
                between the specified atoms.  If group2 is given, then this
                script will return interacions between the two groups of
                atoms.  Else, nonbonded interactions within a single group of
                atoms will be returned.""")
    parser.add_argument("cms", help="Input cms file")
    parser.add_argument("-g1", help="define group 1",
                        metavar="ASL", default=None)
    parser.add_argument("-g2", help="define group 2",
                        metavar="ASL", default=None)
    parser.add_argument("-p", "--plot", help="plot bond nr over time",
                        action="store_true")
    parser.add_argument("-o", help="output base filename (.txt and/or .png extensions are added automatically)", metavar="FILE")
    # parser.add_argument("-em", help="output indexed existence map", action="store_true")
    parser.add_argument('-s', help='slice trajectory', metavar='START:END:STEP', default='::')

    args = parser.parse_args()
    return args


def bond_counter(bonds_list):
    # there is one list for each frame, hence n = frames number
    n = len(bonds_list)
    d = {}
    for bonds in bonds_list:
        for bond in bonds:
            d.setdefault(bond, 0)
            d[bond] += 1 / n
    return d


def get_bonds_list(cms, trj, g1=None, g2=None, ses='::'):

    system_st = cms.extract(cms.select_atom('all'))

    group1 = analyze.evaluate_asl(system_st, g1)
    group2 = analyze.evaluate_asl(system_st, g2) if g2 else g2

    salt_bridges = []
    hydrogen_bonds = []

    fnum = len(trj)
    start, end, step = map(lambda x: int(x) if x else None, ses.split(':'))
    for fr in trj[start:end:step]:
        system_st.setXYZ(fr.pos())
        saltbr = interactions.get_salt_bridges(
            system_st, group1=group1, group2=group2,
            order_by=interactions.OrderBy.InputOrder)
        hbonds = analyze.hbond.get_hydrogen_bonds(
            system_st, atoms1=group1, atoms2=group2)
        salt_bridges.append(saltbr)
        hydrogen_bonds.append(hbonds)

        sys.stderr.write('\rframe %4d of %d' % (int(fr.orig_index), fnum))
        sys.stderr.flush()

    return salt_bridges, hydrogen_bonds


def print_output(salt_bridges, hydrogen_bonds, saltbr_dict, hbonds_dict, trj, em_hb, em_sb, fd):

    def atmfmt(atom):
        return "{}{:<4} {} {:>10}".format(
            atom.pdbres, atom.resnum, atom.chain,
            "(" + atom.element + " " + str(atom.index) + ")")

    fd.write("\n")
    fd.write(26 * "-" + " HBONDS " + 26 * "-" + '\n')
    for i, (bond, freq) in enumerate(sorted(hbonds_dict.items(), key=lambda x: x[1], reverse=True)):
        fd.write("{:4}   {}      ---      {}: {:>7.2%}\n".format(
            i, atmfmt(bond[0]), atmfmt(bond[1]), freq))

    fd.write("\n")
    fd.write(23 * "-" + " SALT BRIDGES " + 23 * "-" + '\n')
    for i, (bond, freq) in enumerate(sorted(saltbr_dict.items(), key=lambda x: x[1], reverse=True)):
        fd.write("{:4}   {}      ---      {}: {:>7.2%}\n".format(
            i, atmfmt(bond[0]), atmfmt(bond[1]), freq))


def makedat(salt_bridges, hydrogen_bonds, em_hb, em_sb, trj, out):
    fd_hbsbt = open(out + '-cumulative.dat', 'w')
    fd_hbsbt.write("# time sbr hb\n")
    for sbr, hb, fr in zip(salt_bridges, hydrogen_bonds, trj):
        fd_hbsbt.write(f"{fr.time} {len(sbr)} {len(hb)}\n")

    fd_hbem = open(out + '-HBEM.dat', 'w')
    fd_hbem.write("# HB Existence Map\n")
    for line in em_hb:
        fd_hbem.write(' '.join(str(i) for i in line))
        fd_hbem.write('\n')

    fd_sbem = open(out + '-SBEM.dat', 'w')
    fd_sbem.write("# SB Existence Map\n")
    for line in em_sb:
        fd_sbem.write(' '.join(str(i) for i in line))
        fd_sbem.write('\n')


def plot_data(salt_bridges, hydrogen_bonds, em_hb, em_sb, trj, out, ses='::'):
    start, end, step = map(lambda x: int(x) if x else None, ses.split(':'))
    fig, axs = plt.subplots(2, 1)
    axs[0].plot([fr.time for fr in trj[start:end:step]], [
                len(bonds) for bonds in salt_bridges])
    axs[0].set_title('salt bridges')
    axs[1].plot([fr.time for fr in trj[start:end:step]], [
                len(bonds) for bonds in hydrogen_bonds])
    axs[1].set_title('hydrogen bonds')
    fig.savefig(out + '.png')

    fig, axs = plt.subplots(2, 1)
    axs[0].imshow(em_sb, aspect='auto')
    axs[0].set_title('salt bridges')

    axs[1].imshow(em_hb, aspect='auto')
    axs[1].set_title('hydrogen bonds')
    fig.set_size_inches(20,30)
    fig.savefig(out + '-emap.png')

def create_existence_map(interaction_list):
    d = {}
    for i, bonds_per_frame in enumerate(interaction_list):
        for bond in bonds_per_frame:
            d.setdefault(bond, [])
            d[bond].append(i)
    em = np.zeros((len(d.keys()), len(interaction_list)))
    for i, pair_frindex in enumerate(sorted(d.items(), key=lambda x: len(x[1]), reverse=True)):
        em[i,pair_frindex[1]] = 1
    return em


def main():
    args = cml()

    _, cms, trj = traj_util.read_cms_and_traj(args.cms)

    salt_bridges, hydrogen_bonds = get_bonds_list(cms, trj, args.g1, args.g2, args.s)

    saltbr_dict = bond_counter(salt_bridges)
    hbonds_dict = bond_counter(hydrogen_bonds)

    em_sb = create_existence_map(salt_bridges)
    em_hb = create_existence_map(hydrogen_bonds)

    fd = open(args.o + '.txt', 'w') if args.o else sys.stdout
    print_output(salt_bridges, hydrogen_bonds, saltbr_dict, hbonds_dict, trj, em_hb, em_sb, fd)
    makedat(salt_bridges, hydrogen_bonds, em_hb, em_sb, trj, args.o)
    fd.close()

    if args.plot:
        plot_data(salt_bridges, hydrogen_bonds, em_hb, em_sb, trj, args.o, args.s)


if __name__ == "__main__":
    main()
