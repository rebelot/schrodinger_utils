import argparse

import matplotlib.pyplot as plt
import numpy as np
from schrodinger.application.desmond.packages import topo, traj, traj_util
from schrodinger.structutils import analyze
from tqdm import tqdm

#TODO: refactor using analyzers!

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("cms", help="input cms file")
    parser.add_argument("-t", help="trajectory", required=True)
    parser.add_argument("atoms", help="ASL specifying atoms for SASA calculations")
    parser.add_argument("-o", help="output filename. extensions are added automatically")
    parser.add_argument('-s', help='slice trj START:END:STEP (e.g.: "::10" will pick every 10th frame)')
    args = parser.parse_args()

    msys, cms = topo.read_cms(args.cms)
    trj = traj.read_traj(args.t)

    start, end, step = args.s.split(':')
    start = int(start) if start else None
    end = int(end) if end else None
    step = int(step) if step else None
    slicer = slice(start, end, step)

    subsys = cms.select_atom(args.atoms)
    subst = cms.extract(subsys)
    subres = list(subst.residue)

    HYDROPHILIC = ["ARG ", "ASP ", "GLU ", "HIS ", "ASN ", "GLN ", "LYS ", "SER ", "THR "]
    HYDROPHOBIC = ["PHE ", "LEU ", "ILE ", "TYR ", "TRP ", "VAL ", "MET ", "PRO ", "CYS ", "ALA "]

    hydrophobic_idx = [i for i, res in enumerate(subres) if res.pdbres in HYDROPHOBIC]
    hydrophilic_idx = [i for i, res in enumerate(subres) if res.pdbres in HYDROPHILIC]

    tot_sasa = []
    hydrophobic_sasa = []
    hydrophilic_sasa = []
    byres_sasa = []
    for fr in tqdm(trj[slicer]):
        cms.setXYZ(fr.pos())
        sasa = analyze.calculate_sasa_by_residue(cms, atoms=subsys, exclude_water=True)
        byres_sasa.append(sasa)
        sasa = np.array(sasa)
        hydrophilic_sasa.append(np.sum(sasa[hydrophilic_idx]))
        hydrophobic_sasa.append(np.sum(sasa[hydrophobic_idx]))
        tot_sasa.append(np.sum(sasa))

    byres_sasa = np.array(byres_sasa).T

    with open(args.o + '.dat', 'w') as fh:
        fh.write('# tot hydrophobic hydrophilic\n')
        for i in range(len(trj[slicer])):
            fh.write(f'{tot_sasa[i]} {hydrophobic_sasa[i]} {hydrophilic_sasa[i]}\n')

    byres_sasa.tofile(args.o + '_byres.dat', sep=' ')

    plt.plot(tot_sasa)
    plt.plot(hydrophilic_sasa)
    plt.plot(hydrophobic_sasa)
    plt.legend(['tot', 'hydrophilic', 'hydrophobic'])
    plt.xlabel('frame')
    plt.ylabel('SASA ($Å^2$)')
    plt.savefig(args.o + '.png')

    plt.figure()
    plt.imshow(byres_sasa, aspect='auto')
    plt.colorbar()
    plt.xlabel('frame')
    plt.ylabel('Residue SASA ($Å^2$)')
    plt.savefig(args.o + '_byres.png')


if __name__ == "__main__":
    main()
