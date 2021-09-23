import argparse
import os
os.environ["SCHRODINGER_ALLOW_UNSAFE_MULTIPROCESSING"] = '1' #FUCK OFF

import matplotlib.pyplot as plt
import numpy as np
from schrodinger.application.desmond.packages import topo, traj, traj_util
from schrodinger.application.desmond.packages import analysis
from schrodinger.structutils import analyze
from schrodinger import structure
from tqdm import tqdm

# TODO: refactor using analyzers!


def main():
    parser = argparse.ArgumentParser(description='''
            Calculate Solvent Accessible Surface Area for a set of atoms over a trajectory.
            output: 
            sasa.dat: tot, hfo, hfi
            sasa_byres.dat: matrix of shape (times, residues)
            ''')

    parser.add_argument("cms", help="input cms file")
    parser.add_argument("asl", help="ASL specifying `atoms` for SASA calculations")
    parser.add_argument("-t", help="trajectory")
    parser.add_argument("-o", help="output filename. extensions are added automatically, default = sasa", default='sasa')
    parser.add_argument('-s', help='slice trj START:END:STEP (e.g.: "::10" will pick every 10th frame)')
    parser.add_argument('-cutoff', help='Atoms within this distance of `atoms` will be considered for occlusion. Requires `atoms` to be specified. (default: 8.0A)', default=8.0, type=float)
    parser.add_argument('-probe_radius', help="Probe radius, in units of angstroms, of the solvent. (default: 1.4A)", default=1.4, type=float)
    parser.add_argument('-resolution', help="Resolution to use. Decreasing this number will yield better results, increasing it will speed up the calculation.", default=0.2 , type=float)
    parser.add_argument('-exclude_water', help="If set to True then explicitly exclude waters in the method. This option is only works when 'atoms' argument is passed.", default=False, action='store_true')
    args = parser.parse_args()

    if not args.t:
        msys, cms, trj = traj_util.read_cms_and_traj(args.cms)
    else:
        msys, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.t)

    slicer = slice(*[int(x) if x else None for x in args.s.split(':')]) if args.s else slice(None)

    sel_idxs = cms.select_atom(args.asl)
    sel_st = cms.extract(sel_idxs)
    sel_res = structure.get_residues_by_connectivity(sel_st)  # maybe an overkill... st.residue should do the same?

    HYDROPHILIC = ['ASN ', 'CYS ', 'GLN ', 'SER ', 'THR ', 'ASP ', 'GLU ', 'ARG ', 'HIS ', 'LYS ']
    HYDROPHOBIC = ["ALA ", 'ILE ', 'LEU ', "MET",  'VAL ', 'PHE ', 'TRP ', 'GLY ', 'PRO ', 'TYR ']

    hydrophobic_idx = [i for i, res in enumerate(sel_res) if res.pdbres in HYDROPHOBIC]
    hydrophilic_idx = [i for i, res in enumerate(sel_res) if res.pdbres in HYDROPHILIC]

    tot_sasa = []
    hydrophobic_sasa = []
    hydrophilic_sasa = []
    byres_sasa = []
    for i, fr in enumerate(tqdm(trj[slicer])):  # type: ignore
        cms.setXYZ(fr.pos(cms.allaid_gids))
        sasa = analyze.calculate_sasa_by_residue(cms, atoms=sel_idxs, cutoff=args.cutoff, probe_radius=args.probe_radius, resolution=args.resolution, exclude_water=args.exclude_water)
        byres_sasa.append(sasa)
        sasa = np.asarray(sasa)
        hydrophilic_sasa.append(np.sum(sasa[hydrophilic_idx]))  # type: ignore
        hydrophobic_sasa.append(np.sum(sasa[hydrophobic_idx]))  # type: ignore
        tot_sasa.append(np.sum(sasa))

    byres_sasa = np.array(byres_sasa)

    with open(args.o + '.dat', 'w') as fh:
        fh.write('# time tot hydrophobic hydrophilic\n')
        for fr, tot, hfo, hfi in zip(trj[slicer], tot_sasa, hydrophobic_sasa, hydrophilic_sasa): # type: ignore
            fh.write(f'{fr.time} {tot} {hfo} {hfi}\n')

    with open(args.o + "_byres.dat", 'w') as fh:
        fh.write('# time ' + ' '.join(f"{r.chain}:{r.pdbres}_{r.resnum} " for r in sel_st.residue) + "\n")
        for fr, row in zip(trj[slicer], byres_sasa): #type: ignore
            fh.write(f'{fr.time} ' + " ".join(str(v) for v in row) + "\n")

    plt.plot(tot_sasa)
    plt.plot(hydrophilic_sasa)
    plt.plot(hydrophobic_sasa)
    plt.legend(['tot', 'hydrophilic', 'hydrophobic'])
    plt.xlabel('frame')
    plt.ylabel('SASA ($Å^2$)')
    plt.savefig(args.o + '.png')

    plt.figure()
    plt.imshow(byres_sasa.T, aspect='auto')
    plt.colorbar()
    plt.xlabel('frame')
    plt.ylabel('Residue SASA ($Å^2$)')
    plt.savefig(args.o + '_byres.png')


if __name__ == "__main__":
    main()
