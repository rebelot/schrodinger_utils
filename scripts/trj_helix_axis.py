from schrodinger.application.desmond.packages import traj, traj_util, topo
from schrodinger.structutils import analyze
import numpy as np
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('cms')
    parser.add_argument('--helix', help='ASL specifying helix residue range')
    args = parser.parse_args()

    msys, cms, trj = traj_util.read_cms_and_traj(args.cms)
    
    for fr in trj:
        cms.setXYZ(fr.pos())
        atoms = analyze.get_atoms_from_asl(cms, args.helix + ' and a.pt CA')
        coords = np.array([a.xyz for a in atoms if a.secondary_structure == 1])
        coords_mean = coords.mean(axis=0)
        uu, dd, vv = np.linalg.svd(coords - coords_mean)


