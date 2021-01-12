from schrodinger.structutils import rmsd, transform, interactions, analyze, measure
from schrodinger.structure import StructureReader, Structure
import numpy as np
import argparse
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial.transform import Rotation


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input', help='input file')
    parser.add_argument('ligand_asl', help='ASL that defines the ligand')
    parser.add_argument('receptor_asl', help='ASL that defines the receptor')
    parser.add_argument('XL-in', help='input file cointaining XL pairs')
    args = parser.parse_args()

    system: Structure = next(StructureReader(filename=args.input))
    ligand = system.extract(analyze.evaluate_asl(system, args.ligand_asl))
    receptor = system.extract(analyze.evaluate_asl(system, args.receptor_asl))

    ligand_xl, receptor_xl = read_xl_list(args.XL_in)

    ligand_xl_idx = np.asarray([analyze.evaluate_asl(ligand, asl)[0] for asl in ligand_xl]) - 1
    receptor_xl_idx = np.asarray([analyze.evaluate_asl(receptor, asl)[0] for asl in receptor_xl]) - 1

    # Initialize coordinates
    transform.translate_center_to_origin(receptor)
    transform.translate_center_to_origin(ligand)

    # Get gridpoints (Nx3)
    grid = init_grid(receptor, ligand, spacing=5)

    # Iterate on each gridpoint
    best_d = np.inf
    d = 0
    refine_radius = 10
    while True:
        best_pose, d, best_T = iterate(grid, ligand, ligand_xl_idx, receptor, receptor_xl_idx)
        grid = refine_grid(best_T, r=refine_radius, spacing=5)
        if d >= best_d:
            break
        best_d = d
        refine_radius -= 2

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    plottamela(grid, best_T, best_pose, ligand_xl_idx, receptor, receptor_xl_idx, f'THE BEST {best_d}')


def iterate(grid, ligand, ligand_xl_idx, receptor, receptor_xl_idx):
    best_pose = None
    best_d = np.inf
    best_T = 0
    for T in grid:
        d, pose = dock(grid, T, ligand, ligand_xl_idx, receptor, receptor_xl_idx)
        plottamela(grid, T, ligand, ligand_xl_idx, receptor, receptor_xl_idx, f'RMSD = {d}')
        if d is None:
            continue
        if d < best_d:
            best_pose, best_d = pose.copy(), d
            best_T = T
    return best_pose, best_d, best_T


def dock(grid, grid_point, ligand, ligand_xl_idx, receptor, receptor_xl_idx):
    # Translate receptor to grid point
    transform.translate_center_to_origin(receptor, grid_point.tolist())

    lig_vec = ligand.getXYZ()[ligand_xl_idx]
    rec_vec = receptor.getXYZ()[receptor_xl_idx]

    # Calculate super transformation matrix
    # (a, b) | a = R @ b
    R, RMSD = Rotation.align_vectors(rec_vec, lig_vec)
    R = R.as_matrix()
    R = np.vstack((R, np.array([0, 0, 0])))
    R = np.hstack((R, np.array([0, 0, 0, 1])[:, None]))
    # Apply transformation (rotation only)
    transform.transform_structure(ligand, R)

    # check that no clashes occur
    if len(list(interactions.clash_iterator(ligand, struc2=receptor))) > 1:
        return None, None
    return RMSD, ligand


def init_grid(receptor: Structure, ligand: Structure, spacing=100):
    receptor_xyz = receptor.getXYZ()
    ligand_xyz = ligand.getXYZ()
    a = np.linspace(receptor_xyz[:, 0].min() + ligand_xyz[:, 0].min(),
                    receptor_xyz[:, 0].max() + ligand_xyz[:, 0].max(), spacing)
    b = np.linspace(receptor_xyz[:, 1].min() + ligand_xyz[:, 1].min(),
                    receptor_xyz[:, 1].max() + ligand_xyz[:, 1].max(), spacing)
    c = np.linspace(receptor_xyz[:, 2].min() + ligand_xyz[:, 2].min(),
                    receptor_xyz[:, 2].max() + ligand_xyz[:, 2].max(), spacing)
    return np.array(np.meshgrid(a, b, c)).T.reshape(-1, 3)


def refine_grid(pos, r=10, spacing=100):
    a = np.linspace(pos[0].min() - r, pos[0].max() + r, spacing)
    b = np.linspace(pos[1].min() - r, pos[1].max() + r, spacing)
    c = np.linspace(pos[2].min() - r, pos[2].max() + r, spacing)
    return np.array(np.meshgrid(a, b, c)).T.reshape(-1, 3)


def read_xl_list(filename):
    ligand_xl = []
    receptor_xl = []
    with open(filename, 'r') as fh:
        lines = fh.readlines()
        ligand_xl = [line.split(';')[0] for line in lines if line[0] != '#']
        receptor_xl = [line.split(';')[1].strip()
                       for line in lines if line[0] != '#']
    return ligand_xl, receptor_xl


def plottamela(grid, grid_point, ligand, ligand_xl_idx, receptor, receptor_xl_idx, d):
    plt.pause(.0000000000001)
    ax = plt.gca()
    ax.clear()
    ligand_xyz = ligand.getXYZ()
    receptor_xyz = receptor.getXYZ()
    lig_ca = ligand_xyz[np.array(analyze.evaluate_asl(ligand, 'a.pt CA')) - 1]
    rec_ca = receptor_xyz[np.array(analyze.evaluate_asl(receptor, 'a.pt CA')) - 1]
    rec_het = receptor_xyz[np.array(analyze.evaluate_asl(receptor, 'r.pt LIP, POPC'))- 1]
    ax.set_xlim([-100,100])
    ax.set_ylim([-100,100])
    ax.set_zlim([-100,100])
    # ax.set_xlim([-2,2])
    # ax.set_ylim([-2,2])
    # ax.set_zlim([-2,2])
    ax.set_title(f'{d}')
    # Plot grid points
    ax.plot3D(grid[:, 0], grid[:, 1], grid[:, 2], 'ko', markersize=.5)
    # Plot current grid point
    ax.plot3D([grid_point[0]], [grid_point[1]], [grid_point[2]], 'ro')
    # Plot ligand coordinates
    ax.plot3D(lig_ca[:, 0], lig_ca[:, 1], lig_ca[:, 2], 'b-', linewidth=2)
    # Plot receptor coordinates
    ax.plot3D(rec_ca[:, 0], rec_ca[:, 1], rec_ca[:, 2], 'r-', linewidth=2)
    ax.plot3D(rec_het[:, 0], rec_het[:, 1], rec_het[:, 2], 'o', color='orange', markersize=1)
    # Plot contact residues
    ax.plot3D(ligand_xyz[ligand_xl_idx][:, 0], ligand_xyz[ligand_xl_idx][:, 1], ligand_xyz[ligand_xl_idx][:, 2], 'bo', markersize=3)
    ax.plot3D(receptor_xyz[receptor_xl_idx][:, 0], receptor_xyz[receptor_xl_idx][:, 1], receptor_xyz[receptor_xl_idx][:, 2], 'ro', markersize=3)
    # Plot contact residue pairs distances
    for a1, a2 in zip(ligand_xl_idx, receptor_xl_idx):
        ax.plot3D([ligand_xyz[a1, 0], receptor_xyz[a2, 0]],
                  [ligand_xyz[a1, 1], receptor_xyz[a2, 1]],
                  [ligand_xyz[a1, 2], receptor_xyz[a2, 2]], 'g--')
    return ax


def test():
    # s1 = analyze.create_new_structure()
    # s1.addAtom('H', 2*.5,0,0)
    # s1.addAtom('H', 2*-.5,0,0)
    # s1.addAtom('H', 0,2*np.sqrt(3)/2,0)
    # s1.addAtom('H', 0,2*(np.sqrt(3)/2)/2,2*np.sqrt(2/3))
    # s1.addBonds([(1,2,1),(1,3,1),(1,4,1),(3,2,1),(3,4,1),(2,4,1)])
    # s2 = s1.copy()
    # ligand = s1
    # receptor = s2
    # ligand_xl_idx = [3]
    # receptor_xl_idx = [3]
    #
    system = next(StructureReader('./rHDL-LCAT_wombppdxa_00_MD-ea-out.cms'))
    ligand = system.extract(analyze.evaluate_asl(system, 'c.n C'))
    receptor = system.extract(analyze.evaluate_asl(system, 'c.n A, B'))
    ligand_xl, receptor_xl = read_xl_list('./xl_new.txt')
    ligand_xl_idx = [analyze.evaluate_asl(ligand, asl)[0] for asl in ligand_xl]
    receptor_xl_idx = [analyze.evaluate_asl( receptor, asl)[0] for asl in receptor_xl]
    transform.translate_center_to_origin(receptor)
    transform.translate_center_to_origin(ligand)
    grid = init_grid(receptor, ligand, spacing=10)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot3D(grid[:, 0], grid[:, 1], grid[:, 2], 'ko', markersize=.5)
    # Iterate on each gridpoint
    best_d = np.inf
    d = 0
    refine_radius = 10
    while True:
        best_pose, d, best_T = iterate(grid, ligand, ligand_xl_idx, receptor, receptor_xl_idx)
        grid = refine_grid(best_T, r=refine_radius, spacing=10)
        if d >= best_d:
            print('Converged')
            break
        best_d = d
        refine_radius -= 2
    else:
        print('Grid collapse')
    plottamela(grid, best_T, best_pose, ligand_xl_idx, receptor, receptor_xl_idx, f'THE BEST {best_d}')
