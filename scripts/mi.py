import argparse
import multiprocessing as mp
import time

import matplotlib.pyplot as plt
import numpy as np
from sklearn.feature_selection import mutual_info_regression


def worker(x):  # -> (na * d)[:i]
    return mutual_info_regression(
        x, x[:, -1], discrete_features=False, n_neighbors=6, random_state=42
    )


def get_MI(x, njobs, dump=False, out="corr"):
    nt, na, d = x.shape
    x = x - x.mean(axis=0)  # center atoms on their mean position
    # normalize fluctuation interval between 0, 1 WARNING: is this required?
    x = (x - x.min(axis=0)) / np.ptp(x, axis=0)
    x = x.reshape(nt, na * d)
    MI = np.zeros((na * d, na * d))

    start_time = time.time()
    with mp.Pool(processes=njobs) as pool:
        res = pool.imap(worker, (x[:, : i + 1] for i in range(na * d)))
        for i, r in enumerate(res):
            MI[i, : i + 1] = r
            MI[: i + 1, i] = r
    print(
        f"finished {na} atoms for {nt} frames in", time.time() - start_time, "seconds"
    )

    if dump:
        import pickle

        with open(out + ".pickle", "wb") as f:
            pickle.dump(MI, f)

    # postprocess
    # (ai, xi, aj, xj)
    MI = MI.reshape(na, d, na, d).mean(axis=(1, 3))
    MI = np.sqrt(1 - np.exp(-MI))
    MI /= MI.max()
    return MI


def get_cormat(x):
    """
    Get covariance matrix of atom fluctuations
    :param x: atom positions over time
    :type  x: np.array of shape (time, atoms, dims)
    :return : the normalized covariance matrix
    :return type: np.array (atoms, atoms)
    """
    # dots = np.zeros((na, na))
    # for i in range(nt):
    #     dot = np.dot(fluct[i], fluct[i].T)
    #     dots = dots + dot
    fluct = x - np.mean(x, axis=0)
    dots = 1 / x.shape[0] * np.einsum("ijk,ilk->jl", fluct, fluct)

    diagonal = np.diag(dots)
    norm_matrix = np.outer(diagonal, diagonal)
    norm_matrix = np.sqrt(np.absolute(norm_matrix))
    corr_matrix = np.divide(dots, norm_matrix)
    return corr_matrix


def preproc_schrodinger(args):
    from schrodinger.application.desmond.packages import topo, traj, traj_util

    align_sel = "a.pt CA" if args.align == "CA" else args.align
    corr_sel = "a.pt CA" if args.asl == "CA" else args.asl

    if args.t:
        msys, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.t)
    else:
        msys, cms, trj = traj_util.read_cms_and_traj(args.cms)

    slicer = (
        [int(i) if i else None for i in args.s.split(":")]
        if args.s
        else slice(None, None)
    )

    trj = trj[slicer]

    if args.align:
        print(f'Aligning trajectory to "{align_sel}"')
        fit_gids = topo.asl2gids(cms, align_sel)
        ref = cms.extract(cms.select_atom(align_sel)).getXYZ()
        trj = topo.superimpose(msys, fit_gids, trj, ref)

    gids = topo.asl2gids(cms, corr_sel)
    x = np.array([fr.pos(gids) for fr in trj])  # type: ignore
    return x


def preproc_mda(args):
    import MDAnalysis as mda
    from MDAnalysis.analysis import align

    align_sel = "name CA" if args.align == "CA" else args.align
    corr_sel = "name CA" if args.asl == "CA" else args.asl

    top = args.cms
    traj = (
        args.t if isinstance(args.t, list) else [args.t]
    )  # better way to handle this using argparse action='append'
    slicer = (
        [int(i) if i else None for i in args.s.split(":")]
        if args.s
        else slice(None, None)
    )

    U = mda.Universe(top, *traj)

    if args.align:
        print(f'Aligning trajectory to "{align_sel}"')
        U.trajectory[0]
        mobile = U.copy()
        align.AlignTraj(mobile, U, select=align_sel, in_memory=True).run()
        U = mobile

    U.trajectory[0]
    atoms = U.select_atoms(corr_sel)
    x = np.array([atoms.positions for _ in U.trajectory[slicer]])
    return x


def main():
    parser = argparse.ArgumentParser(
        description="""
            Calculate generalized correlations based on mutual information of selected atoms fluctuations
            example:
            mi.py my_md.pdb corr_out -t my_md.xtc -backend mda -j 16 -s ::10 -align "name CA and segid A" -asl "name CA"
            """
    )
    parser.add_argument("cms", help="input topology file")
    parser.add_argument("out", help="base name for the output")
    parser.add_argument(
        "-backend",
        help="preprocess trajectory using selected backend: MDAanalysis (mda) or Schrodinger (schr)",
        choices=["mda", "schr"],
        default="schr",
    )
    parser.add_argument(
        "-t",
        help="trj dir, this is mandatory if backend is MDAnalysis",
        metavar="trajectory",
    )
    parser.add_argument(
        "-j",
        help="number of processes (max physical cores), default 1",
        type=int,
        metavar="N cores",
        default=1,
    )
    parser.add_argument("-asl", help="atom selection, default CA", default="CA")
    parser.add_argument("-s", help="slicer", metavar="START:END:STEP")
    parser.add_argument(
        "-pickle", help="dump unmodified MI matrix to pickle file", action="store_true"
    )
    parser.add_argument(
        "-corr", help="calculate the correlation matrix", action="store_true"
    )
    parser.add_argument(
        "-align",
        help="asl to align trajectory, default CA",
        const="CA",
        nargs="?",
        metavar="ASL",
    )
    args = parser.parse_args()

    if args.backend == "schr":
        preproc = preproc_schrodinger
    else:
        preproc = preproc_mda

    x = preproc(args)
    nt, na, d = x.shape

    MI = get_MI(x, args.j, dump=args.pickle, out=args.out)

    with open(args.out + "_MI.dat", "w") as f:
        f.write(f"# MI matrix, {na} x {na}\n")
        for row in MI:
            f.write(" ".join(str(i) for i in row))
            f.write("\n")

    fig = plt.figure(figsize=(15, 15))
    fig.suptitle("Mutual Information")
    plt.imshow(MI, origin="lower")
    plt.colorbar()
    fig.tight_layout()
    plt.savefig(args.out + "_MI.png")

    if args.corr:
        C = get_cormat(x.reshape(nt, na, d))
        with open(args.out + "_Cor.dat", "w") as f:
            f.write(f"# Correlation matrix, {na} x {na}\n")
            for row in C:
                f.write(" ".join(str(i) for i in row))
                f.write("\n")
        fig = plt.figure(figsize=(15, 15))
        fig.suptitle("Correlation Matrix")
        plt.imshow(C, origin="lower")
        plt.colorbar()
        fig.tight_layout()
        plt.savefig(args.out + "_Cor.png")


if __name__ == "__main__":
    main()
