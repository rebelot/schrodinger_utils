#!/opt/schrodinger/suites2019-3/run

import matplotlib

matplotlib.use("Agg")
from schrodinger.application.desmond.packages import traj_util, topo, traj, analysis
from schrodinger.structutils import analyze
from schrodinger.structure import StructureReader
import argparse
import numpy as np
import matplotlib.pyplot as plt
import sys


def main():
    parser = argparse.ArgumentParser(
        description="Calculate RMSD of selected atoms over the course of MD simulation."
    )
    parser.add_argument("cms", help="Input cms file.")
    parser.add_argument("-trj", help="Input trajectory dir", default=None)
    parser.add_argument("-mode", help="RMSD or RMSF, default RMSD", default="RMSD")
    parser.add_argument("-s", help="Slice trajectory START:END:STEP")
    parser.add_argument(
        "-rmsd",
        help="Atom selection for RMSD calculation, default 'a.pt CA'",
        default="a.pt CA",
    )
    parser.add_argument(
        "-fit",
        help="Atom selection for superposition, defaults to 'rmsd' selection",
        default=None,
    )
    parser.add_argument(
        "-ref",
        help="Reference frame (default 0) OR structure filename",
        default="0",
        type=str,
    )
    parser.add_argument("-o", help="Output filename, defaults to stdout.")
    parser.add_argument("-p", help="Plot results", action="store_true")
    args = parser.parse_args()

    if not args.trj:
        msys, cms, trj = traj_util.read_cms_and_traj(args.cms)
    else:
        msys, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.trj)

    if args.s:
        start, end, step = args.s.split(":")
        start = int(start) if start else None
        end = int(end) if end else None
        step = int(step) if step else None
        slicer = slice(start, end, step)
    else:
        slicer = slice(None)

    rmsd_asl = args.rmsd
    rmsd_aids = cms.select_atom(rmsd_asl)
    rmsd_Atoms = list(analyze.get_atoms_from_asl(cms, rmsd_asl))
    rmsd_gids = topo.aids2gids(cms, rmsd_aids, include_pseudoatoms=False)

    fit_asl = args.fit or rmsd_asl
    fit_aids = cms.select_atom(fit_asl)
    fit_gids = topo.aids2gids(cms, fit_aids, include_pseudoatoms=False)

    if args.ref.isdigit():
        ref = int(args.ref)
        rmsd_ref_pos = trj[ref].pos(rmsd_gids)
        fit_ref_pos = trj[ref].pos(fit_gids)
    else:
        ref_st = StructureReader.read(args.ref)
        rmsd_ref_pos = np.array([ref_st.atom[i].xyz for i in analysis.evaluate_asl(ref_st, rmsd_asl)])
        fit_ref_pos = np.array([ref_st.atom[i].xyz for i in analysis.evaluate_asl(ref_st, fit_asl)])

    mode = args.mode.lower()
    if mode == "rmsd":
        analyzer = analysis.RMSD(
            msys,
            cms,
            rmsd_aids,
            rmsd_ref_pos,
            fit_aids=fit_aids,
            fit_ref_pos=fit_ref_pos,
        )
    elif mode == "rmsf":
        analyzer = analysis.RMSF(
            msys, cms, rmsd_aids, fit_aids=fit_aids, fit_ref_pos=fit_ref_pos
        )
    else:
        raise ValueError("Unrecognized mode, specify one of (RMSD, RMSF)")

    res = analysis.analyze(trj[slicer], analyzer)

    fh = open(args.o + ".dat", "w") if args.o else sys.stdout

    if mode == "rmsd":
        for fr, r in zip(trj[slicer], res):
            fh.write(f"{fr.time/1000.0} {r}\n")
    else:
        for a, r in zip(rmsd_Atoms, res):
            fh.write(f"{a.index} {a.resnum} {r}\n")

    if args.p:
        out = args.o + ".png" if args.o else f"{mode}_calc.png"
        if mode == "rmsd":
            plt.plot([fr.time / 1000 for fr in trj[slicer]], res)
            plt.xlabel("time (ns)")
        else:  # assume RMSF
            plt.plot(res)
            n = np.linspace(0, len(res), 10)
            plt.xticks(
                n,
                [a.resnum for a in rmsd_Atoms[:: int(len(res) / 10)]],
                rotation="vertical",
            )
            plt.xlabel("atom index")
        plt.ylabel(f"{mode.upper()} (Å)")
        plt.savefig(out)


if __name__ == "__main__":
    main()
