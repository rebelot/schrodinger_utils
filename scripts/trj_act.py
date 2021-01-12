from schrodinger.application.desmond.packages import traj_util, traj, topo
from schrodinger.application.desmond.packages.msys import pfx as pfx_module
from tqdm import tqdm
import typing
import numpy as np
import argparse
from argparse import RawDescriptionHelpFormatter


def main():
    parser = argparse.ArgumentParser(
        description="""
Align, Center and Translate stuff within a simulation box using expressions.
e.g.: to get something like the following, use -e "protein:center:0,1" -e "membrane:center:2" -off ::0.5
Get this fancy visualization of stuff interacting with membranes:
 - the protein is centered on the XY plane
 - the membrane is centered on the Z axis
 - atoms are then shifted upwards 1/2 of the length of the Z axis

        |||||||||||||||||||
        ooooooooooooooooooo

             ~protein~

        oooooooooooooooooooo
        ||||||||||||||||||||""",
        formatter_class=RawDescriptionHelpFormatter,
    )
    parser.add_argument("cms", help="Input cms file")
    parser.add_argument("trj", help="Input trajectory dir")
    parser.add_argument("out", help="output basename")
    parser.add_argument( "-e", action="append",
        help="""Expression in the form of ASL:MODE[:DIMS or :REF], where:
ASL is the asl expression specifying the selection to acts upon; 
MODE = center or align; DIMS is a comma-separated list of dimension (X = 0, Y = 1, Z = 2) for centering, default is 0,1,2;
REF is the index of the frame used as reference positions for alignment, default is 0.
Multiple expressions can be provided and will be executed sequentially from left to right.
Note that alignment will also automatically perform centering.
e.g: -e "protein:center:0,1" will center the protein on the XY plane;
-e "membrane:center:2" will center the membrane on the Z axis;
-e "protein and a.pt CA:align:100" will align protein Calpha using frame 100 as reference.""")
    parser.add_argument( "-off",
        help='Translate atoms by "fa:fb:fc", where fa, fb, fc are fractions of the simulation box axes; \
            e.g.: -off "::0.5" will translate atoms upwards by half of the simulation box "c" axis', default=None)
    args = parser.parse_args()

    msys, cms = topo.read_cms(args.cms)
    trj: typing.List[traj.Frame] = traj.read_traj(args.trj)

    pfx = pfx_module.Pfx(msys.glued_topology, fixbonds=True)

    for expr in args.e:
        print("Executing", expr)
        expr = expr.split(":")
        opt = None
        if len(expr) == 2:
            asl, mode = expr
        else:
            asl, mode, opt = expr
        gids = topo.aids2gids(cms, cms.select_atom(asl))
        if mode == "align":
            ref_frame = int(opt) if opt else 0
            ref_pos = trj[ref_frame].pos(gids)
            weights = [msys.atom(g).mass for g in gids]
            topo.superimpose(msys, gids, trj, ref_pos, weights)
        elif mode == "center":
            dims = [int(d) for d in opt.split(",")] if opt else [0, 1, 2]
            topo.center(msys, gids, trj, dims=dims)
        else:
            raise ValueError("Unrecognized mode")

    if args.off:
        print("Applying translation offset", args.off)
        foff = np.array([float(i) if i else 0 for i in args.off.split(":")])
        for fr in trj:
            offset = fr.box.diagonal() * foff
            fr.moveby(*offset)
            pfx.apply(
                fr.pos(), fr.box, fr.vel()
            )  # <- Kung Fu. This also makes the system whole (see topo.make_whole docstring)

    print("Writing output...")
    cms_updated = topo.update_cms(cms, trj[-1])
    cms_updated.fix_filenames(args.out + '-out.cms', args.out + '_trj')
    traj.write_traj(trj, args.out + "_trj")
    cms_updated.write(args.out + "-out.cms")
    # cms_updated.property["s_chorus_trajectory_file"] = args.out + "_trj"


if __name__ == "__main__":
    main()
