from schrodinger.application.desmond.packages import topo, traj, traj_util, analysis
import sys
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap


def gromos(rmsd, threshold, nmax):
    # gromos clustering Daura et al. (Angew. Chem. Int. Ed. 1999, 38, pp 236-240)
    threshold_matrix = rmsd < threshold
    number_of_elements = np.count_nonzero(threshold_matrix, axis=0)
    medoid_index = np.argmax(number_of_elements)
    members = threshold_matrix[medoid_index, :]
    not_members = np.invert(members)
    member_indices = members.nonzero()[0]
    not_members_indices = not_members.nonzero()[0]
    result = [(medoid_index, member_indices)]
    if nmax > 1 and len(members) > 1:
        other_results = gromos(rmsd[not_members, :][:, not_members], threshold, nmax-1)
        for o in other_results:
            result.append((not_members_indices[o[0]], not_members_indices[o[1]]))
    return result


def write_log(outname, result, threshold, trajectory):
    short = open(f'{outname}_cluster.txt', "w")
    long = open(f'{outname}_cluster.log', "w")
    short.write(f'Threshold : {threshold}\n')
    long.write(f'Threshold : {threshold}\n')
    for index, item in enumerate(result):
        short.write(f'Cluster # {index + 1:2d}, medoid is frame {trajectory[item[0]].orig_index:5d}, '
                    f'number of frames {len(item[1])}\n')
        long.write(f'Cluster # {index + 1:2d}, medoid is frame {trajectory[item[0]].orig_index:5d}\n')
        for i, frame in enumerate(item[1]):
            long.write(f'{trajectory[frame].orig_index}')
            if (i + 1) % 16 == 0:
                long.write('\n')
            else:
                long.write(', ')
        long.write('\n\n')


def plot(outname, result, trajectory):
    name = "Set2"  # 8 colors colormap
    cmap = get_cmap(name)
    colors = cmap.colors
    ax = plt.gca()
    ax.set_prop_cycle(color=colors)
    ticks = []
    for i, cluster in enumerate(result):
        if i > 7:  # only first 8 are plotted
            continue
        ticks.append(i+1)
        x = [trajectory[frame].time / 1000.0 for frame in cluster[1]]  # nanoseconds
        y = np.full((len(x),), i+1, dtype=float)
        plt.scatter(x, y, alpha=0.2)
    plt.xlabel("Time (ns)")
    plt.ylabel("Cluster #")
    plt.yticks(ticks)
    plt.savefig(f'{outname}_cluster.png')


def read_traj(cmsfile, my_slicer):
    the_msys, the_cms, trajectory = traj_util.read_cms_and_traj(cmsfile)
    trajectory = trajectory[my_slicer]
    the_time = [frame.time for frame in trajectory]
    return the_time, the_msys, the_cms, trajectory


def write_clusters(the_trj, the_cms, clusters, outbase):
    nptrj = np.array(the_trj)
    for i, (med, mem) in enumerate(clusters):
        the_cms.setXYZ(nptrj[med].pos())
        traj.write_traj(nptrj[mem], f'{outbase}_cluster_{i + 1}_trj')
        the_cms.fix_filenames(f'{outbase}_cluster_{i + 1}-out.cms', f'{outbase}_cluster_{i + 1}_trj')
        the_cms.write(f'{outbase}_cluster_{i + 1}-out.cms')


if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Perform RMSD-based clustering of a Desmond MD with gromos algorithm.")
    args.add_argument('-c', '--cms', help="Cms input file")
    args.add_argument('-o', '--output', help="Output base name", required=True)
    args.add_argument('-r', '--rmsd_asl', help="Define ASL for RMSD calculation", default='a.pt CA')
    args.add_argument('-f', '--fit_asl', help="Define ASL for fitting;", default=None)
    args.add_argument('-t', '--threshold', help="RMSD threshold defining clusters (Angstroms)", default=1.5)
    args.add_argument('-m', '--max_clusters', help="Max number of clusters returned", default=8, type=int)
    args.add_argument('-w', '--write_matrix', help="Store intermediate distance matrix to FILE", default=None)
    args.add_argument('-l', '--load_matrix', help="Load intermediate distance matrix from FILE", default=None)
    args.add_argument('-s', '--slice', help='slice trajectory START:END:STEP')
    args = args.parse_args()

    slicer = slice(None)

    if args.slice:
        start, end, step = args.slice.split(':')
        start = int(start) if start else None
        end = int(end) if end else None
        step = int(step) if step else None
        slicer = slice(start, end, step)

    if not args.fit_asl:
        args.fit_asl = args.rmsd_asl

    if args.load_matrix:
        print(f'loading distance matrix from {args.load_matrix}')
        with open(args.load_matrix, "rb") as pk:
            data = pickle.load(pk)
        rmsd_matrix = data['rmsd']
        if args.slice and slicer != data['slicer']:
            print("saved and requested slices differ")
            sys.exit(1)
        slicer = data['slicer']
        if args.cms and args.cms != data['cms']:
            print("saved and requested cms differ")
            sys.exit(1)
        args.cms = data['cms']
        time, msys, cms, trj = read_traj(args.cms, slicer)
    else:
        print('generating distance matrix (RMSD)')
        time, msys, cms, trj = read_traj(args.cms, slicer)
        rmsd_gids = topo.aids2gids(cms, cms.select_atom(args.rmsd_asl), include_pseudoatoms=False)
        fit_gids = topo.aids2gids(cms, cms.select_atom(args.fit_asl), include_pseudoatoms=False)
        rmsd_matrix = analysis.rmsd_matrix(msys, trj, rmsd_gids, fit_gids)
        if args.write_matrix:
            print(f'saving distance matrix to {args.write_matrix}')
            data = {'slicer': slicer,
                    'cms': args.cms,
                    'rmsd': rmsd_matrix}
            with open(args.write_matrix, "wb") as pk:
                pickle.dump(data, pk)

    my_result = gromos(rmsd_matrix, args.threshold, args.max_clusters)
    write_log(args.output, my_result, args.threshold, trj)
    plot(args.output, my_result, trj)
    write_clusters(trj, cms, my_result, args.output)
