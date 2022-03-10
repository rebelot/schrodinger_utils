from schrodinger.application.desmond.packages import topo, traj, traj_util, analysis
import sys
import argparse
import pickle
import numpy as np
import matplotlib.pyplot as plt

sys.setrecursionlimit(20000)


def gromos(rmsd, threshold, nmax):
    def gromos_inner(the_rmsd, the_threshold, n_max):
        # gromos clustering Daura et al. (Angew. Chem. Int. Ed. 1999, 38, pp 236-240)
        threshold_matrix = the_rmsd < the_threshold
        number_of_elements = np.count_nonzero(threshold_matrix, axis=0)
        the_medoid_index = np.argmax(number_of_elements)
        members = threshold_matrix[the_medoid_index, :]
        not_members = np.invert(members)
        the_member_indices = members.nonzero()[0]
        not_members_indices = not_members.nonzero()[0]
        result = [(the_medoid_index, the_member_indices)]
        if n_max > 1 and len(not_members_indices) > 0:
            other_results = gromos_inner(the_rmsd[not_members, :][:, not_members], the_threshold, n_max - 1)
            for o in other_results:
                result.append((not_members_indices[o[0]], not_members_indices[o[1]]))
        return result
    dimension = rmsd.shape[0]
    unassigned = set(range(dimension))
    the_result = gromos_inner(rmsd, threshold, nmax)
    assigned = {}
    for medoid_index, member_indices in the_result:
        assigned[medoid_index] = list(member_indices)
        unassigned.difference_update(set(member_indices))
    return assigned, list(unassigned)


def ucluster(rmsd, threshold, nmax=None):
    def find_medoids(the_rmsd, the_threshold, n_max):
        threshold_matrix = the_rmsd < the_threshold
        number_of_elements = np.count_nonzero(threshold_matrix, axis=0)
        medoid_index = np.argmax(number_of_elements)
        members = threshold_matrix[medoid_index, :]
        not_members = np.invert(members)
        not_members_indices = not_members.nonzero()[0]
        the_medoids = [medoid_index]
        if n_max > 1 and len(not_members_indices) > 0:
            other_medoids = find_medoids(the_rmsd[not_members, :][:, not_members], the_threshold, n_max - 1)
            for o in other_medoids:
                the_medoids.append(not_members_indices[o])
        return the_medoids

    dimension = rmsd.shape[0]
    medoids = find_medoids(rmsd, threshold, nmax)
    frames = list(range(dimension))
    assigned = {medoid: [] for medoid in medoids}
    unassigned = []
    for frame in frames:
        frames_distances = rmsd[:, frame]
        medoids_distances = frames_distances[medoids]
        closest_medoid = np.argmin(medoids_distances)
        if medoids_distances[closest_medoid] <= threshold:
            assigned[medoids[closest_medoid]].append(frame)
        else:
            unassigned.append(frame)
    return assigned, unassigned


def compute_nonlocality(assigned, unassigned):
    result = 0
    for _, cluster in assigned.items():
        arr = np.array(sorted(cluster))
        diff = np.diff(arr,1)
        result += np.sum(diff**2)
    result -= len(unassigned)**2
    result -= len(assigned)**3
    return result


def find_threshold(rmsd, nmax, method='closest'):
    dimension = rmsd.shape[0]
    if method == 'closest':
        find_clusters = ucluster
    else:
        find_clusters = gromos
    raw_thresholds = np.linspace(rmsd.max() * 0.01, rmsd.max() * 0.99, endpoint=True, )[1:]
    raw_clusters_list = [find_clusters(rmsd, threshold, nmax) for threshold in raw_thresholds]
    raw_nonlocality = [compute_nonlocality(assigned, unassigned) for (assigned, unassigned) in raw_clusters_list]
    first = 0
    last = 0
    for i, non_loc in enumerate(raw_nonlocality):
        if non_loc <= 0 and first == 0:
            continue

        if first == 0:
            first = i
            continue
        if non_loc > dimension and last == 0:
            continue
        if non_loc <= dimension :
            last = i
            break
    thresholds = np.linspace(raw_thresholds[first], raw_thresholds[last], endpoint=True, num=200)
    clusters_list = [find_clusters(rmsd, threshold, nmax) for threshold in thresholds]
    nonlocality = [compute_nonlocality(assigned, unassigned) for (assigned, unassigned) in clusters_list]
    selected = np.array(nonlocality).argmax()
    assigned, unassigned = clusters_list[selected]
    threshold = thresholds[selected]
    result = {'threshold': threshold,
              'clusters': assigned,
              'unassigned': unassigned,
              'method': method,
              'thresholds': thresholds,
              'nonlocality': nonlocality,
              'selected': selected}
    return result


def find_cluster(rmsd, nmax, threshold, method='closest'):
    if method == 'closest':
        find_clusters = ucluster
    else:
        find_clusters = gromos
    assigned, unassigned = find_clusters(rmsd, threshold, nmax)
    result = {'threshold': threshold,
              'clusters': assigned,
              'unassigned': unassigned,
              'method': method}
    return result


def write_log(outname, result, trajectory):
    unassigned = result['unassigned']
    assigned = sort_clusters(result['clusters'])
    log = open(f'{outname}_cluster.log', "w")
    log.write(f'Threshold : {result["threshold"]}, ')
    if 'selected' in result.keys():
        log.write('computed.\n')
    else:
        log.write('assigned.\n')
    for index, item in enumerate(assigned):
        log.write(f'Cluster # {index + 1:2d}, medoid is frame {trajectory[item[0]].orig_index:5d}, '
                  f'number of frames {len(item[1])}\n')
    if len(unassigned) > 0:
        log.write(f'{len(unassigned)} frames were further than {result["threshold"]}'
                  ' from each medoid and were unassigned\n')
    if 'selected' in result.keys():
        for th, ncl in zip(result['thresholds'], result['nonlocality']):
            log.write(f'Threshold {th:7.3f} has {ncl:6d} nonlocality.\n')
        log.write('\n')
    for index, item in enumerate(assigned):
        log.write(f'Cluster # {index + 1:2d}, medoid is frame {trajectory[item[0]].orig_index:5d}\n')
        for i, frame in enumerate(item[1]):
            log.write(f'{trajectory[frame].orig_index}')
            if (i + 1) % 16 == 0:
                log.write('\n')
            else:
                log.write(', ')
        log.write('\n\n')


def plot(outname, result, trajectory):
    unassigned = result['unassigned']
    assigned = sort_clusters(result['clusters'])
    xx = [trajectory[frame].time / 1000.0 for frame in unassigned]  # nanoseconds
    yy = np.full((len(xx),), 0, dtype=float)
    fig, ax1 = plt.subplots()
    ax1.scatter(xx, yy, alpha=0.2, color='red')
    if xx:
        labels = ['None']
        labels2 = [str(len(xx))]
    else:
        labels = ['']
        labels2 = ['']
    ticks = [0.0]
    for i, (medoid, cluster) in enumerate(assigned):
        if i > 16:
            continue
        labels.append(str(trajectory[int(medoid)].orig_index))
        labels2.append(str(len(cluster)))
        ticks.append(float(i+1))
        x = [trajectory[frame].time / 1000.0 for frame in cluster]  # nanoseconds
        y = np.full((len(x),), i+1, dtype=float)
        ax1.scatter(x, y, alpha=0.2, color='black')
    ax1.set_xlabel("Time (ns)")
    ax1.set_ylabel("Cluster medoid")
    plt.title(f'Threshold = {result["threshold"]}, method = {result["method"]}')
    ax1.set_yticks(ticks)
    ax1.set_yticklabels(labels)
    if xx:
        ax1.set_ylim(-1, len(labels))
    else:
        ax1.set_ylim(0, len(labels))
    ax2 = ax1.twinx()
    ax2.set_yticks(ticks)
    ax2.set_yticklabels(labels2)
    if xx:
        ax2.set_ylim(-1, len(labels))
    else:
        ax2.set_ylim(0, len(labels))
    ax2.set_ylabel("Number of frames")
    plt.savefig(f'{outname}_cluster.png')
    if 'selected' in result.keys():
        plt.figure()
        x = np.array(result['thresholds'], dtype=float)
        y = np.array(result['nonlocality'], dtype=int)
        plt.scatter(x, y)
        x = x[result['selected']]
        y = y[result['selected']]
        plt.scatter(x, y, color='red')
        plt.xlabel("Threshold (Angstroms)")
        plt.ylabel("Nonlocality")
        plt.title(f'Threshold = {result["threshold"]}, method = {result["method"]}')
        plt.savefig(f'{outname}_threshold.png')


def sort_clusters(cluster_dict):
    cluster_list = [(medoid, members) for medoid, members in cluster_dict.items()]
    return sorted(cluster_list, key=lambda x: len(x[1]), reverse=True)


def read_traj(cmsfile, my_slicer):
    the_msys, the_cms, trajectory = traj_util.read_cms_and_traj(cmsfile)
    trajectory = trajectory[my_slicer]
    the_time = [frame.time for frame in trajectory]
    return the_time, the_msys, the_cms, trajectory


def write_clusters(outbase, result, the_trj, the_cms, split):
    clusters = sort_clusters(result['clusters'])
    nptrj = np.array(the_trj)
    for i, (med, mem) in enumerate(clusters):
        the_cms.setXYZ(nptrj[med].pos())
        if split:
            traj.write_traj(nptrj[mem], f'{outbase}_cluster_{i + 1}_trj')
        the_cms.fix_filenames(f'{outbase}_cluster_{i + 1}-out.cms', f'{outbase}_cluster_{i + 1}_trj')
        the_cms.write(f'{outbase}_cluster_{i + 1}-out.cms')


if __name__ == "__main__":
    args = argparse.ArgumentParser(description="Perform RMSD-based clustering of a Desmond MD with gromos algorithm.")
    args.add_argument('cms', help="Cms input file")
    args.add_argument('output', help="Output base name")
    args.add_argument('-rmsd_asl', help="Define ASL for RMSD calculation", default='a.pt CA')
    args.add_argument('-fit_asl', help="Define ASL for fitting;", default=None)
    args.add_argument('-t', help="RMSD threshold defining clusters (Angstroms)", type=float)
    args.add_argument('-n', help="Max number of clusters returned", default=None, type=int)
    args.add_argument('-w', '--write_matrix', help="Store intermediate distance matrix to FILE", default=None)
    args.add_argument('-l', '--load_matrix', help="Load intermediate distance matrix from FILE", default=None)
    args.add_argument('-s', '--slice', help='slice trajectory START:END:STEP')
    args.add_argument('-split_trj', help='save cluster frames in trajectories', action='store_true')
    args.add_argument('-m', '--method', choices=['closest', 'largest'],
                      help='crietrion to assign frames to medoids. Default "closest" (ucluster)', default='closest')
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
        print('loading trajectory')
        time, msys, cms, trj = read_traj(args.cms, slicer)
    else:
        print('loading trajectory')
        time, msys, cms, trj = read_traj(args.cms, slicer)
        print('generating distance matrix (RMSD)')
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

    if args.n is None:
        args.n = rmsd_matrix.shape[0]//2

    print("clustering")
    if args.t:
        my_result = find_cluster(rmsd_matrix, args.n, args.t, args.method)
    else:
        my_result = find_threshold(rmsd_matrix, args.n, args.method)
    print("saving outputs")
    write_log(args.output, my_result, trj)
    plot(args.output, my_result, trj)
    if args.split_trj:
        write_clusters(args.output, my_result, trj, cms, args.split_trj)
