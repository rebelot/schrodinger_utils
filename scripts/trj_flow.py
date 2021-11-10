"""
File: trj_flow.py
Author: Tommaso laurenzi
Email: tommaso.laurenzi@unimi.it
Github: https://github.com/rebelot
Description: Compute and visualize any molecule flow within Molecular Dynamic simulations. 
Date: 2021-11-01
"""

import argparse
import multiprocessing as mp
from typing import List
import pickle

import numpy as np
from matplotlib.colors import to_rgb
from schrodinger.application.desmond.packages import topo, traj, traj_util
from schrodinger.structutils import analyze
from scipy.spatial import ConvexHull, Delaunay
from scipy import interpolate
from tqdm import tqdm
from scipy.signal import savgol_filter
from sklearn.cluster import DBSCAN

# https://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl
# TODO:
# - [ ] FIXME: extra care in handling empty paths
# - [x] : store/load data
# - [ ] generate some statistics log (i.e n molecules, n paths, mean elapsed time in site)
# - [ ] mdanalysis 
# - [ ] maybe add support for variable n_atom_per_mol (multiple obj types)
# - [x] clustering
# - [-] path smoothing
# - [ ] multiprocessing
# - [x] draw volumes FEAT: make color ramp

class Cgo:
# http://pymol.sourceforge.net/newman/user/S0500cgo.html
    COLOR     = 6.0
    BEGIN     = 2.0
    END       = 3.0
    LINE_LOOP = 2.0
    SPHERE    = 7.0
    VERTEX    = 4.0
    LINEWIDTH = 10.0

    def __init__(self):
        self.obj = []

    def scatter(self,
        points,
        linecolor=[1.0, 1.0, 1.0],
        line=True,
        marker=False,
        radius=0.1,
        markercolor=False,
        linewidth=1.0,
    ):
        obj = []
        if line:
            obj.extend([
                    self.BEGIN, self.LINE_LOOP,
                    self.COLOR, *linecolor,
                    self.LINEWIDTH, linewidth,
                ])
            for xyz in points:
                obj.extend([self.VERTEX, *xyz])
            obj.extend([self.END])

        if marker:
            obj.extend([self.COLOR, *(markercolor or linecolor)])
            for xyz in points:
                obj.extend([self.SPHERE, *xyz, radius])

        self.obj.extend([float(i) for i in obj])

    @classmethod
    def Scatter(cls, *args, **kwargs):
        cgo = Cgo()
        cgo.scatter(*args, **kwargs)
        return cgo


class PathCollection:
    def __init__(self, paths, interpf=None, smoothf=0, splorder=3):
        self.paths = paths
        self.interpf = interpf
        self.smoothf = smoothf
        self.splorder = splorder

    @property
    def smooth(self):
        return self.interpf is not None

    def _collect(self, attr):
        return [pos for path in self.paths for pos in getattr(path, attr) if len(pos)]

    def _savgol(self, pos, win=5, poly=2):
        try:
            return savgol_filter(pos, axis=0, window_length=win, polyorder=poly)
        except (TypeError, ValueError):
            return pos

    def _interpsmooth(self, pos):
        N = len(pos) * self.interpf # type: ignore maybe set a maximum value..
        try:
            tck, u = interpolate.splprep(pos.T.tolist(), s=self.smoothf, k=self.splorder)
            xyz = interpolate.splev(np.linspace(u.min(), u.max(), N), tck)
            return np.array(xyz).T
        except TypeError:
            return pos

    def _smoothdec(func): #type: ignore
        def wrapper(self, *args, **kwargs):
            res = func(self, *args, **kwargs) # type: ignore
            if not self.smooth:
                return res

            new_res = []
            for pos in res:
                new_res.append(self._interpsmooth(pos))
            return new_res
        return wrapper

    @property
    def inlets(self):
        items = self._collect("inlets")
        return np.vstack(items) if items else []

    @property
    def outlets(self):
        items = self._collect("outlets")
        return np.vstack(items) if items else []

    @property
    @_smoothdec # type: ignore
    def pos(self):
        return self._collect("pos")

    @property
    @_smoothdec # type: ignore
    def incoming(self):
        return self._collect("incoming")

    @property
    @_smoothdec # type: ignore
    def outgoing(self):
        return self._collect("outgoing")

    @property
    @_smoothdec # type: ignore
    def collat(self):
        return self._collect("collat")

    @property
    @_smoothdec # type: ignore
    def inside(self):
        return self._collect("inside")



class SubPath:
    r"""
    This class stores a subpath, defined as the contiguous set of points for
    which object is in scope and it has passed at least once from site.

    - A SubPath always begins in scope, either at a boundary or at any point (if it is already there)
    - A SubPath always ends in scope, either at a boundary or at any point
      (if it is still there at the end of the simulation)
    - A SubPath has at least one point within site and scope
    
        -----------------------
        |                   .2|...    Path      = 1 to 4
        |        _____   ../  |   \   SubPath_1 = 1 to 2
        |   ....|.....|./   ..|.../   SubPath_2 = 3 to 4
    0...|1./  4.|..   |    /  |
        |       |__\__|    :  |       From 2 to 3 object enters scope again, but
        |           :      :  |       it does not pass from site, thus this not a
        |           \.3    \..|..     valid subpath.
        ---------------\-------  \
                        \......../
    """
    def __init__(self, pos, in_site):
        self.pos = pos
        self.in_site = in_site
        self.inside = pos[in_site]
        inlets_i = np.where(np.diff(in_site.astype(int)) == 1)[0] + 1
        outlets_i = np.where(np.diff(in_site.astype(int)) == -1)[0] + 1
        self.inlets = pos[inlets_i]
        self.outlets = pos[outlets_i]
        self.incoming = pos[: inlets_i[0] + 1] if len(inlets_i) else []
        self.outgoing = pos[outlets_i[-1] :] if len(outlets_i) else []
        collat_mask = np.zeros(len(pos), dtype=bool)
        if len(outlets_i) and len(inlets_i):
            collat_mask[outlets_i[0] : inlets_i[-1] + 1] = True
        self.collat = pos[collat_mask & ~in_site]


class Path:
    """
    - pos:      set of coordinates containing all and only the points listed below
    - inlets:   the first point inside site each time obj enters site from scope
    - outlets:  the first point outside site (and in scope) each time obj leavs site
    - incoming: the sets of points from when object enters scope to the first inlet.
    - outgoing: the sets of points from the last outlet to when object leaves scope.
    - inside:   the sets of points when object is in site (and in scope).
    - collat:   the sets of point from the first outlet the last inlet (object leaves site but not scope).
    - subpaths: list of SubPath objects generated from this Path. See `trj_flow.SubPath`.

    types
    -----

    - pos:  n_subpaths lists of array (n_sub_points, 3)
    - inlets/outlets:   [] or array (n_subpaths, 3)
    - incoming/outgoing/inside/collat: [] or n_subpath lists of array (n_sub_points, 3)
    - subpaths: list of SubPath

    NOTE: all points that are not in scope are ignored, even if they are
    within site. Although this not enforced, site *should* always remain within scope.
    """

    def __init__(self, gids, pos, in_site, in_scope):
        """
        :gids: gid(s) of atom(s) of object molecule, array (n_atoms,)
        :pos: coordinates of obj center-of-mass during trajectory, array(n_times, 3) of folat
        :in_site: bitmask pos[in_site] will return all positions within site hull, array(n_times,) of bool
        :in_scope: bitmask pos[in_scope] will return all positions within scope hull, array(n_times,) of bool
        """
        self.gids = gids
        self._pos = pos
        self._in_site = in_site
        self._in_scope = in_scope
        self.subpaths = []
        self._make_subpaths()

    def _collect(self, attr):
        return [getattr(sub, attr) for sub in self.subpaths if len(getattr(sub, attr))]

    @property
    def inlets(self):
        items = self._collect("inlets")
        return np.vstack(items) if items else []

    @property
    def outlets(self):
        items = self._collect("outlets")
        return np.vstack(items) if items else []

    @property
    def pos(self):
        return self._collect("pos")

    @property
    def incoming(self):
        return self._collect("incoming")

    @property
    def outgoing(self):
        return self._collect("outgoing")

    @property
    def collat(self):
        return self._collect("collat")

    @property
    def inside(self):
        return self._collect("inside")

    def _make_subpaths(self):
        """
        scope = np.array([0,0,1,1,1,0,0,1,1,1,1,1,1,1,1,1,1,0,0])
        site  = np.array([0,0,0,0,0,0,0,0,1,1,0,0,1,1,1,0,0,0,0])
        chunks = [np.array(0,1,1,0,0,1,1,1,0,0)]
        """

        # break in_site into chunks corresponding to in_scope contiguous blocks
        pos_chunks = self._chunkize(self._pos, self._in_scope)
        in_site_chunks = self._chunkize(self._in_site, self._in_scope)
        for i in range(len(in_site_chunks)):
            if np.any(in_site_chunks[i]):
                self.subpaths.append(SubPath(pos_chunks[i], in_site_chunks[i]))

    @staticmethod
    def _chunkize(path, mask, all=False):
        """
        Break paths into time-contiguous chunks according to a time mask

        ```python
        mask = np.array([1,0,0,1,1,1,0,0,1,1,0,1]).astype(bool)
        np.split(mask.astype(int), np.nonzero(np.diff(mask))[0] + 1)[int(not mask[0])::2]
        [array([1]), array([1, 1, 1]), array([1, 1]), array([1])]
        ```

        :return: list of paths
        """
        gaps = np.nonzero(np.diff(mask))[0]
        s = slice(None, None) if all else slice(int(not mask[0]), None, 2)
        return np.split(path, gaps + 1)[s]

    @staticmethod
    def make_paths(trj, obj_gids_bymol, obj_in_site_matrix, obj_in_scope_matrix):
        """
        :obj_gids_bymol: obj gids grouped by molecule, array (n_mol, n_atoms_per_molecule) of int
        :obj_in_site_matrix: bitmask, array (n_times, n_obj_groups) when each obj is in site.
        :obj_in_scope_matrix: bitmask, array (n_times, n_obj_groups) when each obj is in scope.
        :return: list of Path objects
        """
        # mol_gids must be (1, n_atoms), but iterating on (n_obj, n_atoms) will yield (n_atoms,)
        # adding one dimension will do the trick (or use list-indexing)
        getpos = getposfunc(obj_gids_bymol.shape[-1])
        paths = []
        for mol_gids, in_site, in_scope in tqdm(zip( 
            obj_gids_bymol[:, None],
            obj_in_site_matrix.T,
            obj_in_scope_matrix.T
        ), total=len(obj_gids_bymol)):
            pos = np.array([getpos(fr, mol_gids).flatten() for fr in trj], dtype=np.float32)
            path = Path(mol_gids, pos, in_site, in_scope)
            paths.append(path)
        return paths


def getposfunc(size: int):
    """
    :return: center of mass of each object molecule, array (n_obj_molecules, 3)
    """

    def gidgroups2coms(fr: traj.Frame, groups):
        return np.array([fr.pos(group).mean(axis=0) for group in groups])

    def gidgroups2pos(fr: traj.Frame, groups):
        return fr.pos(groups.flatten())

    if size > 1:
        return gidgroups2coms
    else:
        return gidgroups2pos


def find_objects_in_hull(
    cms: topo.cms.Cms, trj: List[traj.Frame], obj_gids_bymol, hull_asl: str
):
    """
    Compute, for each frame, wether center of mass of `object` is within the
    convex hull of an atom `selection`.

    :ojb_gids_bymol: obj gids grouped by molecule, array (n_mol, n_atoms_per_molecule) of int
    :return 1: array (n_times, n_obj) of bool: M[i,j] = True if center of mass of
        obj `j` is within hull at time `i`.
    :return 2: list of hull objects 
    """

    getpos = getposfunc(obj_gids_bymol.shape[-1])

    hull_gids = topo.asl2gids(cms, hull_asl)  # variable (na_stite!,)
    is_dynamic = topo.is_dynamic_asl(cms, hull_asl)

    hulls = []
    is_in_hull_matrix = np.empty((len(trj), len(obj_gids_bymol)), dtype=bool)
    for i, fr in enumerate(tqdm(trj)):
        p_obj = getpos(fr, obj_gids_bymol)  # (na_obj, 3)

        # dynamic asl
        if is_dynamic:
            cms.setXYZ(fr.pos(cms.allaid_gids))
            hull_gids = topo.asl2gids(cms, hull_asl)  # variable (na_stite!,)

        points = fr.pos(hull_gids)  # ndarray(na_site!, 3)
        hulls.append(ConvexHull(points))

        # find which object atoms are found within site hull
        is_in_hull_matrix[i] = (
            Delaunay(points).find_simplex(p_obj) >= 0
        )  # bool ndarray(na_obj,)

    return is_in_hull_matrix, hulls

def get_centroids(points, eps=0.5, min_samples=5):
    db =  DBSCAN(eps=eps, min_samples=min_samples, metric='euclidean').fit(points)
    labels = db.labels_
    cc = np.array([points[labels == v].mean(axis=0) for v in set(labels) if v != -1])
    return cc

def make_hulls(cms, trj, hull_asl):
    hull_gids = topo.asl2gids(cms, hull_asl)
    is_dynamic = topo.is_dynamic_asl(cms, hull_asl)

    # allpos = []
    hulls = []
    for fr in tqdm(trj):
        if is_dynamic:
            cms.setXYZ(fr.pos(cms.allaid_gids))
            hull_gids = topo.asl2gids(cms, hull_asl)

        hull_pos = fr.pos(hull_gids)
        hull = ConvexHull(hull_pos)
        # allpos.append(hull_pos)
        hulls.append(hull)
    return hulls # ConvexHull(np.vstack(allpos))

def main():
    parser = argparse.ArgumentParser(description="wip")
    parser.add_argument("cms", help="input topology file")
    parser.add_argument("out", help="base name for the output")
    parser.add_argument(
        "-t",
        help="trj dir",
        metavar="trajectory",
    )
    parser.add_argument(
        "-j",
        help="number of processes (max physical cores), default 1",
        type=int,
        metavar="N cores",
        default=1,
    )
    parser.add_argument(
        "-site_asl",
        help="atom selection, default Ca within 6 Å from ligand",
        default="a.pt CA and (fillres (within 6 (ligand)))",
    )
    parser.add_argument(
        "-scope_asl", help="atom selection, default Ca", default="a.pt CA"
    )
    parser.add_argument(
        "-obj_asl",
        help="atom selection, default water Oxygen",
        default="water and a.e O",
    )
    parser.add_argument("-s", help="slicer", metavar="START:END:STEP")
    parser.add_argument('-smooth', nargs='?', default=False, const='10,0,2', help='''Smooth obj positions using scipy wrappers of FITPACK.
            It can be optionally followed by a comma-spearated list of values N,S,O
            indicating the new number of points factor for interpolation (len(xi) == N*len(x)),
            the smooth factor and spline order. Default 10,0,2''')
    parser.add_argument('-save', nargs='?', const=None, default=False, help='store paths data to FILE', metavar='FILE')
    parser.add_argument('-load', help='load paths data from FILE. Useful for changing visualization parameters', metavar='FILE')
    parser.add_argument('-clust', default='1,5', help='comma-separated list of DBSCAN clustering parameters: eps,min_samples. Default 1,5')
    parser.add_argument('-ppa', default=2, type=int, help='points-per-Ångstrom: determine the volume resolution for displaying.')
    args = parser.parse_args()
    # args = parser.parse_args(["aligned-out.cms", "flowtest", "-obj_asl", "water"])

    scope_asl = args.scope_asl
    obj_asl = args.obj_asl
    site_asl = args.site_asl

    if args.t:
        msys, cms = topo.read_cms(args.cms)
        trj = traj.read_traj(args.t)
    else:
        msys, cms, trj = traj_util.read_cms_and_traj(args.cms)

    slicer = (
        slice(*[int(i) if i else None for i in args.s.split(":")])
        if args.s
        else slice(None, None)
    )

    trj = trj[slicer]

    if args.load:
        print(f'Intermediate data loaded from {args.load}')
        with open(args.load, 'rb') as f:
            paths, site_hulls, scope_hulls = pickle.load(f)
    else:
        obj_aids_bymol = analyze.group_by_connectivity(cms, cms.select_atom(obj_asl))  # (n_mol, n_at_per_mol)
        obj_gids_bymol = np.array([topo.aids2gids(cms, group) for group in obj_aids_bymol], dtype=int, ndmin=2)

        print("Selecting objects to track...", flush=True)
        obj_in_site_matrix, site_hulls = find_objects_in_hull(cms, trj, obj_gids_bymol, site_asl)

        # find indices of the objects that enter site at least once
        obj_ok_i = np.nonzero(np.any(obj_in_site_matrix, axis=0))[0]
        obj_gids_bymol_ok = obj_gids_bymol[obj_ok_i, :]
        print("Done.")

        print("Checking scope...", flush=True)
        obj_in_scope_matrix, scope_hulls = find_objects_in_hull(cms, trj, obj_gids_bymol_ok, scope_asl)
        print("Done.")

        print("Building Paths...", flush=True)
        paths = Path.make_paths(
            trj,
            obj_gids_bymol_ok,
            obj_in_site_matrix[:, obj_ok_i],
            obj_in_scope_matrix,
        )

        if args.save != False:
            fname = args.save or args.out
            print(f'Intermediate data stored to {fname}')
            with open(fname, 'wb') as f:
                pickle.dump((paths, site_hulls, scope_hulls), f)

    if args.smooth:
        _interpf, _smoothf, _splorder = args.smooth.split(',')
        _interpf, _smoothf, _splorder = int(_interpf), float(_smoothf), int(_splorder)
    else:
        _interpf, _smoothf, _splorder = [None] * 3
    paths = PathCollection(paths, interpf=_interpf, smoothf=_smoothf, splorder=_splorder)
    print("Done.")

    print("Clustering...", end=' ', flush=True)
    _eps, _min_samples = args.clust.split(',')
    _eps, _min_samples = float(_eps), int(_min_samples)
    inlets_c = get_centroids(paths.inlets, eps=_eps, min_samples=_min_samples) if len(paths.inlets) else []
    outlets_c = get_centroids(paths.outlets, eps=_eps, min_samples=_min_samples) if len(paths.outlets) else []
    print("Done.")

    # print("Calculating convex hulls for displaying...", flush=True)
    # scope_hulls = make_hulls(cms, trj, scope_asl)
    # site_hulls = make_hulls(cms, trj, site_asl)
    # print("Done.")
    # scope_hull = ConvexHull(trj[0].pos(topo.asl2gids(cms, scope_asl)))
    # site_hull = ConvexHull(trj[0].pos(topo.asl2gids(cms, site_asl)))
    # scope_p = trj[0].pos(topo.asl2gids(cms, scope_asl))

    print("Rendering...", flush=True)
    script = render_pymol(paths, site_hulls, scope_hulls, inlets_c, outlets_c, ppa=args.ppa, name=args.out)
    with open(args.out + "_pymol.py", "w") as f:
        f.write(script)
        print("  " + args.out + "_pymol.py", 'written.')
    print('Done.')

    print("All done.\n\n")
    print("""
      .-``'.     S T A Y   W I T H     .'''-.
    .`   .`~           T H E            ~`.   '.
_.-'     '._        ~ F L O W ~         _.'     '-._
  ,(   ,(   ,(   ,(   ,(   ,(   ,(   ,(   ,(   ,(   ,(
`-'  `-'  `-'  `-'  `-'  `-'  `-'  `-'  `-'  `-'  `-' `
""")
    return


def render_pymol(paths, site_hulls, scope_hulls, inlets_c, outlets_c, ppa=10, name='flow'):

    def plot_hull(hull, *args, **kwargs):
        hull_cgo = Cgo()
        for simplex in hull.simplices:  # type: ignore
            facet = np.hstack((simplex, simplex[0]))
            hull_cgo.scatter(hull.points[facet], *args, **kwargs)
        return hull_cgo

    def plot_paths(paths, attr, *args, **kwargs):
        path_cgo = Cgo()
        # for path in paths:
        for p in getattr(paths, attr):
            # if len(p):
            path_cgo.scatter(p, *args, **kwargs)
        return path_cgo


    site_hull_cgos = [plot_hull(hull, linecolor=to_rgb("cyan"), marker=True, radius=0.08) for hull in site_hulls]
    scope_hull_cgos = [plot_hull(hull, linecolor=to_rgb("gray"), marker=True, radius=0.08) for hull in scope_hulls]
    inlets_cgo = Cgo.Scatter(paths.inlets, linecolor=to_rgb("mediumslateblue"), line=False, marker=True, radius=.1,)
    outlets_cgo = Cgo.Scatter(paths.outlets, linecolor=to_rgb("tomato"), line=False, marker=True, radius=.1,)
    inlets_c_cgo = Cgo.Scatter(inlets_c, linecolor=to_rgb("mediumslateblue"), line=False, marker=True, radius=.5,)
    outlets_c_cgo = Cgo.Scatter(outlets_c, linecolor=to_rgb("tomato"), line=False, marker=True, radius=.5,)
    incoming_cgo = plot_paths(paths, 'incoming', linecolor=to_rgb("mediumslateblue"))
    outgoing_cgo = plot_paths(paths, 'outgoing', linecolor=to_rgb("tomato"))
    inside_cgo = plot_paths(paths, 'inside', linecolor=to_rgb("green"))
    collat_cgo = plot_paths(paths, 'collat', linecolor=to_rgb("gold"))

    allpos = np.vstack(paths.pos)
    extent = np.ptp(allpos, axis=0)
    nbins = np.round(extent * ppa).astype(int)
    # density, edges = np.histogramdd(allpos, bins=nbins)
    density, edges = np.histogramdd(allpos, bins=nbins, density=True)
    density = -np.log(np.maximum(density, 1e-16))
    spacing = np.ptp(allpos, axis=0) / nbins
    minmax = np.array([(e.min(), e.max()) for e in edges]).T

    data = {
        "name": name,
        "site_hull": [cgo.obj for cgo in site_hull_cgos],
        "scope_hull": [cgo.obj for cgo in scope_hull_cgos],
        "inlets": inlets_cgo.obj,
        "outlets": outlets_cgo.obj,
        "inlets_c": inlets_c_cgo.obj,
        "outlets_c": outlets_c_cgo.obj,
        "incoming": incoming_cgo.obj,
        "outgoing": outgoing_cgo.obj,
        "inside": inside_cgo.obj,
        "collat": collat_cgo.obj,
        "density": density,
        "gridspacing": spacing,
        "minmax": minmax
        }

    with open(name + '_pymoldata.pickle', 'wb') as f:
        pickle.dump(data, f)
        print("  " + name + "_pymoldata.pickle", 'written.')

    script = f"""
# Generated by flow.py

import pickle
import numpy as np
from pymol import cmd
from chempy.brick import Brick

with open(f'{name}_pymoldata.pickle', 'rb') as f:
    data = pickle.load(f)

name = data['name']

cmd.load_cgo(data["incoming"], name + "_incoming", 1)
cmd.load_cgo(data["outgoing"], name + "_outgoing", 1)
cmd.load_cgo(data["inside"], name + "_inside", 1)
cmd.load_cgo(data["collat"], name + "_collat", 1)
cmd.load_cgo(data["inlets"], name + "_inlets", 1)
cmd.load_cgo(data["outlets"], name + "_outlets", 1)
cmd.load_cgo(data["inlets_c"], name + "_inlets_centroids", 1)
cmd.load_cgo(data["outlets_c"], name + "_outlets_centroids", 1)

for i, (site_hull, scope_hull) in enumerate(zip(data['site_hull'], data['scope_hull'])):
    cmd.load_cgo(site_hull, name + "_site_hull", i + 1)
    cmd.load_cgo(scope_hull, name + "_scope_hull", i + 1)

density = data['density']
gridspacing = data['gridspacing']
mn, mx = data['minmax']
brick = Brick.from_numpy(density, gridspacing)
brick.origin = list(mn)
brick.range = list(mx - mn)
cmd.load_brick(brick, name + '_map')

alpha = 0.2
colors = reversed([
    [0.00, 0.00, 1.00],
    [0.00, 1.00, 1.00],
    [0.00, 1.00, 0.00],
    [1.00, 1.00, 0.00],
    [1.00, 0.50, 0.00],
    [1.00, 0.00, 0.00],
])
# dnz = np.log(density[np.nonzero(density)])
dnz = density[density > 1e-16]
spacing = np.ptp(dnz)/6
w = spacing/4
ramp_iso = np.linspace(dnz.min() + w, dnz.max() - w, 6)

ramp = []
for iso, c in zip(ramp_iso, colors):
    ramp += [iso - w] + c + [0]
    ramp += [iso]     + c + [alpha]
    ramp += [iso + w] + c + [0]

cmd.volume_ramp_new('flow_ramp', ramp)
cmd.volume(name + "_volume", name + "_map", "flow_ramp")
"""
    return script


if __name__ == "__main__":
    def testplot(path, win, poly):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        fig = plt.figure(figsize=plt.figaspect(0.5))
        
        axx = fig.add_subplot(3,2,1)
        axy = fig.add_subplot(3,2,3)
        axz = fig.add_subplot(3,2,5)
        ax3 = fig.add_subplot(1,2,2, projection='3d')

        axx.plot(path[:,0])
        axy.plot(path[:,1])
        axz.plot(path[:,2])
        ax3.plot3D(*path.T)

        path_s = savgol_filter(path, axis=0, window_length=win, polyorder=poly)

        axx.plot(path_s[:,0])
        axy.plot(path_s[:,1])
        axz.plot(path_s[:,2])
        ax3.plot3D(*path_s.T)
        ax3.set_title(f'{win=} {poly=}')
        plt.tight_layout()

    def testplot(path, s=2, N=1000):
        # https://stackoverflow.com/questions/18962175/spline-interpolation-coefficients-of-a-line-curve-in-3d-space
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from scipy import interpolate
        fig = plt.figure(figsize=plt.figaspect(0.5))
        
        axx = fig.add_subplot(3,2,1)
        axy = fig.add_subplot(3,2,3)
        axz = fig.add_subplot(3,2,5)
        ax3 = fig.add_subplot(1,2,2, projection='3d')

        x, y, z = path.T
        axx.plot(np.linspace(x.min(), x.max(), len(x)), x)
        axy.plot(np.linspace(y.min(), y.max(), len(y)), y)
        axz.plot(np.linspace(z.min(), z.max(), len(z)), z)
        ax3.plot3D(*path.T, '-o')

        tck, u = interpolate.splprep(path.T.tolist(), s=s)
        x, y, z = interpolate.splev(np.linspace(u.min(),u.max(),N), tck)


        axx.plot(np.linspace(x.min(), x.max(), len(x)), x)
        axy.plot(np.linspace(y.min(), y.max(), len(y)), y)
        axz.plot(np.linspace(z.min(), z.max(), len(z)), z)
        ax3.plot3D(x, y, z)
        ax3.set_title(f'{s=} {N=}')
        plt.tight_layout()

    def testclust(points, eps=0.5, min_samples=5):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D
        from sklearn.metrics import silhouette_score

        db = DBSCAN(eps=eps, min_samples=min_samples, metric='euclidean').fit(points)
        labels = db.labels_
        core_samples_mask = np.zeros_like(labels, dtype=bool)
        core_samples_mask[db.core_sample_indices_] = True
        n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
        n_noise_ = list(labels).count(-1)
        print("Estimated number of clusters: %d" % n_clusters_)
        print("Estimated number of noise points: %d" % n_noise_)
        print("Silhouette Coefficient: %0.3f" % silhouette_score(points, labels))

        ax = Axes3D(plt.figure())
        ax.scatter3D(*points[core_samples_mask].T, c=labels[core_samples_mask], s=50)
        ax.scatter3D(*points[~core_samples_mask].T, c=labels[~core_samples_mask], s=1)



    
    def test():
        """
        doctest
            >>> test()
            pos
            [4 5 6 7] [ 9 10 11 12 13 14 15 16] [21 22] [24] [26 27]
            inlets
            [5] [10] [14] [27]
            outlets
            [7] [12] [16] [22]
            incoming
            [4 5] [ 9 10] [26 27]
            outgoing
            [7] [16] [22]
            inside
            [5 6] [10 11 14 15] [21] [24] [27]
            collat
            [12 13 14]
        """
        scope = np.array([0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 1, 1, 0], dtype=bool,)
        site =  np.array([0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1], dtype=bool,)
        pos =   np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28], dtype=int,)
        p = Path(0, pos[:, None], site, scope)
        pc = PathCollection([p])
        pp = lambda pos: print(" ".join("".join(str(l.flatten())) for l in pos))
        print("pos")
        pp(pc.pos)
        print("inlets")
        pp(pc.inlets)
        print("outlets")
        pp(pc.outlets)
        print("incoming")
        pp(pc.incoming)
        print("outgoing")
        pp(pc.outgoing)
        print("inside")
        pp(pc.inside)
        print("collat")
        pp(pc.collat)
        print()
        for i, sp in enumerate(p.subpaths):
            print("---", "subpath", i, "---")
            print("pos")
            pp(sp.pos)
            print("inside")
            pp(sp.inside)
            print("inlets")
            pp(sp.inlets)
            print("outlets")
            pp(sp.outlets)
            print("incoming")
            pp(sp.incoming)
            print("outgoing")
            pp(sp.outgoing)
            print("collat")
            pp(sp.collat)
            print()

    # test()
    main()
