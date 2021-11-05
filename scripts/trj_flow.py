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
from dataclasses import dataclass
from typing import List

import numpy as np
import plotly.graph_objects as go
from matplotlib.colors import to_rgb
from schrodinger.application.desmond.packages import topo, traj, traj_util
from schrodinger.structutils import analyze
from scipy.spatial import ConvexHull, Delaunay
from tqdm import tqdm

# https://stackoverflow.com/questions/16750618/whats-an-efficient-way-to-find-if-a-point-lies-in-the-convex-hull-of-a-point-cl


# def runner(fr, obj_gids, hull, hull_asl=None, cms=None, hull_gids=None):
#     p_obj = fr.pos(obj_gids) # (na_obj, 3)
#
#     # dynamic asl
#     if hull_asl:
#         cms.setXYZ(fr.pos(cms.allaid_gids))
#         hull_gids = topo.asl2gids(cms, hull_asl) # variable (na_stite!,)
#
#     points = fr.pos(hull_gids) # ndarray(na_site!, 3)
#     hull.add_points(points)
#
#     # find which object atoms are found within site hull
#     return Delaunay(points).find_simplex(p_obj) >= 0 # bool ndarray(na_obj,)


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
    :return: array (n_times, n_obj) of bool: M[i,j] = True if center of mass of
    obj `j` is within hull at time `i`.
    """

    getpos = getposfunc(obj_gids_bymol.shape[-1])

    # initialize site hull
    hull_gids = topo.asl2gids(cms, hull_asl)  # variable (na_stite!,)

    is_in_hull_matrix = np.empty((len(trj), len(obj_gids_bymol)), dtype=bool)
    for i, fr in enumerate(tqdm(trj)):
        p_obj = getpos(fr, obj_gids_bymol)  # (na_obj, 3)

        # dynamic asl
        if topo.is_dynamic_asl(cms, hull_asl):
            cms.setXYZ(fr.pos(cms.allaid_gids))
            hull_gids = topo.asl2gids(cms, hull_asl)  # variable (na_stite!,)

        points = fr.pos(hull_gids)  # ndarray(na_site!, 3)

        # find which object atoms are found within site hull
        is_in_hull_matrix[i] = (
            Delaunay(points).find_simplex(p_obj) >= 0
        )  # bool ndarray(na_obj,)

    return is_in_hull_matrix


class SubPath:
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
        self.collat = (
            pos[outlets_i[0] : inlets_i[-1] + 1]
            if len(outlets_i) and len(inlets_i)
            else []
        )


class Path:
    """
    - pos:      set of coordinates containing all and only the points listed below
    - inlets:   the first point inside site each time obj enters site from scope
    - outlets:  the first point outside site (and in scope) each time obj leavs site
    - incoming: the points from when object enters scope to the first inlet.
    - outgoing: the points from the last outlet to when object leaves scope.
    - inside:   the points when object is in site (and in scope).
    - collat:   the point from the first outlet the last inlet (object leaves site but not scope).

    types
    -----

    - pos:  n_subpaths lists of array (n_sub_points, 3)
    - inlets/outlets:   [] or array (n_subpaths, 3)
    - incoming/outgoing/inside/collat: [] or n_subpath lists of array (n_sub_points, 3)

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
        getpos = getposfunc(obj_gids_bymol.shape[-1])
        paths = []
        for mol_gids, in_site, in_scope in zip(
            obj_gids_bymol, obj_in_site_matrix.T, obj_in_scope_matrix.T
        ):
            pos = np.array([getpos(fr, [ mol_gids ]).flatten() for fr in trj])
            path = Path(mol_gids, pos, in_site, in_scope)
            paths.append(path)
        return paths


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

    # TODO: maybe add support for variable n_atom_per_mol (multiple obj types)
    obj_aids_bymol = analyze.group_by_connectivity(
        cms, cms.select_atom(obj_asl)
    )  # (n_mol, n_at_per_mol)
    obj_gids_bymol = np.array(
        [topo.aids2gids(cms, group) for group in obj_aids_bymol], dtype=int, ndmin=2
    )

    print("Selecting objects to track...")
    obj_in_site_matrix = find_objects_in_hull(cms, trj, obj_gids_bymol, site_asl)
    # obj_in_site_matrix = pickle.load(open('in_site_test', 'rb'))

    # find indices of the objects that enter site at least once
    obj_ok_i = np.nonzero(np.any(obj_in_site_matrix, axis=0))[0]
    obj_gids_bymol_ok = obj_gids_bymol[obj_ok_i, :]
    print("Done.")

    print("Checking scope...")
    obj_in_scope_matrix = find_objects_in_hull(cms, trj, obj_gids_bymol_ok, scope_asl)
    # obj_in_scope_matrix = pickle.load(open('in_scope_test', 'rb'))
    print("Done.")

    print("Building Paths...")
    paths = Path.make_paths(
        trj,
        obj_gids_bymol_ok,
        obj_in_site_matrix[:, obj_ok_i],
        obj_in_scope_matrix,
    )
    print("Done.")

    scope_hull = ConvexHull(trj[0].pos(topo.asl2gids(cms, scope_asl)))
    site_hull = ConvexHull(trj[0].pos(topo.asl2gids(cms, site_asl)))
    scope_p = trj[0].pos(topo.asl2gids(cms, scope_asl))

    print("Rendering...")
    # fig = plot(paths, scope_p, scope_hull, site_hull)
    # fig.show()

    script = plotcgo(paths, site_hull, scope_hull)
    with open(args.out + "_pymol.py", "w") as f:
        f.write(script)

    return


def plotcgo(paths, site_hull, scope_hull):
    @dataclass
    class Cgo:
        BEGIN: float
        END: float
        LINE_LOOP: float
        SPHERE: float
        VERTEX: float
        LINEWIDTH: float
        COLOR: float

    cgo = Cgo(
        **{
            "COLOR": 6.0,
            "BEGIN": 2.0,
            "END": 3.0,
            "LINE_LOOP": 2.0,
            "SPHERE": 7.0,
            "VERTEX": 4.0,
            "LINEWIDTH": 10.0,
        }
    )

    def scatter(
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
            obj.extend(
                [
                    cgo.BEGIN,
                    cgo.LINE_LOOP,
                    cgo.COLOR,
                    *linecolor,
                    cgo.LINEWIDTH,
                    linewidth,
                ]
            )
            for xyz in points:
                obj.extend([cgo.VERTEX, *xyz])
            obj.extend([cgo.END])

        if marker:
            obj.extend([cgo.COLOR, *(markercolor or linecolor)])
            for xyz in points:
                obj.extend([cgo.SPHERE, *xyz, radius])

        return [float(i) for i in obj]

    site_hull_obj = []
    for simplex in site_hull.simplices:  # type: ignore
        facet = np.hstack((simplex, simplex[0]))
        site_hull_obj += scatter(
            site_hull.points[facet], linecolor=to_rgb("cyan"), marker=True, radius=0.05
        )

    scope_hull_obj = []
    for simplex in scope_hull.simplices:  # type: ignore
        facet = np.hstack((simplex, simplex[0]))
        scope_hull_obj += scatter(
            scope_hull.points[facet], linecolor=to_rgb("gray"), marker=True, radius=0.05
        )

    inlets_obj = scatter(
        np.vstack([p.inlets for p in paths if len(p.inlets)]),
        linecolor=to_rgb("lightblue"),
        line=False,
        marker=True,
        radius=1,
    )
    outlets_obj = scatter(
        np.vstack([p.outlets for p in paths if len(p.outlets)]),
        linecolor=to_rgb("purple"),
        line=False,
        marker=True,
        radius=1,
    )

    incoming_obj, outgoing_obj = [], []
    inside_obj, collat_obj = [], []
    for path in paths:
        for p in path.incoming:
            if len(p):
                incoming_obj += scatter(p, linecolor=to_rgb("lightblue"))
        for p in path.outgoing:
            if len(p):
                outgoing_obj += scatter(p, linecolor=to_rgb("purple"))
        for p in path.inside:
            if len(p):
                inside_obj += scatter(p, linecolor=to_rgb("green"))
        for p in path.collat:
            if len(p):
                collat_obj += scatter(p, linecolor=to_rgb("orange"))

    script = f"""
from pymol import cmd

cmd.load_cgo({site_hull_obj}, "site_hull", 1)
cmd.load_cgo({scope_hull_obj}, "scope_hull", 1)
cmd.load_cgo({incoming_obj}, "incoming", 1)
cmd.load_cgo({outgoing_obj}, "outgoing", 1)
cmd.load_cgo({inside_obj}, "inside", 1)
cmd.load_cgo({collat_obj}, "collat", 1)
cmd.load_cgo({inlets_obj}, "inlets", 1)
cmd.load_cgo({outlets_obj}, "outlets", 1)
    """
    return script


def plot(paths, p_protein, scope_hull, site_hull):
    density, (x, y, z) = np.histogramdd(
        np.vstack([chunk for path in paths for chunk in path.pos]), bins=105, density=1
    )
    X, Y, Z = np.meshgrid(x[1:], y[1:], z[1:])
    vol = go.Volume(
        x=X.flatten(),
        y=Y.flatten(),
        z=Z.flatten(),
        value=density.flatten(),
        isomin=0,
        isomax=density.max(),
        opacity=0.1,  # needs to be small to see through all surfaces
        surface_count=17,  # needs to be a large number for good volume rendering
    )

    scope = go.Scatter3d(
        **{k: v for (k, v) in zip("xyz", p_protein.T)},
        line=dict(color="darkblue", width=5),
        marker=dict(size=1),
    )

    scope_facets = []
    for simplex in scope_hull.simplices:  # type: ignore
        facet = np.hstack((simplex, simplex[0]))
        scope_facets.append(
            go.Scatter3d(
                **{k: v for (k, v) in zip("xyz", scope_hull.points[facet].T)},
                line=dict(color="black", width=0.5),
                marker=dict(size=0.5),
            )
        )
    site_facets = []
    for simplex in site_hull.simplices:  # type: ignore
        facet = np.hstack((simplex, simplex[0]))
        site_facets.append(
            go.Scatter3d(
                **{k: v for (k, v) in zip("xyz", site_hull.points[facet].T)},
                line=dict(color="blue", width=5),
                marker=dict(size=2),
            )
        )
    inlets = go.Scatter3d(
        **{
            k: v
            for (k, v) in zip(
                "xyz", np.vstack([p.inlets for p in paths if len(p.inlets)]).T
            )
        },
        marker=dict(color="green", size=2),
        mode="markers",
    )
    outlets = go.Scatter3d(
        **{
            k: v
            for (k, v) in zip(
                "xyz", np.vstack([p.outlets for p in paths if len(p.outlets)]).T
            )
        },
        marker=dict(color="red", size=2),
        mode="markers",
    )

    lines = []
    for path in paths:
        for p in path.incoming:
            if len(p):
                lines.append(
                    go.Scatter3d(
                        **{k: v for (k, v) in zip("xyz", p.T)},
                        line=dict(color="blue", width=0.5),
                        mode="lines",
                    )
                )
        for p in path.outgoing:
            if len(p):
                lines.append(
                    go.Scatter3d(
                        **{k: v for (k, v) in zip("xyz", p.T)},
                        line=dict(color="red", width=0.5),
                        mode="lines",
                    )
                )
        for p in path.collat:
            if len(p):
                lines.append(
                    go.Scatter3d(
                        **{k: v for (k, v) in zip("xyz", p.T)},
                        line=dict(color="orange", width=0.5),
                        mode="lines",
                    )
                )
        for p in path.inside:
            if len(p):
                lines.append(
                    go.Scatter3d(
                        **{k: v for (k, v) in zip("xyz", p.T)},
                        line=dict(color="green", width=0.5),
                        mode="lines",
                    )
                )
    fig = go.Figure(data=(scope, *scope_facets, *site_facets, inlets, outlets, *lines))
    fig.update_layout(showlegend=False)
    return fig


if __name__ == "__main__":

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
        scope = np.array(
            [
                0,
                1,
                1,
                0,
                1,
                1,
                1,
                1,
                0,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                0,
                0,
                0,
                0,
                1,
                1,
                0,
                1,
                0,
                1,
                1,
                0,
            ],
            dtype=bool,
        )
        site = np.array(
            [
                0,
                0,
                0,
                0,
                0,
                1,
                1,
                0,
                0,
                0,
                1,
                1,
                0,
                0,
                1,
                1,
                0,
                0,
                1,
                1,
                0,
                1,
                0,
                0,
                1,
                0,
                0,
                1,
                1,
            ],
            dtype=bool,
        )
        pos = np.array(
            [
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9,
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                28,
            ],
            dtype=int,
        )
        p = Path(0, pos[:, None], site, scope)
        pp = lambda pos: print(" ".join("".join(str(l.flatten())) for l in pos))
        print("pos")
        pp(p.pos)
        print("inlets")
        pp(p.inlets)
        print("outlets")
        pp(p.outlets)
        print("incoming")
        pp(p.incoming)
        print("outgoing")
        pp(p.outgoing)
        print("inside")
        pp(p.inside)
        print("collat")
        pp(p.collat)
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
