__doc__ = """

Map data to atom colors

"""
# Name: Data to atom color
# Command: pythonrun maestro_applycscheme map2color

import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from schrodinger.structutils import analyze


def from_file(datafile):
    data = []
    chainresnum = []
    with open(datafile, "r") as f:
        lines = f.readlines()
    for line in lines:
        chain, resnum, val = line.split()
        data.append(float(val))
        chainresnum.append((chain, resnum))

    return data, chainresnum


def map2color(st, data, chainresnum, cmap):
    """
    :st: maestro.Structure
    :data: 1D np.array of length st.residue
    :cmap: string (available matplotlib colormaps)
    """

    data = np.array(data)
    vmin = data.min()
    vmax = data.max()
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm, cmap)
    for res, val in zip(chainresnum, data):
        c_alpha = list(
            analyze.get_atoms_from_asl(st, f"c.n {res[0]} and r.n {res[1]} and a.pt CA")
        )
        c_alpha = c_alpha[0] if c_alpha else None
        if c_alpha is not None:
            c_alpha.color = [int(255 * c) for c in sm.to_rgba(val)[:-1]]
    return st
