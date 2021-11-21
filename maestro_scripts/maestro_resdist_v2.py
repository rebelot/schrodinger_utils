import numpy as np
from schrodinger import maestro
from schrodinger.structutils import analyze
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize

def read_align(fasta):
    fasta = open(fasta, 'r')
    sequences = []
    for line in fasta.readlines():
        if line[0] == '>':
            seq = []
            sequences.append(seq)
        else:
            seq += list(line.strip())
    return np.array(sequences)

def seq2atom(seq, Ca):
    atoms = []
    rindex = 0
    for res in seq:
        if res != '-':
            atoms.append(Ca[rindex])
            rindex +=1
        else:
            atoms.append(None)
    return atoms

def color(sm, dists, atoms, row):
    row.includeOnly()
    for a, dist in zip(atoms, dists):
        if a:
            rgba = sm.to_rgba(dist) if dist else [1, 0, 0, 0]
            rgb = [int(255*c) for c in rgba[:3]]
            command = f'coloratomrgb red={rgb[0]} green={rgb[1]} blue={rgb[2]} protein and r.n {a.resnum}'
            maestro.command(command)
            maestro.redraw_request()

def calcavgdist(atoms):
    dists = np.empty_like(atoms, dtype=float)

    for i, block in enumerate(atoms.T):
        if None in block:
            dists[:,i] = None
        else:
            xyz = np.array([a.xyz for a in block])
            avg = xyz.mean(axis=0)
            dists[:, i] = np.linalg.norm(avg - xyz, axis=1)
    return dists

def prd(fasta):
    """
    Color residues of multiple entries by their distances to an average position.
    Residues across entries are paired according to an alignment file.

    Entries must be selected in the project table and their order
    be the same in which they appear in the alignment file.

    Works only with single-chain strucures.
    """
    sequences = read_align(fasta)

    pt = maestro.project_table_get()
    rows = list(pt.selected_rows)
    Ca = [list(analyze.get_atoms_from_asl(row.getStructure(), 'a.pt CA')) for row in rows]
    atoms = np.array([seq2atom(seq, ca) for (seq, ca) in zip(sequences, Ca)])
    dists = calcavgdist(atoms)

    vmin = dists[~np.isnan(dists)].min()
    vmax = dists[~np.isnan(dists)].max()
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm, 'RdYlGn_r')
    for d, a, r in zip(dists, atoms, rows):
        color(sm, d, a, r)


