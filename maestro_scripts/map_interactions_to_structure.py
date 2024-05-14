import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from schrodinger.graphics3d import common, lines


class Interactions:
    def __init__(self, labels, matrix):
        self.labels: list[str] = labels
        self.matrix = matrix

    @classmethod
    def read(cls, *files):
        def readinter(file):
            labels = np.genfromtxt(
                file, usecols=(0), delimiter=",", skip_header=1, dtype=str
            )
            matrix = np.genfromtxt(file, delimiter=",", skip_header=1)
            matrix = np.array(matrix, ndmin=2)[:, 1:].astype(int)
            return labels, matrix

        multi_labels = []
        multi_matrix = []
        for file in files:
            labels, matrix = readinter(file)
            multi_labels.extend(labels)
            multi_matrix.extend(matrix)
        return cls(multi_labels, np.array(multi_matrix, ndmin=2))

    def split_labels_and_merge(self):
        interdict = {}
        for label, row in zip(self.labels, self.matrix.astype(bool)):
            r1, r2 = label.split(" - ")
            interdict.setdefault(r1, np.zeros(len(row)).astype(bool))
            interdict.setdefault(r2, np.zeros(len(row)).astype(bool))
            interdict[r1] |= row
            interdict[r2] |= row

        for k, val in interdict.items():
            interdict[k] = val.mean()

        return interdict

    def get_inter_dict(self):
        interdict = {}
        for label, row in zip(self.labels, self.matrix.astype(bool)):
            interdict.setdefault(label, np.zeros(len(row)).astype(bool))
            interdict[label] |= row

        for k, val in interdict.items():
            interdict[k] = val.mean()
        return interdict


def key2asl(key):
    chain, resname, resnum = key.split()
    return f'r.n {resnum} and c.n {chain} and r.pt "{resname:>3s}"'


def color_atoms(data: dict, cmap):
    from schrodinger import maestro

    ws = maestro.workspace_get()

    vmin = min(data.values())
    vmax = max(data.values())
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm, cmap)

    for key, val in data.items():
        res_atoms = list(maestro.analyze.get_atoms_from_asl(ws, key2asl(key)))
        if not res_atoms:
            raise ValueError(f"no atoms selected for {key} -> {key2asl(key)}")
        color = [int(255 * c) for c in sm.to_rgba(val)[0:3]]
        for a in res_atoms:
            a.color = color

    maestro.workspace_set(ws)
    return ws


def parse_residue_interactions_and_color_residues(
    *filenames, cmap="Greens", thresh=0.3, color="red"
):
    inter = Interactions.read(*filenames)
    data = inter.split_labels_and_merge()
    interdict = inter.get_inter_dict()
    color_atoms(data, cmap)
    linegrp = make_interaction_lines(interdict, thresh=thresh, color=color)
    return inter, linegrp


def atoms2pos(atoms):
    atoms = list(atoms)
    # try:
    #     return atoms[0].getResidue().getAlphaCarbon().xyz
    # except AttributeError:
    #     return np.array([a.xyz for a in atoms]).mean(axis=0)
    return np.array([a.xyz for a in atoms]).mean(axis=0)


def write_interdict(filename, interdict):
    with open(filename, "w") as f:
        for k, v in interdict.items():
            f.write(k + "," + str(v) + "\n")


def make_interaction_lines(interdict, thresh=0.3, color="red"):
    from schrodinger import maestro

    ws = maestro.workspace_get()
    linegrp = common.Group()
    wmin = min(interdict.values())
    wmax = max(interdict.values())
    # segments = []
    for label, val in interdict.items():
        key1, key2 = label.split(" - ")
        res_1 = maestro.analyze.get_atoms_from_asl(ws, key2asl(key1))
        res_2 = maestro.analyze.get_atoms_from_asl(ws, key2asl(key2))
        segment = [atoms2pos(res_1), atoms2pos(res_2)]
        # segments.append(segment)
        if val >= thresh:
            w = 10 * (val - wmin) / (wmax - wmin) + 0.01
            linegrp.add(lines.MaestroLines([segment], color=color, width=w))
    return linegrp
