import re

import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from schrodinger import maestro

CT = "m_psp_residue_interaction_energies"


class InteractionEnergies:
    inter_ene_keys = [
        "s_psp_Prime_Coulomb",
        "s_psp_Prime_Solv_GB",
        "s_psp_Prime_Covalent",
        "s_psp_Prime_vdW",
        "s_psp_Prime_Hbond",
        "s_psp_Prime_Lipo",
        "s_psp_Prime_Packing",
        "s_psp_Prime_SelfCont",
    ]
    all_keys = [
        "s_psp_res1",
        "s_psp_res2",
        "s_psp_Prime_Energy",
        "s_psp_Prime_Covalent",
        "s_psp_Prime_Coulomb",
        "s_psp_Prime_vdW",
        "s_psp_Prime_Solv_GB",
        "s_psp_Prime_Solv_SA",
        "s_psp_Prime_Lipo",
        "s_psp_Prime_Hbond",
        "s_psp_Prime_Packing",
        "s_psp_Prime_SelfCont",
        "s_psp_Prime_Entropy",
    ]
    ene_key = [
        "s_psp_Prime_Energy",
    ]

    def __init__(self, header, data):
        self.header = header
        self.data = data

    @classmethod
    def read(cls, file):
        with open(file) as f:
            lines = [l.strip() for l in f.readlines()]

        block_header_start = 0
        block_header_end = 0
        data_block_end = 0

        for i, line in enumerate(lines):
            if re.match(CT, line):
                block_header_start = i
            if (block_header_start and not block_header_end) and re.match(":::", line):
                block_header_end = i
            if block_header_start and block_header_end and re.match(":::", line):
                data_block_end = i

        header = lines[block_header_start + 1 : block_header_end]
        data = lines[block_header_end + 1 : data_block_end]
        return cls(header, [d.replace('"', "").split()[1:] for d in data])

    def get_inter_chain_inter(self):
        rows = []
        for i, row in enumerate(self.data):
            r1, r2 = row[0:2]
            c1, c2 = (r.split(":")[0] for r in (r1, r2))
            if c1 != c2:
                rows.append(i)

        return [self.data[i] for i in rows]

    def header2indices(self, keys=None):
        keys = keys or self.inter_ene_keys
        indices = []
        for key in keys:
            indices.append(self.header.index(key))
        return indices


def merge_res(data, indices):
    resinterdict = {}
    for i, row in enumerate(data):
        key1, key2 = row[0:2]
        e = sum(float(row[i]) for i in indices)
        resinterdict.setdefault(key1, 0)
        resinterdict.setdefault(key2, 0)
        resinterdict[key1] += e
        resinterdict[key2] += e
    return resinterdict


def key2asl(key):
    chain, name_num = key.split(":")
    resname, resnum = name_num.split("_")
    resnum = re.match(r"\d*", resnum).group()
    return f"r.n {resnum} and c.n {chain} and r.pt {resname}"


def color_atoms(data: dict, cmap):
    ws = maestro.workspace_get()

    vmin = min(data.values())
    vmax = max(data.values())
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm, cmap)

    for key, val in data.items():
        res_atoms = maestro.analyze.get_atoms_from_asl(ws, key2asl(key))
        color = [int(255 * c) for c in sm.to_rgba(val)[0:3]]
        for a in res_atoms:
            a.color = color

    maestro.workspace_set(ws)
    return ws


def parse_residue_interaction_energies_and_color_residues(filename, cmap="Greens_r"):
    ie = InteractionEnergies.read(filename)
    inter = ie.get_inter_chain_inter()
    indices = ie.header2indices()
    data = merge_res(inter, indices)
    color_atoms(data, cmap)
