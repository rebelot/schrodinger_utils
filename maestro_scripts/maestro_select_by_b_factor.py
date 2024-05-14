__doc__ = """

Select residues with highest b-factor within selection.
Either return a list of residues which B-factor is above cutoff (cutoff),
or a number residues with the highest B-factor (highest).

"""
# Name: B-Select
# Command: pythonrun maestro_select_by_b_factor.main

import numpy as np
from schrodinger import get_maestro
from schrodinger.structutils import analyze

maestro = get_maestro()


def select(asl, cutoff=40, highest=False):
    cutoff = float(cutoff)
    highest = int(highest) if highest is not False else False
    st = maestro.workspace_get()
    atoms = analyze.get_atoms_from_asl(st, asl)
    residues = list(dict.fromkeys(a.getResidue() for a in atoms))
    residues = np.array(residues)
    B = np.array([r.temperature_factor for r in residues])
    res = residues[B >= cutoff] if not highest else residues[np.argsort(B)[-highest:]]

    if len(res) > 0:
        maestro.command(
            "workspaceselectionreplace " + " OR ".join(r.getAsl() for r in res)
        )
    else:
        maestro.warning("Empty selection. Atoms are cool.")

    return res


def main():
    selected_atoms = maestro.selected_atoms_get_asl()
    if selected_atoms:
        select(selected_atoms)
    else:
        maestro.warning(
            'Select something or use: "maestro_select_by_b_factor.select asl cutoff highest" in the command line.'
        )


if __name__ == "__main__":
    main()
