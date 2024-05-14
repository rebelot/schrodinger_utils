__doc__ = """

Select interacting residues

"""
#Name: Interacting
#Command: pythonrun maestro_interacting.main

from schrodinger import structutils
from schrodinger import maestro
from schrodinger.protein.analysis import analyze
from schrodinger.structutils import interactions


def get_interacting(asl_g1, asl_g2=None, structure=None):
    # who does not define a structure when out of Maestro?!
    st = maestro.workspace_get() if not structure else structure
    g1 = list(analyze.get_atoms_from_asl(st, asl_g1))
    g2 = list(analyze.get_atoms_from_asl(st, asl_g2)) if asl_g2 else asl_g2

    saltbr = interactions.get_salt_bridges(
        st, group1=g1, group2=g2)

    hbond = analyze.hbond.get_hydrogen_bonds(st, atoms1=g1, atoms2=g2)

    pipi = interactions.find_pi_pi_interactions(st, atoms1=g1, atoms2=g2)
    picat = interactions.find_pi_cation_interactions(st, atoms1=g1, atoms2=g2)

    atoms = []
    for bond in saltbr:
        atoms.append(bond[0])
        atoms.append(bond[1])

    for bond in hbond:
        atoms.append(bond[0])
        atoms.append(bond[1])

    for bond in pipi:
        atoms.append(bond[0].index)
        atoms.append(bond[1].index)

    for bond in picat:
        atoms.append(bond[0].index)
        atoms.append(bond[1].index)

    in_maestro = True   # TODO: check if inside maestro
    if in_maestro and len(atoms) > 0:
        maestro.command('workspaceselectionadd fillres a.n ' +
                        ', '.join(str(atom.index) for atom in atoms))
    return atoms


def main():
    selected_atoms = maestro.selected_atoms_get_asl()
    if selected_atoms:
        get_interacting(selected_atoms)
    else:
        maestro.warning(
            'Select something or use:\n"pythonrun maestro_interacting.get_interacting \"group1\" \"group2\""\nfrom the command prompt')


if __name__ == "__main__":
    main()
