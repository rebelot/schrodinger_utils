__doc__ = """

Place a dummy atom at the center of mass of selected atoms.

"""

#Name: Center of Mass
#Command: pythonrun maestro_dummy_at_center_of_mass.dummy_center_of_mass

from schrodinger import maestro


def dummy_center_of_mass():
    """
    Place a dummy atom at the center of mass of selected atoms

    """

    # get objects in workspace as a single structure
    st = maestro.workspace_get()
    # get selected atom indices (n,)
    sel = maestro.selected_atoms_get()

    # get center of mass of selected atoms
    com = maestro.analyze.center_of_mass(st, atom_indices=list(sel))
    # create dummy atom at com coordinates
    dummySt = maestro.analyze.create_new_structure()
    dummySt.addAtom("P", com[0], com[1], com[2], atom_type=150)

    # add dummy atom to pt and include dummy atom in workspace
    pt = maestro.project_table_get()
    row = pt.importStructure(dummySt, name="COM")
    row.title = 'COM'
    pt.includeRows([int(row.entry_id)], exclude_others=False)

    # set sphere representation for dummy atom
    maestro.command('repatom rep=cpk entry.id ' + row.entry_id)


if __name__ == "__main__":
    dummy_center_of_mass()
