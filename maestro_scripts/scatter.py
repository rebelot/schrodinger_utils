from schrodinger import maestro
from schrodinger.structutils import transform
from schrodinger.structutils import analyze
from schrodinger import structure
import numpy as np


def clash(st, structs):
    if not structs:
        return False
    com = analyze.center_of_mass(st)
    coms = np.array([analyze.center_of_mass(st) for st in structs])
    d = np.linalg.norm(coms - com, axis=1)
    r = analyze.radius_of_gyration(st)
    radii = np.array([analyze.radius_of_gyration(st) for st in structs])

    if np.any(d <= r + radii):
        return True
    return False


def scatter(boxsize, nmols):
    pt = maestro.project_table_get()
    structs = [row.structure for row in list(pt.selected_rows)]

    n = list(range(len(structs)))
    np.random.shuffle(n)

    if nmols < len(structs):
        idxs = np.random.choice(n, nmols, replace=False)
    else:
        idxs = n + np.random.choice(n, nmols - len(n)).tolist()

    new_structs = []

    for i in idxs:
        st = structs[i].copy()

        while clash(st, new_structs):
            transform.translate_centroid_to_origin(st)

            r = np.random.random(3) * 2 * np.pi
            transform.rotate_structure(st, *r)

            t = np.random.random(3) * boxsize
            transform.translate_structure(st, *t)

        new_structs.append(st)

    final = structure.create_new_structure()
    for st in new_structs:
        final = final.merge(st)

    pt.importStructure(final)
    pt.update()
    return final
