__doc__ = """

Draw principal axes of inertia

"""

# Name: MOI
# Command: pythonrun maestro_moi.App

import numpy as np
from schrodinger import maestro
from schrodinger.graphics3d import arrow, common
from schrodinger.Qt import PyQt6
from PyQt6.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QWidget,
)
from schrodinger.structutils import analyze


class App(QWidget):
    def __init__(self):
        super().__init__()
        self.axes = None
        self.labels = None
        self.generate_ui()

    def generate_ui(self):
        gen_button = QPushButton("Generate from selection")
        gen_button.clicked.connect(self.on_clicked_gen_button)
        del_button = QPushButton("Hide")
        del_button.clicked.connect(self.on_clicked_del_button)
        layout = QHBoxLayout()
        layout.addWidget(gen_button)
        layout.addWidget(del_button)
        self.setLayout(layout)
        self.setWindowTitle("MOI")
        self.show()

    def on_clicked_gen_button(self):
        self.axes, self.labels = draw()

    def on_clicked_del_button(self):
        self.axes.hide()
        for lab in self.labels:
            maestro.remove_object(lab)


def draw():
    ws = maestro.workspace_get()
    atom_i = maestro.selected_atoms_get()
    st_sel = ws.extract(atom_i)
    com = np.array(analyze.center_of_mass(st_sel))
    moi, pax = analyze.calculate_principal_moments(struct=st_sel)
    pax = np.array(pax)

    arrow_grp = common.Group()
    lab_handles = []
    for ax, c, lab in zip(com - pax * 10, ("red", "green", "blue"), ("a", "b", "c")):
        arrow_grp.add(arrow.MaestroArrow(*ax, *com, color=c))
        lab_pos = (ax - com) / 2 + com + 0.5
        lab_handles.append(maestro.create_single_line_z_buffered_text(lab, *lab_pos))
    arrow_grp.show()
    return arrow_grp, lab_handles
