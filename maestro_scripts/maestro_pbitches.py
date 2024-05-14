__doc__ = """

Fuck with PBCs!

"""

# Name: PBitChes
# Command: pythonrun maestro_pbitches.App

import numpy as np
from schrodinger import maestro
from schrodinger.application.desmond.packages import topo
from schrodinger.application.desmond.packages.pfx import pfx_module
from schrodinger.application.matsci.nano import xtal
from schrodinger.infra.structure import PBC
from schrodinger.Qt import PyQt6
from PyQt6.QtWidgets import (
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QVBoxLayout,
    QWidget,
)

# FIXME: It does not work on cms files, probably because of s_ffio_ct_type equal to full system property,
# which is ignored by desmond msys. Possible workarounds:
# - read structure from workspace, and props from pt
# - strip s_ffio_ct_type property or change it to solute
# option 1 is probably the best, so we do not fuck around with props,
# since geometry has to be rebuilt anyways


class App(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()
        self.generateUI()
        self.set_defaults()

    def initUI(self):
        self.lattice_a_label = QLabel("a")
        self.lattice_b_label = QLabel("b")
        self.lattice_c_label = QLabel("c")
        self.lattice_alpha_label = QLabel("alpha")
        self.lattice_beta_label = QLabel("beta")
        self.lattice_gamma_label = QLabel("gamma")

        self.lattice_a_le = QLineEdit()
        self.lattice_b_le = QLineEdit()
        self.lattice_c_le = QLineEdit()
        self.lattice_alpha_le = QLineEdit()
        self.lattice_beta_le = QLineEdit()
        self.lattice_gamma_le = QLineEdit()
        self.lattice_sync_button = QPushButton("Get from chorus")
        self.lattice_sync_button.clicked.connect(self.on_clicked_lattice_sync_button)

        self.chorus_ax_le = QLineEdit()
        self.chorus_ay_le = QLineEdit()
        self.chorus_az_le = QLineEdit()
        self.chorus_bx_le = QLineEdit()
        self.chorus_by_le = QLineEdit()
        self.chorus_bz_le = QLineEdit()
        self.chorus_cx_le = QLineEdit()
        self.chorus_cy_le = QLineEdit()
        self.chorus_cz_le = QLineEdit()
        self.chorus_sync_button = QPushButton("Get from lattice")
        self.chorus_sync_button.clicked.connect(self.on_clicked_chorus_sync_button)

        self.center_box_rb = QRadioButton("Box centered on molecules")
        self.anchor_box_rb = QRadioButton("Anchor box on:")
        self.center_box_rb.toggled.connect(self.on_toggled_position_rb)
        self.anchor_box_rb.toggled.connect(self.on_toggled_position_rb)
        self.anchor_x_label = QLabel("X")
        self.anchor_y_label = QLabel("Y")
        self.anchor_z_label = QLabel("Z")
        self.anchor_x_le = QLineEdit()
        self.anchor_y_le = QLineEdit()
        self.anchor_z_le = QLineEdit()
        self.anchor_center_to_origin_button = QPushButton(
            "Box center to origin (Desmond)"
        )
        self.anchor_center_to_origin_button.clicked.connect(
            self.on_clicked_center_to_origin_button
        )

        self.set_prop_button = QPushButton("Set Properties")
        self.set_prop_button.clicked.connect(self.on_clicked_set_prop_button)

        self.reload_button = QPushButton("Reload Structure")
        self.reload_button.clicked.connect(self.on_clicked_reload_button)

        self.apply_pfx_button = QPushButton("Apply PeriodicFix")
        self.apply_pfx_button.clicked.connect(self.on_clicked_apply_pfx_button)

    def generateUI(self):
        """
        |-------------------------|
        |LatticeGroup_____________|
        |            |            |
        |ChorusGroup |PositionGroup
        |            |            |
        |            |            |
        |____________|____________|
        |ButtonHBox_______________|

        """
        LatticeGroup = QGroupBox("Lattice Properties")
        LatticeHBox = QHBoxLayout()
        LatticeHBox.addWidget(self.lattice_a_label)
        LatticeHBox.addWidget(self.lattice_a_le)
        LatticeHBox.addWidget(self.lattice_b_label)
        LatticeHBox.addWidget(self.lattice_b_le)
        LatticeHBox.addWidget(self.lattice_c_label)
        LatticeHBox.addWidget(self.lattice_c_le)
        LatticeHBox.addWidget(self.lattice_alpha_label)
        LatticeHBox.addWidget(self.lattice_alpha_le)
        LatticeHBox.addWidget(self.lattice_beta_label)
        LatticeHBox.addWidget(self.lattice_beta_le)
        LatticeHBox.addWidget(self.lattice_gamma_label)
        LatticeHBox.addWidget(self.lattice_gamma_le)
        LatticeVBox = QVBoxLayout()
        LatticeVBox.addLayout(LatticeHBox)
        LatticeVBox.addWidget(self.lattice_sync_button)
        LatticeGroup.setLayout(LatticeVBox)

        ChorusGroup = QGroupBox("Chorus Porperties")
        ChorusGrid = QGridLayout()
        ChorusGrid.addWidget(self.chorus_ax_le, 1, 1)
        ChorusGrid.addWidget(self.chorus_ay_le, 1, 2)
        ChorusGrid.addWidget(self.chorus_az_le, 1, 3)
        ChorusGrid.addWidget(self.chorus_bx_le, 2, 1)
        ChorusGrid.addWidget(self.chorus_by_le, 2, 2)
        ChorusGrid.addWidget(self.chorus_bz_le, 2, 3)
        ChorusGrid.addWidget(self.chorus_cx_le, 3, 1)
        ChorusGrid.addWidget(self.chorus_cy_le, 3, 2)
        ChorusGrid.addWidget(self.chorus_cz_le, 3, 3)
        ChorusVBox = QVBoxLayout()
        ChorusVBox.addLayout(ChorusGrid)
        ChorusVBox.addWidget(self.chorus_sync_button)
        ChorusGroup.setLayout(ChorusVBox)

        PositionHBox = QHBoxLayout()
        PositionHBox.addWidget(self.anchor_x_label)
        PositionHBox.addWidget(self.anchor_x_le)
        PositionHBox.addWidget(self.anchor_y_label)
        PositionHBox.addWidget(self.anchor_y_le)
        PositionHBox.addWidget(self.anchor_z_label)
        PositionHBox.addWidget(self.anchor_z_le)
        PositionVBox = QVBoxLayout()
        PositionVBox.addWidget(self.center_box_rb)
        PositionVBox.addWidget(self.anchor_box_rb)
        PositionVBox.addLayout(PositionHBox)
        PositionVBox.addWidget(self.anchor_center_to_origin_button)

        MidHBox = QHBoxLayout()
        MidHBox.addWidget(ChorusGroup)
        MidHBox.addStretch(1)
        MidHBox.addLayout(PositionVBox)

        ButtonHbox = QHBoxLayout()
        ButtonHbox.addWidget(self.reload_button)
        ButtonHbox.addWidget(self.set_prop_button)
        ButtonHbox.addWidget(self.apply_pfx_button)

        UILayout = QVBoxLayout()
        UILayout.addWidget(LatticeGroup)
        UILayout.addLayout(MidHBox)
        UILayout.addLayout(ButtonHbox)

        self.setLayout(UILayout)
        self.setGeometry(145, 215, 606, 283)
        self.setWindowTitle("PBitChes")
        self.show()

    def set_defaults(self):
        self.get_structure()
        try:
            lattice_props, chorus_props = self.read_props()
            lattice_props = [str(prop) for prop in lattice_props]
            chorus_props = [str(prop) for prop in chorus_props]
            self.lattice_a_le.setText(lattice_props[0])
            self.lattice_b_le.setText(lattice_props[1])
            self.lattice_c_le.setText(lattice_props[2])
            self.lattice_alpha_le.setText(lattice_props[3])
            self.lattice_beta_le.setText(lattice_props[4])
            self.lattice_gamma_le.setText(lattice_props[5])
            self.chorus_ax_le.setText(chorus_props[0])
            self.chorus_ay_le.setText(chorus_props[1])
            self.chorus_az_le.setText(chorus_props[2])
            self.chorus_bx_le.setText(chorus_props[3])
            self.chorus_by_le.setText(chorus_props[4])
            self.chorus_bz_le.setText(chorus_props[5])
            self.chorus_cx_le.setText(chorus_props[6])
            self.chorus_cy_le.setText(chorus_props[7])
            self.chorus_cz_le.setText(chorus_props[8])
            try:
                if self.row.property[xtal.PBC_POSITION_KEY].split("_")[0] == "anchor":
                    self.anchor_box_rb.setChecked(True)
                    _, x, y, z = self.row.property[xtal.PBC_POSITION_KEY].split("_")
                    self.anchor_x_le.setText(x)
                    self.anchor_y_le.setText(y)
                    self.anchor_z_le.setText(z)
                else:
                    self.center_box_rb.setChecked(True)
                    self.anchor_x_le.setText("0")
                    self.anchor_y_le.setText("0")
                    self.anchor_z_le.setText("0")
            except KeyError:
                self.center_box_rb.setChecked(True)
                self.anchor_x_le.setText("0")
                self.anchor_y_le.setText("0")
                self.anchor_z_le.setText("0")
        except KeyError:
            MBox = QMessageBox()
            MBox.setText(
                "The selected structure does not have PBC properties, create periodic system?"
            )
            MBox.setStandardButtons(QMessageBox.Ok)
            MBox.buttonClicked.connect(self.on_clicked_mbox_ok)
            MBox.setIcon(QMessageBox.Warning)
            MBox.exec_()

    def on_clicked_apply_pfx_button(self):
        self.read_props()
        apply_pfx(self.st, self.chorus_box)
        # self.row.setStructure(self.st)
        maestro.workspace_set(self.st)

    def on_clicked_reload_button(self):
        self.set_defaults()

    def on_clicked_mbox_ok(self):
        pbc = PBC(1, 1, 1)
        pbc.applyToStructure(self.row)
        # self.row.setStructure(self.st)
        self.set_defaults()

    def on_clicked_set_prop_button(self):
        lattice_properties = [
            self.lattice_a_le.text() or 0,
            self.lattice_b_le.text() or 0,
            self.lattice_c_le.text() or 0,
            self.lattice_alpha_le.text() or 90,
            self.lattice_beta_le.text() or 90,
            self.lattice_gamma_le.text() or 90,
        ]
        chorus_properties = [
            self.chorus_ax_le.text() or 0,
            self.chorus_ay_le.text() or 0,
            self.chorus_az_le.text() or 0,
            self.chorus_bx_le.text() or 0,
            self.chorus_by_le.text() or 0,
            self.chorus_bz_le.text() or 0,
            self.chorus_cx_le.text() or 0,
            self.chorus_cy_le.text() or 0,
            self.chorus_cz_le.text() or 0,
        ]
        lattice_properties = [float(prop) for prop in lattice_properties]
        chorus_properties = [float(prop) for prop in chorus_properties]
        # FIXME: row workaround does not work, have to reimplement set_props
        # xtal.set_lattice_properties(self.row, lattice_properties)
        # xtal.set_pbc_properties(self.row, chorus_properties)
        chorus_prop_keys = [
            "r_chorus_box_ax",
            "r_chorus_box_ay",
            "r_chorus_box_az",
            "r_chorus_box_bx",
            "r_chorus_box_by",
            "r_chorus_box_bz",
            "r_chorus_box_cx",
            "r_chorus_box_cy",
            "r_chorus_box_cz",
        ]
        lattice_prop_keys = [
            "r_pdb_PDB_CRYST1_a",
            "r_pdb_PDB_CRYST1_b",
            "r_pdb_PDB_CRYST1_c",
            "r_pdb_PDB_CRYST1_alpha",
            "r_pdb_PDB_CRYST1_beta",
            "r_pdb_PDB_CRYST1_gamma",
        ]
        space_group_prop_key = "s_pdb_PDB_CRYST1_Space_Group"

        for prop, key in zip(chorus_properties, chorus_prop_keys):
            self.row.property[key] = prop
        for prop, key in zip(lattice_properties, lattice_prop_keys):
            self.row.property[key] = prop
        self.row.property[space_group_prop_key] = "P 1"

        if self.anchor_box_rb.isChecked():
            self.row.property[xtal.PBC_POSITION_KEY] = xtal.ANCHOR_PBC_POSITION % (
                self.anchor_x_le.text(),
                self.anchor_y_le.text(),
                self.anchor_z_le.text(),
            )
        elif self.center_box_rb.isChecked():
            self.row.property.pop(xtal.PBC_POSITION_KEY, None)

        self.pt.update()

    def on_clicked_chorus_sync_button(self):
        lattice_properties = [
            self.lattice_a_le.text() or 0,
            self.lattice_b_le.text() or 0,
            self.lattice_c_le.text() or 0,
            self.lattice_alpha_le.text() or 90,
            self.lattice_beta_le.text() or 90,
            self.lattice_gamma_le.text() or 90,
        ]
        lattice_properties = [float(prop) for prop in lattice_properties]
        chorus_props = xtal.get_chorus_from_params(lattice_properties)
        chorus_props = [str(prop) for prop in chorus_props]
        self.chorus_ax_le.setText(chorus_props[0])
        self.chorus_ay_le.setText(chorus_props[1])
        self.chorus_az_le.setText(chorus_props[2])
        self.chorus_bx_le.setText(chorus_props[3])
        self.chorus_by_le.setText(chorus_props[4])
        self.chorus_bz_le.setText(chorus_props[5])
        self.chorus_cx_le.setText(chorus_props[6])
        self.chorus_cy_le.setText(chorus_props[7])
        self.chorus_cz_le.setText(chorus_props[8])

    def on_clicked_lattice_sync_button(self):
        chorus_properties = [
            self.chorus_ax_le.text() or 0,
            self.chorus_ay_le.text() or 0,
            self.chorus_az_le.text() or 0,
            self.chorus_bx_le.text() or 0,
            self.chorus_by_le.text() or 0,
            self.chorus_bz_le.text() or 0,
            self.chorus_cx_le.text() or 0,
            self.chorus_cy_le.text() or 0,
            self.chorus_cz_le.text() or 0,
        ]
        chorus_properties = [float(prop) for prop in chorus_properties]
        lattice_props = xtal.get_params_from_chorus(chorus_properties)
        lattice_props = [str(prop) for prop in lattice_props]
        self.lattice_a_le.setText(lattice_props[0])
        self.lattice_b_le.setText(lattice_props[1])
        self.lattice_c_le.setText(lattice_props[2])
        self.lattice_alpha_le.setText(lattice_props[3])
        self.lattice_beta_le.setText(lattice_props[4])
        self.lattice_gamma_le.setText(lattice_props[5])

    def on_clicked_center_to_origin_button(self):
        x, y, z = (
            self.lattice_a_le.text(),
            self.lattice_b_le.text(),
            self.lattice_c_le.text(),
        )
        x, y, z = -0.5 * float(x), -0.5 * float(y), -0.5 * float(z)
        self.anchor_x_le.setText(str(x))
        self.anchor_y_le.setText(str(y))
        self.anchor_z_le.setText(str(z))

    def on_toggled_position_rb(self):
        if self.anchor_box_rb.isChecked():
            self.anchor_x_le.setDisabled(False)
            self.anchor_y_le.setDisabled(False)
            self.anchor_z_le.setDisabled(False)
        else:
            self.anchor_x_le.setDisabled(True)
            self.anchor_y_le.setDisabled(True)
            self.anchor_z_le.setDisabled(True)

    def get_structure(self):
        pt = maestro.project_table_get()
        row = list(pt.selected_rows)[0]
        self.pt = pt
        self.row = row
        self.st = maestro.workspace_get()

    def read_props(self):
        lattice_props = xtal.get_lattice_param_properties(self.row)
        chorus_props = xtal.get_chorus_properties(self.row)
        self.chorus_box = np.array(chorus_props).reshape((3, 3))
        return lattice_props, chorus_props


def get_topo(st):
    sys = topo.msys.LoadMAE(
        buffer=st.writeToString("maestro").encode(), structure_only=True
    )
    return sys.glued_topology


def apply_pfx(st, box):
    pos = st.getXYZ()
    glued_topo = get_topo(st)
    pfx = pfx_module.Pfx(glued_topo)
    pfx.apply(pos, box)
    st.setXYZ(pos)


def stretch_box(st, amount):
    box = get_box(st)
    box[-1, -1] = amount
    xtal.set_pbc_properties(st, box.flatten())
    pos = st.getXYZ()
    pos[:, -1] += 0.5 * box[-1, -1]
    st.setXYZ(pos)
    apply_pfx(st)
