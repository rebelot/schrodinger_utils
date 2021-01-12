__doc__ = """

Rotate along specified axis

"""
#Name: Transform
#Command: pythonrun maestro_transform.App

import numpy as np
from schrodinger import maestro
from schrodinger.Qt.PyQt5.QtWidgets import (QWidget, QGroupBox,
                                            QCheckBox, QRadioButton,
                                            QGridLayout, QLabel, QPushButton,
                                            QLineEdit, QHBoxLayout,
                                            QVBoxLayout)
from schrodinger.Qt.PyQt5.QtCore import pyqtSlot

class App(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()
        self.generateUI()
        self.set_defaults()

    def initUI(self):
        self.translate_rb = QRadioButton("Translate", self)
        self.rotate_rb = QRadioButton("Rotate", self)
        self.rotate_rb.toggled.connect(self.on_toggled_rotate_rb)

        self.hhelp = QLabel("Help ?")
        self.hhelp.setToolTip("Suka")
        self.hhelp.setDisabled(True)

        self.self_cb = QCheckBox("Self", self)
        self.self_cb.toggled.connect(self.on_toggled_self_cb)

        self.ori_cb = QCheckBox("Self")
        self.ori_cb.toggled.connect(self.on_toggled_ori_cb)

        self.abs_cb = QCheckBox("Absolute")
        self.abs_cb.toggled.connect(self.on_toggled_abs_cb)

        self.aalabel = QLabel("Move:")
        self.oalabel = QLabel("Centered on:")
        self.ralabel = QLabel("With respect to:")
        self.aale = QLineEdit("all")
        self.oale = QLineEdit("all")
        self.rale = QLineEdit("all")
        self.aale.textChanged.connect(self.on_aale_changed_text)

        self.axlabel = QLabel("Axis:")
        self.xlabel = QLabel("X")
        self.ylabel = QLabel("Y")
        self.zlabel = QLabel("Z")
        self.xle = QLineEdit()
        self.yle = QLineEdit()
        self.zle = QLineEdit()

        self.cax_cb = QCheckBox("Define W axis")
        self.cax_label1 = QLabel("from: ")
        self.cax_label2 = QLabel("to: ")
        self.cax_val_label = QLabel("W")
        self.cax_le1 = QLineEdit()
        self.cax_le2 = QLineEdit()
        self.cax_val_le = QLineEdit()
        self.cax_cb.toggled.connect(self.on_toggled_cax_cb)

        self.basel = QPushButton("Get selected", self)
        self.bosel = QPushButton("Get selected", self)
        self.brsel = QPushButton("Get selected", self)
        self.bcaxfsel = QPushButton("Get selected", self)
        self.bcaxtsel = QPushButton("Get selected", self)
        self.basel.clicked.connect(self.on_clicked_basel)
        self.bosel.clicked.connect(self.on_clicked_bosel)
        self.brsel.clicked.connect(self.on_clicked_brsel)
        self.bcaxfsel.clicked.connect(self.on_clicked_bcaxfsel)
        self.bcaxtsel.clicked.connect(self.on_clicked_bcaxtsel)

        self.bundo = QPushButton("Undo", self)
        self.bundo.setToolTip("wtf?!")
        self.bundo.clicked.connect(lambda: maestro.command("undo "))

        self.tbutton = QPushButton('Transform', self)
        self.tbutton.clicked.connect(self.on_click_tbutton)

        self.bdraw = QPushButton('Draw Axes')
        self.bdel = QPushButton('Delete')
        self.bdraw.clicked.connect(self.on_clicked_bdraw)
        self.bdel.clicked.connect(self.on_clicked_bdel)

        self.bset = QPushButton('Set')
        self.bset.clicked.connect(self.on_clicked_bset)

    def axis_groupUI(self):
        # Axis Group {{{
        axis_group = QGroupBox("Axis")

        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.axlabel)
        hbox1.addWidget(self.xlabel)
        hbox1.addWidget(self.xle)
        hbox1.addWidget(self.ylabel)
        hbox1.addWidget(self.yle)
        hbox1.addWidget(self.zlabel)
        hbox1.addWidget(self.zle)
        hbox1.addWidget(self.cax_val_label)
        hbox1.addWidget(self.cax_val_le)

        hbox2 = QHBoxLayout()
        hbox2.addStretch(1)
        hbox2.addWidget(self.bdraw)
        hbox2.addWidget(self.bdel)
        hbox2.addWidget(self.bset)

        hbox3 = QHBoxLayout()
        hbox3.addWidget(self.cax_cb)
        hbox3.addStretch(1)

        hbox4 = QHBoxLayout()
        hbox4.addWidget(self.cax_label1)
        hbox4.addWidget(self.cax_le1)
        hbox4.addWidget(self.bcaxfsel)

        hbox5 = QHBoxLayout()
        hbox5.addWidget(self.cax_label2)
        hbox5.addWidget(self.cax_le2)
        hbox5.addWidget(self.bcaxtsel)

        vbox1 = QVBoxLayout()
        vbox1.addLayout(hbox1)
        vbox1.addLayout(hbox3)
        vbox1.addLayout(hbox4)
        vbox1.addLayout(hbox5)
        vbox1.addLayout(hbox2)

        axis_group.setLayout(vbox1)
        # }}}
        return axis_group

    def aro_groupUI(self):
        # ARO Group
        aro_group = QGroupBox("Atoms")

        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.aalabel)
        hbox1.addWidget(self.aale)
        hbox1.addWidget(self.basel)

        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.oalabel)
        hbox2.addWidget(self.oale)
        hbox2.addWidget(self.bosel)

        hbox3 = QHBoxLayout()
        hbox3.addWidget(self.ralabel)
        hbox3.addWidget(self.rale)
        hbox3.addWidget(self.brsel)

        hbox4 = QHBoxLayout()
        hbox4.addWidget(self.self_cb)
        hbox4.addWidget(self.abs_cb)
        hbox4.addStretch(1)

        hbox5 = QHBoxLayout()
        hbox5.addWidget(self.ori_cb)
        hbox5.addStretch(1)

        vbox1 = QVBoxLayout()
        vbox1.addLayout(hbox1)
        vbox1.addLayout(hbox2)
        vbox1.addLayout(hbox5)
        vbox1.addLayout(hbox3)
        vbox1.addLayout(hbox4)

        aro_group.setLayout(vbox1)
        # }}}
        return aro_group

    def transform_groupUI(self):
        # Transform Group {{{
        tra_group = QGroupBox("Transform")

        vbox1 = QVBoxLayout()
        vbox1.addWidget(self.translate_rb)
        vbox1.addWidget(self.rotate_rb)

        hbox1 = QHBoxLayout()
        hbox1.addLayout(vbox1)
        hbox1.addStretch(1)
        hbox1.addWidget(self.tbutton)

        tra_group.setLayout(hbox1)
        # }}}
        return tra_group

    def bottom_lineUI(self):
        # Bottom Line
        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.hhelp)
        hbox1.addStretch(1)
        hbox1.addWidget(self.bundo)
        # }}}
        return hbox1

    def generateUI(self):
        gridUI = QGridLayout()
        gridUI.setSpacing(10)
        gridUI.addWidget(self.axis_groupUI(), 1, 1)
        gridUI.addWidget(self.aro_groupUI(), 2, 1)
        gridUI.addWidget(self.transform_groupUI(), 3, 1)
        gridUI.addLayout(self.bottom_lineUI(), 4, 1)

        self.setLayout(gridUI)

        self.setGeometry(300, 300, 400, 150)
        self.setWindowTitle('Transform')
        self.show()

    def set_defaults(self):
        self.self_cb.toggle()
        self.translate_rb.toggle()
        self.ori_cb.toggle()
        self.xyzList = []
        self.on_toggled_cax_cb()

    @pyqtSlot()
    def on_click_tbutton(self):
        x = float(self.xle.text()) if self.xle.text() else 0
        y = float(self.yle.text()) if self.yle.text() else 0
        z = float(self.zle.text()) if self.zle.text() else 0
        w = float(self.cax_val_le.text()) if self.cax_val_le.text() else 0
        aasl = self.aale.text() if self.aale.text() else 'all'
        rasl = self.rale.text() if self.rale.text() else aasl
        oasl = self.oale.text() if self.oale.text() else rasl
        fasl = self.cax_le1.text() if self.cax_cb.isChecked() else ''
        tasl = self.cax_le2.text() if self.cax_cb.isChecked() else ''
        abso = self.abs_cb.isChecked()
        if self.translate_rb.isChecked():
            translate(x, y, z, w, atoms_asl=aasl, origin_asl=oasl,
                      reference_asl=rasl, absolute=abso,
                      from_asl=fasl, to_asl=tasl)
        if self.rotate_rb.isChecked():
            rotate(x, y, z, w, atoms_asl=aasl, reference_asl=rasl,
                   absolute=abso, from_asl=fasl, to_asl=tasl)

    @pyqtSlot()
    def on_toggled_self_cb(self):
        if self.self_cb.isChecked():
            self.rale.setDisabled(True)
            self.ralabel.setDisabled(True)
            self.brsel.setDisabled(True)
            self.rale.setText(self.aale.text())
        else:
            self.rale.setDisabled(False)
            self.ralabel.setDisabled(False)
            self.brsel.setDisabled(False)

    @pyqtSlot()
    def on_toggled_abs_cb(self):
        if self.abs_cb.isChecked():
            self.rale.setDisabled(True)
            self.ralabel.setDisabled(True)
            self.brsel.setDisabled(True)
            self.self_cb.setDisabled(True)
        else:
            self.self_cb.setDisabled(False)
            self.on_toggled_self_cb()

    @pyqtSlot()
    def on_toggled_rotate_rb(self):
        if self.rotate_rb.isChecked():
            self.oale.setDisabled(True)
            self.oalabel.setDisabled(True)
            self.bosel.setDisabled(True)
            self.ori_cb.setDisabled(True)
        else:
            self.oale.setDisabled(False)
            self.oalabel.setDisabled(False)
            self.bosel.setDisabled(False)
            self.ori_cb.setDisabled(False)
            self.on_toggled_ori_cb()

    @pyqtSlot()
    def on_toggled_ori_cb(self):
        if self.ori_cb.isChecked():
            self.oalabel.setDisabled(True)
            self.oale.setDisabled(True)
            self.bosel.setDisabled(True)
            self.oale.setText(self.aale.text())
        else:
            self.oalabel.setDisabled(False)
            self.oale.setDisabled(False)
            self.bosel.setDisabled(False)

    @pyqtSlot()
    def on_aale_changed_text(self):
        if self.self_cb.isChecked() and not self.abs_cb.isChecked():
            self.rale.setText(self.aale.text())
        if self.ori_cb.isChecked():
            self.oale.setText(self.aale.text())

    @pyqtSlot()
    def on_toggled_cax_cb(self):
        if self.cax_cb.isChecked():
            self.cax_label1.setDisabled(False)
            self.cax_label2.setDisabled(False)
            self.cax_le1.setDisabled(False)
            self.cax_le2.setDisabled(False)
            self.cax_val_label.setDisabled(False)
            self.cax_val_le.setDisabled(False)
            self.bcaxfsel.setDisabled(False)
            self.bcaxtsel.setDisabled(False)
        else:
            self.cax_label1.setDisabled(True)
            self.cax_label2.setDisabled(True)
            self.cax_le1.setDisabled(True)
            self.cax_le2.setDisabled(True)
            self.cax_val_label.setDisabled(True)
            self.cax_val_le.setDisabled(True)
            self.bcaxfsel.setDisabled(True)
            self.bcaxtsel.setDisabled(True)

    @pyqtSlot()
    def on_clicked_basel(self):
        sel = maestro.selected_atoms_get_asl()
        self.aale.setText(sel)

    @pyqtSlot()
    def on_clicked_brsel(self):
        sel = maestro.selected_atoms_get_asl()
        self.rale.setText(sel)

    @pyqtSlot()
    def on_clicked_bosel(self):
        sel = maestro.selected_atoms_get_asl()
        self.oale.setText(sel)

    @pyqtSlot()
    def on_clicked_bcaxfsel(self):
        sel = maestro.selected_atoms_get_asl()
        self.cax_le1.setText(sel)

    @pyqtSlot()
    def on_clicked_bcaxtsel(self):
        sel = maestro.selected_atoms_get_asl()
        self.cax_le2.setText(sel)

    @pyqtSlot()
    def on_clicked_bdraw(self):
        x = float(self.xle.text()) if self.xle.text() else 0
        y = float(self.yle.text()) if self.yle.text() else 0
        z = float(self.zle.text()) if self.zle.text() else 0
        handle = maestro.create_model_xyz_axes(x, y, z, 100, 100, 100, 1, 0, 0)
        self.xyzList.append(handle)

    @pyqtSlot()
    def on_clicked_bdel(self):
        if len(self.xyzList) > 0:
            maestro.remove_object(self.xyzList.pop())

    @pyqtSlot()
    def on_clicked_bset(self):
        ws = maestro.workspace_get()
        sel = list(maestro.selected_atoms_get())
        if len(sel) > 0:
            com = maestro.analyze.center_of_mass(ws, atom_indices=sel)
            self.xle.setText(str(com[0]))
            self.yle.setText(str(com[1]))
            self.zle.setText(str(com[2]))
        else:
            self.xle.setText("0")
            self.yle.setText("0")
            self.zle.setText("0")


class xyzGraphicalObj():
    def __init__(self):
        self.handle = ''
        self.isShown = False

    def draw(self, x, y, z):
        self.handle = maestro.create_model_xyz_axes(
            x, y, z, 100, 100, 100, 1, 0, 0)
        self.isShown = True

    def cancel(self):
        maestro.remove_object(self.handle)
        self.isShown = False


def translate(x, y, z, w=0, atoms_asl='all', reference_asl='', origin_asl='', absolute=False, to_asl='', from_asl=''):
    reference_asl = reference_asl if reference_asl else atoms_asl
    origin_asl = origin_asl if origin_asl else atoms_asl

    ws = maestro.workspace_get()
    atoms = maestro.analyze.get_atoms_from_asl(ws, atoms_asl)

    origin_atoms = maestro.analyze.get_atoms_from_asl(ws, origin_asl)
    origin_indexes = [a.index for a in origin_atoms]
    origin_com = maestro.analyze.center_of_mass(
        ws, atom_indices=origin_indexes)

    if absolute:
        reference_com = np.array([0, 0, 0])
    else:
        reference_atoms = maestro.analyze.get_atoms_from_asl(ws, reference_asl)
        reference_indexes = [a.index for a in reference_atoms]
        reference_com = maestro.analyze.center_of_mass(
            ws, atom_indices=reference_indexes)

    T = np.array([x, y, z])
    if from_asl:
        u = get_axis(ws, from_asl, to_asl)
        T = T + w * u

    maestro.command('beginundoblock ')

    for at in atoms:
        at.xyz = np.array(at.xyz) + T + (reference_com - origin_com)

    maestro.workspace_set(ws)

    maestro.command('endundoblock ')


def rotate(x, y, z, w=0, atoms_asl='all', reference_asl='', absolute=False, from_asl='', to_asl=''):
    x = x / 57.29578
    y = y / 57.29578
    z = z / 57.29578
    w = w / 57.29578

    reference_asl = reference_asl if reference_asl else atoms_asl

    ws = maestro.workspace_get()
    atoms = maestro.analyze.get_atoms_from_asl(ws, atoms_asl)

    if absolute:
        com = np.array([0, 0, 0])
    else:
        reference_atoms = maestro.analyze.get_atoms_from_asl(ws, reference_asl)
        reference_indexes = [a.index for a in reference_atoms]
        com = maestro.analyze.center_of_mass(
            ws, atom_indices=reference_indexes)

    R = get_rotation_matrix(x, y, z)

    if from_asl:
        u = get_axis(ws, from_asl, to_asl)
        R = get_rotation_matrix_along_axis(w, u).dot(R)

    maestro.command('beginundoblock ')

    for at in atoms:
        zero = np.array(at.xyz) - com
        at.xyz = zero.dot(R) + com

    maestro.workspace_set(ws)

    maestro.command('endundoblock ')


def get_rotation_matrix(x, y, z):
    cosx = np.cos(x)
    sinx = np.sin(x)
    cosy = np.cos(y)
    siny = np.sin(y)
    cosz = np.cos(z)
    sinz = np.sin(z)
    Rx = np.array([[1, 0, 0],
                   [0, cosx, -sinx],
                   [0, sinx, cosx]])
    Ry = np.array([[cosy, 0, siny],
                   [0, 1, 0],
                   [-siny, 0, cosy]])
    Rz = np.array([[cosz, -sinz, 0],
                   [sinz, cosz, 0],
                   [0, 0, 1]])
    R = Rz.dot(Ry).dot(Rx)

    return R


def get_rotation_matrix_along_axis(theta, u):
    ux, uy, uz = u
    R = np.array([[np.cos(theta) + ux**2 * (1 - np.cos(theta)),
                   ux * uy * (1 - np.cos(theta)) - uz * np.sin(theta),
                   ux * uz * (1 - np.cos(theta)) + uy * np.sin(theta)],
                  [uy * ux * (1 - np.cos(theta)) + uz * np.sin(theta),
                   np.cos(theta) + uy**2 * (1 - np.cos(theta)),
                   uy * uz * (1 - np.cos(theta)) - ux * np.sin(theta)],
                  [uz * ux * (1 - np.cos(theta)) - uy * np.sin(theta),
                   uz * uy * (1 - np.cos(theta)) + ux * np.sin(theta),
                   np.cos(theta) + uz**2 * (1 - np.cos(theta))]])
    return R


def get_unit_vector(a, b):
    b = b - a
    u = b / np.sqrt(np.sum(np.square(b)))
    return u


def get_axis(ws, asl_a, asl_b):
    atoms_a = maestro.analyze.get_atoms_from_asl(ws, asl_a)
    atoms_b = maestro.analyze.get_atoms_from_asl(ws, asl_b)
    index_a = [a.index for a in atoms_a]
    index_b = [a.index for a in atoms_b]
    com_a = maestro.analyze.center_of_mass(ws, atom_indices=index_a)
    com_b = maestro.analyze.center_of_mass(ws, atom_indices=index_b)

    u = get_unit_vector(com_a, com_b)
    return u
