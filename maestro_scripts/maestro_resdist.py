__doc__ = """

Calculate per-residue distance across multiple images

"""
#Name: Per-residue distance
#Command: pythonrun maestro_resdist.App

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
from schrodinger import maestro
from schrodinger.Qt.PyQt5.QtCore import pyqtSlot
from schrodinger.Qt.PyQt5.QtWidgets import (QCheckBox, QGridLayout, QGroupBox,
                                            QHBoxLayout, QLabel, QLineEdit,
                                            QMainWindow, QPushButton,
                                            QRadioButton, QVBoxLayout, QWidget)
from schrodinger.structutils import analyze, measure


class App(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.generateUI()
        self.set_defaults()
        maestro.project_update_callback_add(self.maestro_callback)

    def initUI(self):
        self.asl_le = QLineEdit('protein')
        self.asl_lb = QLabel('Selection ASL:')
        self.main_b = QPushButton("Go!")
        self.main_b.clicked.connect(self.on_clicked_main_b)
        self.avg_b = QPushButton("Create average structure")
        self.avg_b.clicked.connect(self.on_clicked_avg_b)
        self.cmap_le = QLineEdit('summer')
        self.cmap_lb = QLabel('colormap')
        self.cmap_lb.setToolTip("""Possible values are:
                Accent, Accent_r, Blues, Blues_r, BrBG, BrBG_r, BuGn, BuGn_r,
                BuPu, BuPu_r, CMRmap, CMRmap_r, Dark2, Dark2_r, GnBu, GnBu_r,
                Greens, Greens_r, Greys, Greys_r, OrRd, OrRd_r, Oranges,
                Oranges_r, PRGn, PRGn_r, Paired, Paired_r, Pastel1, Pastel1_r,
                Pastel2, Pastel2_r, PiYG, PiYG_r, PuBu, PuBuGn, PuBuGn_r,
                PuBu_r, PuOr, PuOr_r, PuRd, PuRd_r, Purples, Purples_r, RdBu,
                RdBu_r, RdGy, RdGy_r, RdPu, RdPu_r, RdYlBu, RdYlBu_r, RdYlGn,
                RdYlGn_r, Reds, Reds_r, Set1, Set1_r, Set2, Set2_r, Set3,
                Set3_r, Spectral, Spectral_r, Wistia, Wistia_r, YlGn, YlGnBu,
                YlGnBu_r, YlGn_r, YlOrBr, YlOrBr_r, YlOrRd, YlOrRd_r, afmhot,
                afmhot_r, autumn, autumn_r, binary, binary_r, bone, bone_r,
                brg, brg_r, bwr, bwr_r, cividis, cividis_r, cool, cool_r,
                coolwarm, coolwarm_r, copper, copper_r, cubehelix, cubehelix_r,
                flag, flag_r, gist_earth, gist_earth_r, gist_gray, gist_gray_r,
                gist_heat, gist_heat_r, gist_ncar, gist_ncar_r, gist_rainbow,
                gist_rainbow_r, gist_stern, gist_stern_r, gist_yarg,
                gist_yarg_r, gnuplot, gnuplot2, gnuplot2_r, gnuplot_r, gray,
                gray_r, hot, hot_r, hsv, hsv_r, inferno, inferno_r, jet, jet_r,
                magma, magma_r, nipy_spectral, nipy_spectral_r, ocean, ocean_r,
                pink, pink_r, plasma, plasma_r, prism, prism_r, rainbow,
                rainbow_r, seismic, seismic_r, spring, spring_r, summer,
                summer_r, tab10, tab10_r, tab20, tab20_r, tab20b, tab20b_r,
                tab20c, tab20c_r, terrain, terrain_r, viridis, viridis_r,
                winter, winter_r
                """)
        self.status = QLabel('')
        self.samecb = QCheckBox('Structures are the same')

    def generateUI(self):
        # gridUI = QGridLayout()

        hbox1 = QHBoxLayout()
        hbox1.addWidget(self.asl_lb)
        hbox1.addWidget(self.asl_le)

        hbox2 = QHBoxLayout()
        hbox2.addWidget(self.main_b)
        hbox2.addWidget(self.avg_b)
        hbox2.addWidget(self.cmap_lb)
        hbox2.addWidget(self.cmap_le)

        vbox = QVBoxLayout()
        vbox.addLayout(hbox1)
        vbox.addLayout(hbox2)
        vbox.addWidget(self.samecb)
        vbox.addWidget(self.status)

        self.setLayout(vbox)
        # self.setGeometry()
        self.setWindowTitle('Per-residue distances')
        self.show()

    def set_defaults(self):
        self.avg_b.setDisabled(True)
        self.maestro_callback()
        self.samecb.setChecked(True)

    def maestro_callback(self):
        pt = maestro.project_table_get()
        nsel = len(list(pt.selected_rows))
        nsel = nsel if nsel else 'WARNING, 0 entries selected'
        self.status.setText(f'Selected rows: {nsel}')

    @pyqtSlot()
    def on_clicked_main_b(self):
        if self.samecb.isChecked():
            self.dist = do_the_stuff2(self.asl_le.text())
        else:
            self.dist = do_the_suff(self.asl_le.text())
        make_plot(self.dist, self.cmap_le.text())
        self.avg_b.setDisabled(False)

    @pyqtSlot()
    def on_clicked_avg_b(self):
        dist2_avg_ribbon(self.dist, self.asl_le.text(), self.cmap_le.text())


def do_the_suff(asl):
    """
    Compute the distances between the center of mass of residues
    across copies of the same structure.
    Each structure has to be in a separate entry,
    structures must have the same number of residues
    and residues must be in the same order.
    Structures must be already superposed.

    :asl: ASL string
    :returns: ndarray of shape (entries, residues)

    """
    pt = maestro.project_table_get()
    rows = list(pt.selected_rows)

    ecoms = []
    for e in rows:
        st = e.getStructure()
        atom_indices = analyze.evaluate_asl(st, asl)
        st = st.extract(atom_indices)
        coms = []
        for res in st.residue:
            rcom = analyze.center_of_mass(st, res.getAtomIndices())
            coms.append(rcom)
        ecoms.append(np.array(coms))
    ecoms = np.array(ecoms)

    avg = np.mean(ecoms, axis=0)
    dist = np.sqrt(np.sum(np.square(ecoms - avg), axis=2))
    return dist


def do_the_stuff2(asl):
    """Get the substructure specified by asl, get the mean xyz coords for each
    residue, compute the RMSD for each residue with respect to the mean
    positions.

    :asl: string, asl expression
    :returns: ndarray of shape (entries, residues)
    """

    pt = maestro.project_table_get()
    rows = list(pt.selected_rows)

    XYZ = []

    for e in rows:
        st = e.getStructure()
        atom_indices = analyze.evaluate_asl(st, asl)
        st = st.extract(atom_indices)
        XYZ.append(st.getXYZ())
        n_atoms_per_residue = [len(r.getAtomIndices()) for r in st.residue]

    XYZ = np.array(XYZ)
    meanXYZ = XYZ.mean(axis=0)
    D = np.sqrt(np.sum(np.square(XYZ - meanXYZ), axis=2))

    RMSD = []
    for e in range(len(rows)):
        a, b = 0, 0
        rmsd = []
        for i in n_atoms_per_residue:
            a, b = b, b + i
            r = np.sqrt(np.sum(np.square(D[e, a:b])) / i)
            rmsd.append(r)
        RMSD.append(rmsd)

    RMSD = np.array(RMSD)

    return RMSD


def dist2_avg_ribbon(dist, asl, cmap):
    pt = maestro.project_table_get()
    sts = [row.getStructure() for row in pt.selected_rows]
    sts_e = [st.extract(analyze.evaluate_asl(st, asl)) for st in sts]
    avg_st = analyze.get_average_structure(sts_e)
    avg_row = pt.importStructure(avg_st)
    avg_row.includeOnly()
    avg_row.title = 'AVG'
    pt.update()

    mean_dist = dist.mean(axis=0)
    vmin = mean_dist.min()
    vmax = mean_dist.max()
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm, cmap)
    map2ribbon(avg_st, asl, sm, mean_dist)


def make_plot(dist, cmap):
    pt = maestro.project_table_get()

    f = plt.figure(constrained_layout=True)
    gs = GridSpec(2, 2, figure=f, height_ratios=[
                  3/4, 1/4], width_ratios=[4/5, 1/5])
    ax1 = f.add_subplot(gs[0, 0])
    ax2 = f.add_subplot(gs[1, 0], sharex=ax1)
    ax3 = f.add_subplot(gs[0, 1])
    ax3.set_aspect(20)

    p1 = ax1.pcolormesh(dist, cmap=cmap)
    ax1.set_ylabel('entry_id')
    plt.setp(ax1.get_xticklabels(), visible=False)
    ax1.set_yticks(range(len(pt.selected_rows)))
    ax1.set_yticklabels(row.entry_id for row in pt.selected_rows)
    ax1.grid(b=True, axis='y')

    plt.colorbar(p1, ax3)
    ax3.set_ylabel('distance (Å) to average')

    # ax2.pcolormesh(dist.mean(axis=0).reshape(1, dist.shape[1]), cmap=cmap)
    mean_dist = dist.mean(axis=0)
    vmin = mean_dist.min()
    vmax = mean_dist.max()
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm, cmap)
    ax2.bar(range(dist.shape[1]), mean_dist, color=sm.to_rgba(mean_dist), align='edge', width=1)
    ax2.set_xlabel('residue n')
    ax2.set_xlim([0, dist.shape[1]])
    ax2.set_title('average')
    ax2.set_yticks([])
    plt.show()


def map2ribbon(st, asl, sm, data):
    """
    Filter data through a ScalarMappable and map values to ribbon color
    representation of st.residues specified by asl
    """
    st = st.extract(analyze.evaluate_asl(st, asl))
    for res, col in zip(st.residue, sm.to_rgba(data)):
        res = f"protein and res.num {res.resnum} and chain.name {res.chain}"
        col = [int(255*c) for c in col[:4]]
        command = f'ribbon scheme=constant red={col[0]} green={col[1]} blue={col[2]} {res}'
        maestro.command(command)
        maestro.redraw_request()
