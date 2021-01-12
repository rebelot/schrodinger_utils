__doc__ = """

Map data to atom colors

"""
#Name: Data to atom color
#Command: pythonrun maestro_applycscheme map2color

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

def map2color(st, data, cmap):
    """
    Filter data through a normalized ScalarMappable and map values to atom colors
    data are assumed to map to individual structure residues:
    make sure that the order of `data` and `st.residue` is consistent

    :st: maestro.Structure
    :data: 1D np.array of length st.residue
    :cmap: string (available matplotlib colormaps)
    """
    vmin = data.min()
    vmax = data.max()
    norm = Normalize(vmin=vmin, vmax=vmax)
    sm = ScalarMappable(norm, cmap)
    for res, col in zip(st.residue, sm.to_rgba(data)):
        col = [int(255*c) for c in col[:4]]
        command = f"coloratomrgb red={col[0]} green={col[1]} blue={col[2]} a.n {','.join(str(i) for i in res.getAtomIndices())}"
        maestro.command(command)
    maestro.redraw_request()
