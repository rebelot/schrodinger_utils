__doc__ = """

Renumber ANY residue

"""
# Name: Renumber
# Command: pythonrun maestro_renumber.App

from schrodinger import maestro
from schrodinger.Qt import PyQt6
from PyQt6.QtWidgets import (
    QHBoxLayout,
    QLineEdit,
    QListWidget,
    QPushButton,
    QRadioButton,
    QVBoxLayout,
    QWidget,
)


class App(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.generateUI()
        self.st = None
        self.set_defaults()

    def initUI(self):
        self.number_le = QLineEdit()
        self.start_from_rb = QRadioButton("Start from:")
        self.fix_number_rb = QRadioButton("Fixed number:")
        self.renumber_btn = QPushButton("Renumber")
        self.report_list = QListWidget()
        self.report_list.itemDoubleClicked.connect(self.on_dclicked_list_item)
        self.renumber_btn.clicked.connect(self.on_clicked_renumber_btn)

    def generateUI(self):
        MainHBox = QHBoxLayout()
        MainHBox.addWidget(self.report_list)

        ModeVBox = QVBoxLayout()
        ModeVBox.addWidget(self.fix_number_rb)
        ModeVBox.addWidget(self.start_from_rb)
        ModeVBox.addWidget(self.number_le)
        ModeVBox.addWidget(self.renumber_btn)

        MainHBox.addLayout(ModeVBox)

        self.setLayout(MainHBox)
        # self.setGeometry()
        self.setWindowTitle("Notes")
        self.show()

    def set_defaults(self):
        self.st = maestro.workspace_get()
        self.get_res_list()

    def get_res_list(self):
        resdict = {}
        for res in self.st.residue:
            resdict.setdefault(res.pdbres, 0)
            resdict[res.pdbres] += 1
        for res, count in resdict.items():
            self.report_list.addItem(f"{res}: {count}")

    def on_dclicked_list_item(self):
        item = self.report_list.currentItem()
        resname = item.text().split(":")[0]
        maestro.command(f"workspaceselectionreplace r.pt {resname}")

    def on_clicked_renumber_btn(self):
        item = self.report_list.currentItem()
        resname = item.text().split(":")[0]
        n = int(self.number_le.text())
        residues = [r for r in self.st.residue if r.pdbres == resname]
        for i, res in enumerate(residues):
            res.resnum = i + n if self.start_from_rb.isChecked() else n
        maestro.workspace_set(self.st)
