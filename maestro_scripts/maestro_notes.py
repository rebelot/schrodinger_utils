__doc__ = """

Write awesome notes associated to entries

"""
# Name: Notes
# Command: pythonrun maestro_notes.App

from schrodinger import maestro
from schrodinger.Qt import PyQt6
from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)


class App(QWidget):
    def __init__(self):
        super().__init__()
        self.initUI()
        self.generateUI()
        self.set_defaults()
        # maestro.project_update_callback_add(self.maestro_callback)

    def initUI(self):
        self.editor = QTextEdit()
        self.entrylist = QListWidget()
        self.entrylist.itemDoubleClicked.connect(self.on_dclicked_list_item)
        self.entrylist.itemClicked.connect(self.on_clicked_list_item)
        self.entrylist.itemSelectionChanged.connect(self.on_selected_list_items)
        self.save_button = QPushButton("Save")
        self.save_button.clicked.connect(self.on_clicked_save_note_button)
        self.list_label = QLabel("Notes:")

    def generateUI(self):
        MainHBox = QHBoxLayout()

        ListVBox = QVBoxLayout()
        ListVBox.addWidget(self.list_label)
        ListVBox.addWidget(self.entrylist)

        NoteVBox = QVBoxLayout()
        NoteVBox.addWidget(self.editor)
        NoteVBox.addWidget(self.save_button)

        MainHBox.addLayout(ListVBox)
        MainHBox.addLayout(NoteVBox)

        self.setLayout(MainHBox)
        # self.setGeometry()
        self.setWindowTitle("Notes")
        self.show()

    def set_defaults(self):
        self.pt = maestro.project_table_get()
        self.make_list_items()
        self.maestro_callback()

    def maestro_callback(self):
        pass

    def make_list_items(self):
        self.entrylist.clear()
        for row in self.pt.all_rows:
            try:
                if len(row.property["s_user_notes"]) > 0:
                    pt_row_list_item = QListWidgetItem()
                    pt_row_list_item.setText(
                        f"({row.property['s_m_entry_id']}) {row.property['s_m_title']}"
                    )
                    pt_row_list_item.setData(Qt.UserRole, row)
                    self.entrylist.addItem(pt_row_list_item)
            except KeyError:
                continue

    def on_dclicked_list_item(self):
        item = self.entrylist.currentItem()
        row = item.data(Qt.UserRole)
        self.pt.selectRows(entry_ids=[row.entry_id])
        note = row.property["s_user_notes"]
        self.editor.setText(note)

    def on_clicked_list_item(self):
        item = self.entrylist.currentItem()
        row = item.data(Qt.UserRole)
        # self.pt.selectRows(entry_ids=[row.entry_id])
        note = row.property["s_user_notes"]
        self.editor.setText(note)

    def on_selected_list_items(self):
        pass

    def on_clicked_save_note_button(self):
        for row in self.pt.selected_rows:
            try:
                row.property["s_user_notes"] = self.editor.toPlainText()
            except KeyError:
                row.property.setdefault("s_user_notes", self.editor.toPlainText())
        self.pt.update()
        self.make_list_items()
