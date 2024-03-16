from os import path
import sys
from pop_cnv.ui import ui_popcnv_main_form
from PySide6.QtWidgets import QWidget, QFileDialog, QMessageBox
from PySide6.QtCore import QCoreApplication, QEventLoop


class GPopCNV(QWidget):
    def __init__(self):
        super(GPopCNV, self).__init__()
        self.ui = None
        self.graph_scene = None

        self.__init_ui()

    def __init_ui(self):
        self.ui = ui_popcnv_main_form.Ui_frmMain()
        self.ui.setupUi(self)

        # Disable controls which cannot be used before run popCNV
        self.ui.grpSettings.setEnabled(False)
        self.ui.btn_run.setEnabled(False)
        self.ui.grpOutput.setEnabled(False)
