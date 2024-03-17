from os import path, listdir
from multiprocessing import cpu_count
from pop_cnv.ui import ui_popcnv_main_form
from pop_cnv.io import loader
from pop_cnv.worker.gui_popcnv_worker import Worker
from PySide6.QtWidgets import QWidget, QFileDialog, QMessageBox
from PySide6.QtCore import QThread, Signal
from traceback import format_exc
import gzip


class GPopCNV(QWidget):
    __run_signal = Signal(list)

    def __init__(self):
        super(GPopCNV, self).__init__()
        self.ui = None
        self.graph_scene = None

        # variables for popCNV
        self.__grp_db = {}
        self.__genome_file = ""
        self.__mosdepth_dir = ""
        self.__seq_depth_file = ""
        self.__smp_file = ""
        self.__gene_list_file = ""
        self.__win_size = 0
        self.__wrkdir = ""
        self.__threads = 0
        self.__chr_set = set()

        self.__worker = Worker()
        self.__worker_thread = QThread()

        self.__init_ui()

    def __init_ui(self):
        self.ui = ui_popcnv_main_form.Ui_frmMain()
        self.ui.setupUi(self)

        # Disable controls which cannot be used before run popCNV
        self.ui.grpSettings.setEnabled(False)
        self.ui.btn_run.setEnabled(False)
        self.ui.grpOutput.setEnabled(False)

        # Set button connects
        self.ui.btn_load_genome.clicked.connect(self.__get_genome_path)
        self.ui.btn_load_seqdepth.clicked.connect(self.__get_seqdepth_path)
        self.ui.btn_load_mosdepth.clicked.connect(self.__get_mosdepth_path)
        self.ui.btn_load_smp_list.clicked.connect(self.__get_smplist_path)
        self.ui.btn_load_gene_list.clicked.connect(self.__get_genelist_path)
        self.ui.btn_load_wrkdir.clicked.connect(self.__get_work_path)

        self.ui.btn_check.clicked.connect(self.__check_data)
        self.ui.btn_run.clicked.connect(self.__run_popcnv)

    def __get_genome_path(self):
        genome_file_path = QFileDialog.getOpenFileName(self, "Load Genome")[0]
        if genome_file_path:
            self.ui.text_genome_file_path.setText(genome_file_path)

    def __get_mosdepth_path(self):
        mosdepth_path = QFileDialog.getExistingDirectory(self, "Load Mosdepth")
        if mosdepth_path:
            self.ui.text_mosdepth_folder_path.setText(mosdepth_path)

    def __get_seqdepth_path(self):
        seqdepth_path = QFileDialog.getOpenFileName(self, "Load Seq depth")[0]
        if seqdepth_path:
            self.ui.text_seq_depth_file_path.setText(seqdepth_path)

    def __get_smplist_path(self):
        smplist_path = QFileDialog.getOpenFileName(self, "Load sample list")[0]
        if smplist_path:
            self.ui.text_smp_list_file_path.setText(smplist_path)

    def __get_genelist_path(self):
        genelist_path = QFileDialog.getOpenFileName(self, "Load gene list")[0]
        if genelist_path:
            self.ui.text_gene_list_file_path.setText(genelist_path)

    def __get_work_path(self):
        work_path = QFileDialog.getExistingDirectory(self, "Select work directory")
        if work_path:
            self.ui.text_wrkdir.setText(work_path)

    def __notify_with_title(self, info=""):
        if info:
            self.setWindowTitle("popCNV - %s" % info)
        else:
            self.setWindowTitle("popCNV")

    def __check_data(self):
        self.__genome_file = self.ui.text_genome_file_path.text()
        self.__mosdepth_dir = self.ui.text_mosdepth_folder_path.text()
        self.__seq_depth_file = self.ui.text_seq_depth_file_path.text()
        self.__smp_file = self.ui.text_smp_list_file_path.text()
        self.__gene_list_file = self.ui.text_gene_list_file_path.text()

        if not self.__genome_file or not path.isfile(self.__genome_file):
            QMessageBox.critical(self, "Loading error", "Genome file not found")
            return

        if not self.__mosdepth_dir or not path.isdir(self.__mosdepth_dir):
            QMessageBox.critical(self, "Loading error", "Mosdepth folder not found")
            return

        if not self.__seq_depth_file or not path.isfile(self.__seq_depth_file):
            QMessageBox.critical(self, "Loading error", "Seq depth file not found")
            return

        if not self.__smp_file or not path.isfile(self.__smp_file):
            QMessageBox.critical(self, "Loading error", "Sample list not found")
            return

        if not self.__gene_list_file or not path.isfile(self.__gene_list_file):
            QMessageBox.critical(self, "Loading error", "Gene list not found")
            return

        # Set threads list
        threads_list = [str(_) for _ in range(1, cpu_count() + 1)]
        self.ui.cbox_threads.addItems(threads_list)
        self.ui.cbox_threads.setCurrentText(threads_list[-1])

        # Get chromosomes list
        self.__get_chr_ids()

        # Load sample list
        self.__notify_with_title("Loading samples")
        try:
            grp_loader = loader.GRPLoader()
            grp_loader.load(self.__smp_file)
            self.__grp_db = grp_loader.grp_db
            grp_set = set(self.__grp_db[_] for _ in self.__grp_db)
            self.ui.cbox_wild_group.addItems(grp_set)
            self.__notify_with_title()
        except Exception as e:
            QMessageBox.critical(self, "Loading failed", format_exc())
            self.__notify_with_title(" Failed with: %s" % repr(e))
            return

        # Enable controls
        self.ui.grpSettings.setEnabled(True)
        self.ui.btn_run.setEnabled(True)

    def __get_win(self):
        for fn in listdir(self.__mosdepth_dir):
            if not fn.endswith(".regions.bed.gz"):
                continue
            mos_file = path.join(self.__mosdepth_dir, fn)
            with gzip.open(mos_file, 'rt') as fin:
                for line in fin:
                    data = line.strip().split()
                    sp = int(data[1])
                    ep = int(data[2])
                    return ep - sp

    def __get_chr_ids(self):
        with open(self.__gene_list_file) as fin:
            for line in fin:
                data = line.strip().split()
                self.__chr_set.add(data[0])

    def __popcnv_completed(self, v):
        self.__worker_thread.quit()
        # Enable controls before popcnv finished
        self.ui.btn_run.setEnabled(True)
        self.ui.btn_check.setEnabled(True)
        self.ui.grpInput.setEnabled(True)
        self.ui.grpSettings.setEnabled(True)
        self.ui.grpOutput.setEnabled(True)

        # Set output options
        self.ui.cbox_pic_smp.addItems(self.__grp_db.keys())
        self.ui.cbox_table_smp.addItems(self.__grp_db.keys())
        self.ui.cbox_pic_chr.addItems(sorted(self.__chr_set))
        self.ui.cbox_table_chr.addItems(sorted(self.__chr_set))

    def __run_popcnv(self):
        # Check workdir available
        self.__wild_group = self.ui.cbox_wild_group.currentText()
        self.__wrkdir = self.ui.text_wrkdir.text()
        self.__threads = int(self.ui.cbox_threads.currentText())
        if not self.__wrkdir or not path.isdir(self.__wrkdir):
            QMessageBox.critical(self, "Running error", "Work directory not found")
            return

        # Disable controls before popcnv finished
        self.ui.btn_run.setEnabled(False)
        self.ui.btn_check.setEnabled(False)
        self.ui.grpInput.setEnabled(False)
        self.ui.grpSettings.setEnabled(False)

        self.__win_size = self.__get_win()
        self.__worker.progress.connect(self.__notify_with_title)
        self.__worker.completed.connect(self.__popcnv_completed)
        self.__worker.set_param(self,
                                self.__genome_file,
                                self.__mosdepth_dir,
                                self.__seq_depth_file,
                                self.__smp_file,
                                self.__gene_list_file,
                                self.__win_size,
                                self.__wrkdir,
                                self.__wild_group,
                                self.__threads)
        self.__run_signal.connect(self.__worker.run_popcnv)
        self.__worker.moveToThread(self.__worker_thread)
        self.__run_signal.emit(1)
        self.__worker_thread.start()
