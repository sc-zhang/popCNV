import matplotlib

matplotlib.use("qt5agg")
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.backends.backend_pdf
import matplotlib.pyplot as plt
import bioplotz as bp
from numpy import log2


class PlotCanvas(FigureCanvas):

    def __init__(self, width=20, height=8, dpi=300):
        fig = plt.figure(figsize=(width, height), dpi=dpi, tight_layout=True)
        self.plt = plt
        super(PlotCanvas, self).__init__(fig)


class PlotContent:

    def __init__(self):
        self.figure_content = PlotCanvas()

    def gene_line_graph(self, data_db):
        self.figure_content.plt.clf()
        try:
            xticks = []
            xlabels = []
            is_first = True
            for smp in data_db:
                off_set = 0
                X = []
                Y = []
                data = []
                for chrn in sorted(data_db[smp]):
                    for idx in range(len(data_db[smp][chrn]["X"])):
                        x = data_db[smp][chrn]["X"][idx]
                        y = data_db[smp][chrn]["Y"][idx]
                        data.append([x + off_set, y])
                    if is_first:
                        xticks.append(off_set + max(data_db[smp][chrn]["X"]) / 2.)
                        xlabels.append(chrn)
                    off_set += max(data_db[smp][chrn]["X"])
                is_first = False
                X = [v[0] for v in sorted(data)]
                Y = [v[1] for v in sorted(data)]
                self.figure_content.plt.plot(X, Y, label=smp, lw=1)
            self.figure_content.plt.legend()
            self.figure_content.plt.xticks(xticks, xlabels, rotation=-45)
        except Exception as e:
            return str(e)
        return "Success"

    def gene_manhattan_graph(self, data_db, is_single=False):
        self.figure_content.plt.clf()
        try:
            if is_single:
                plot_data = {}
                for smp in data_db:
                    for chrn in data_db[smp]:
                        plot_data[chrn] = [data_db[smp][chrn]["X"], data_db[smp][chrn]["Y"]]
                bp.manhattan(plot_data, threshold=[log2(0.05), log2(0.01)], color=["steelblue", "darkorange"],
                             threshold_line_color=["red", "blue"], log_base=2, reverse=True)
            else:
                plot_data = {}
                nrow = len(data_db)
                idx = 1
                for smp in data_db:
                    self.figure_content.plt.subplot(nrow, 1, idx)
                    idx += 1
                    for chrn in data_db[smp]:
                        plot_data[chrn] = [data_db[smp][chrn]["X"], data_db[smp][chrn]["Y"]]
                    bp.manhattan(plot_data, threshold=[log2(0.05), log2(0.01)], color=["steelblue", "darkorange"],
                                 threshold_line_color=["red", "blue"], log_base=2, reverse=True)
        except Exception as e:
            return str(e)

        return "Success"
