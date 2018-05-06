from PyQt5.uic import loadUiType

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas,
                                                NavigationToolbar2QT as NavigationToolbar)

Ui_MainWindow, QMainWindow = loadUiType('GUI.ui')

# todo: get the plotting into a function
# todo: have a button to select all the filenames to be treated
# todo: have a text field where the fitting parameters will be written
# todo: have a list where all the filenames will be displayed



class Main(QMainWindow, Ui_MainWindow):
    def __init__(self):
        super(Main, self).__init__()
        self.setupUi(self)

    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)
        self.canvas.draw()

if __name__ == '__main__':
    import sys
    from PyQt5 import QtWidgets
    from RheoFCClass import Fitter
    import Settings
    import numpy as np

    settings = Settings.Settings()
    filename = 'CF_Sac50-3--0.csv'
    fit = Fitter(filename, settings)
    #fit.plot_error_graphs()

    fig1 = Figure()
    ax1 = fig1.add_subplot(211)
    #ax2f1 = fig1.add_subplot(212)
    ax1.set_title(fit.model)
    x = np.logspace(np.log10(fit.GP[0]), np.log10(fit.GP[-1]))
    y, yerr = fit.carr_uncertainty(x, *fit.params, *fit.param_errs)
    ax1.plot(fit.GP, fit.Eta, linewidth=0, marker='o')
    ax1.errorbar(x, y, yerr=yerr)
    ax1.set_xscale('log')
    ax1.set_xlabel('$\dot{\gamma}/s^{-1}$')
    ax1.set_yscale('log')
    ax1.set_ylabel('$\eta_0$')
    ax1.annotate(str(fit.nl_first_point + 1), (fit.GP[fit.nl_first_point], fit.Eta[fit.nl_first_point]),
                 color='red')
    if fit.nl_last_point == -1:
        ax1.annotate(str(len(fit.GP)), (fit.GP[fit.nl_last_point], fit.Eta[fit.nl_last_point]),
                     color='red')  # todo: check this function
    else:
        ax1.annotate(str(fit.nl_last_point), (fit.GP[fit.nl_last_point], fit.Eta[fit.nl_last_point]), color='red')

    y_l, yerr_l = np.ones(len(x)) * fit.int, np.ones(len(x)) * fit.int_err

    ax2 = fig1.add_subplot(212)
    ax2.set_title('Linear')
    ax2.plot(fit.GP, fit.Eta, linewidth=0, marker='o')
    ax2.errorbar(x, y_l, yerr=yerr_l)
    ax2.set_xscale('log')
    ax2.set_xlabel('$\dot{\gamma}/s^{-1}$')
    ax2.set_yscale('log')
    ax2.set_ylabel('$\eta_0$')
    ax2.annotate(str(fit.l_first_point + 1), (fit.GP[fit.l_first_point], fit.Eta[fit.l_first_point]), color='red')
    if fit.l_last_point == -1:
        ax2.annotate(str(len(fit.GP)), (fit.GP[fit.l_last_point], fit.Eta[fit.l_last_point]),
                     color='red')  # todo: check this function
    else:
        ax2.annotate(str(fit.l_last_point), (fit.GP[fit.l_last_point], fit.Eta[fit.l_last_point]), color='red')

    fig1.tight_layout()

    app = QtWidgets.QApplication(sys.argv)
    main = Main()
    main.addmpl(fig1)
    main.show()
    sys.exit(app.exec_())
