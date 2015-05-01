"""
Author: Raphael.
"""

__version__ = '0.4'

import sys, os
from PyQt4 import QtCore, QtGui

import traceback
from run_script import *

class stderr:
    log = []
    def __init__(self):
        stderr.log = []
    def write(self, message):
        stderr.log.append(message)

from Ui_Form_TextEdit import Ui_Form as Ui_Form
from Ui_Form_Plot import Ui_Form as Ui_Form2
from Ui_Form_Log import Ui_Form as Ui_Form3
from Ui_Form_Log import stdProxy
from Ui_Form_Param import Ui_Form as Ui_Form4

from syntax import *


__stdout__ = sys.stdout
__stderr__ = sys.stderr

DEBUG = False


   

sys.path.append(os.path.realpath('.'))

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName( "MainWindow" )
        MainWindow.resize(1024, 800)
        self.MainWindow = MainWindow
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName( "centralwidget" )
        self.verticalLayout = QtGui.QVBoxLayout(self.centralwidget)
        self.verticalLayout.setObjectName( "verticalLayout" )
        self.mdiArea = QtGui.QMdiArea(self.centralwidget)
        self.mdiArea.setVerticalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.mdiArea.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAsNeeded)
        self.mdiArea.setActivationOrder(QtGui.QMdiArea.CreationOrder)
        #self.mdiArea.setViewMode(QtGui.QMdiArea.TabbedView)
        #self.mdiArea.setTabsClosable(True)
        #self.mdiArea.setTabsMovable(True)
        self.mdiArea.setObjectName( "mdiArea" )
        self.verticalLayout.addWidget(self.mdiArea)
        self.verticalLayout.setMargin(0)
        self.verticalLayout.setSpacing(0)

        MainWindow.setCentralWidget(self.centralwidget)

        if 'CST' in globals():
            self.windows_cst = {}


        self.setupMenuAction(MainWindow)

        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName( "statusbar" )
        MainWindow.setStatusBar(self.statusbar)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)
        self.statusbar.showMessage("Ready")

        self.windows = {}
        self.context = {}
        self.current_path = '.'

        # record globals
        self._globals = [k for k, v in globals().iteritems()]

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(   "QTLayout %s"%__version__ )

    def setupMenuAction(self, MainWindow):
        self.menubar = QtGui.QMenuBar(MainWindow)
        MainWindow.setMenuBar(self.menubar)
        
        self.menuFileAction = QtGui.QMenu(self.menubar)
        self.menuFileAction.setTitle( "&File" )
        self.menuFileAction.addAction( "&New", self.new_file, "Ctrl+N" )
        self.menuFileAction.addAction( "&Open", self.load_file, "Ctrl+L" )
        self.menuFileAction.addAction( "&Save", self.save_file, "Ctrl+S" )
        self.menuFileAction.addAction( "&Save As", self.saveas_file, "Ctrl+Shift+S" )
        self.menuFileAction.addAction( "E&xit", QtGui.qApp.quit, "Ctrl+Q" )
        self.menubar.addAction(self.menuFileAction.menuAction())
        
        self.menuRunAction = QtGui.QMenu(self.menubar)
        self.menuRunAction.setTitle( "&Run" )
        self.menuRunAction.addAction( "&Run Module", self.run_script, "F5" )
        self.menuRunAction.addAction( "&Clear Cache", self.clear_cache, "F9" )
        self.menuRunAction.addAction( "&Parameters Box", self.extract_parameters, "F6" )
        self.menuRunAction.addAction( "&Plot", self.plot, "F7" )
        self.menubar.addAction(self.menuRunAction.menuAction())
        
        if 'gdsii' in globals():
            try:
                gdsii.Cell
                self.menuGdsAction = QtGui.QMenu(self.menubar)
                self.menuGdsAction.setTitle( "&Gdsii" )
                self.menuGdsAction.addAction( "&Extract", self.extract_gds, "" )
                self.menubar.addAction(self.menuGdsAction.menuAction())
            except:
                pass

        if 'CST' in globals():
            self.menuCstAction = QtGui.QMenu(self.menubar)
            self.menuCstAction.setTitle( "&CST" )
            self.menuCstAction.addAction( "&New", self.new_cst, "" )
            self.menuCstAction.addAction( "&Open", self.open_cst, "" )
            self.menuCstAction.addAction( "&Save As", self.saveas_cst, "" )
            self.menuCstAction.addAction( "&Run", self.extract_cst, "F8" )
            self.menubar.addAction(self.menuCstAction.menuAction())
        

        self.menuWindowsAction = QtGui.QMenu(self.menubar)
        self.menuWindowsAction.setTitle( "&Windows" )
        self.menuWindowsAction.addAction( "&Cascade", self.cascade, "F2"  )
        self.menubar.addAction(self.menuWindowsAction.menuAction())

        self.menuHelpAction = QtGui.QMenu(self.menubar)
        self.menuHelpAction.setTitle( "&Help" )
        self.menuHelpAction.addAction( "&About", self.on_about, "F1"  )
        self.menubar.addAction(self.menuHelpAction.menuAction())

            
            
    def cascade(self):
        self.mdiArea.cascadeSubWindows()
        

    def on_about(self):
        msg = __doc__
        QtGui.QMessageBox.about(None, "About the QTLayout", msg.strip())

    def new_file(self):
        widget = QtGui.QWidget()
        self.subwin_abq = Ui_Form()
        self.subwin_abq.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( 'Untitled' )
        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()
        self.windows[self.subwindow] = self.subwin_abq.plainTextEdit, './Untitled.py'
      

    def load_file(self):
        print self.current_path
        self.filename = QtGui.QFileDialog.getOpenFileName(None,
            'Open a data file', self.current_path, 'python files (*.py);;All Files (*.*)')
        if self.filename:
            self.statusbar.showMessage("Loaded " + self.filename)
        else:
            self.statusbar.showMessage("Failed to load a file...")
            return
        self.filename = str(self.filename)

        # add path for import module
        sys.path.append(os.path.dirname(self.filename))


        widget = QtGui.QWidget()
        self.subwin_abq = Ui_Form()
        self.subwin_abq.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( self.filename.split('/')[-1].split('\\')[-1] + ' - ' + self.filename )
        self.subwindow.resize(700, 500)

        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()

        f = open(self.filename)
        string = f.read()
        self.subwin_abq.plainTextEdit.setText(string)
        f.close()

        self.windows[self.subwindow] = self.subwin_abq.plainTextEdit, self.filename

        #
        self.extract_parameters()

        self.current_path = os.path.dirname(self.filename)
                        

    def extract_parameters(self):
        active_sub_window = self.mdiArea.activeSubWindow()
        if not active_sub_window :
            QtGui.QMessageBox.about(None, "Alert", 'No window is selected.')
            return
        if not active_sub_window in self.windows:
            QtGui.QMessageBox.about(None, "Alert", 'The selected window is not correct.')
            return
        string = self.windows[active_sub_window][0].edit.toPlainText()
        string = str(string)

        obj = RunScript(string, self.context).run()
        l1 = obj.globals['Parameter'].list
        if len(l1):
            self.subwin_parameters = self.create_subwin_parameter(l1)
            obj.globals['Parameter'].list = list()
        else:
            self.subwin_parameters = None 
        


    def saveas_file(self):
        active_sub_window = self.mdiArea.activeSubWindow()
        if not active_sub_window :
            QtGui.QMessageBox.about(None, "Alert", 'No window is selected.')
            return
        if not active_sub_window in self.windows:
            QtGui.QMessageBox.about(None, "Alert", 'The selected window is not correct.')
            return

        self.filename = QtGui.QFileDialog.getSaveFileName(None,
            'Save a data file', self.current_path, 'python files (*.py);;All Files (*.*)')
        if self.filename:
            self.statusbar.showMessage("Saved " + self.filename)
        else:
            self.statusbar.showMessage("Failed to save to a file...")
            return

        f = open(self.filename, 'w')
        #string = self.subwin_abq.plainTextEdit.edit.toPlainText()
        string = self.windows[active_sub_window][0].edit.toPlainText()
        string = str(string)
        f.write(string)
        f.close()

        active_sub_window.setWindowTitle( self.filename.split('/')[-1].split('\\')[-1] + ' - ' + self.filename )
        self.windows[active_sub_window] = self.windows[active_sub_window][0], self.filename

        self.current_path = os.path.dirname(self.filename)
        
                    
    def save_file(self):
        active_sub_window = self.mdiArea.activeSubWindow()
        if not active_sub_window :
            QtGui.QMessageBox.about(None, "Alert", 'No window is selected.')
            return
        if not active_sub_window in self.windows:
            QtGui.QMessageBox.about(None, "Alert", 'The selected window is not correct.')
            return

        self.filename = self.windows[active_sub_window][1]

        f = open(self.filename, 'w')
        #string = self.subwin_abq.plainTextEdit.edit.toPlainText()
        string = self.windows[active_sub_window][0].edit.toPlainText()
        string = str(string)
        f.write(string)
        f.close()

        active_sub_window.setWindowTitle( self.filename.split('/')[-1].split('\\')[-1] + ' - ' + self.filename )
        
        

    def close(self):
        pass

    def run_script(self):

        active_sub_window = self.mdiArea.activeSubWindow()
        if not active_sub_window :
            QtGui.QMessageBox.about(None, "Alert", 'No window is selected.')
            return
        if not active_sub_window in self.windows:
            QtGui.QMessageBox.about(None, "Alert", 'The selected window is not correct.')
            return
        string = self.windows[active_sub_window][0].edit.toPlainText()
        string = str(string)

        self.statusbar.showMessage("Processing...")
        if not DEBUG:
            sys.stdout = stdProxy(self.create_subwin_stdout())
            sys.stderr = stdProxy(self.create_subwin_stderr())
        if self.subwin_parameters:
            self.context = self.subwin_parameters.getParam()            
        obj = RunScript(string, self.context).run()
        obj.globals['Parameter'].list = list()
        sys.stdout = __stdout__
        sys.stderr = __stderr__

        
        self.statusbar.showMessage("Done.")
        


    def create_subwin_stdout(self):
        # Standard Output Redirection Window
        widget = QtGui.QWidget()
        self.subwin_abq_stdout = Ui_Form3()
        self.subwin_abq_stdout.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( 'Standard Output Redirection' )
        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()
        return self.subwin_abq_stdout.plainTextEdit.append
        
    def create_subwin_stderr(self):
        # Standard Error Redirection Window
        widget = QtGui.QWidget()
        self.subwin_abq_stderr = Ui_Form3()
        self.subwin_abq_stderr.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( 'Standard Error Redirection' )
        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()
        self.subwin_abq_stderr.plainTextEdit.setTextColor(QtGui.QColor( 0, 1, 0 ))
        def wrap(text):
            self.subwin_abq_stderr.plainTextEdit.append("<span style=\"color: red\">" + text.replace('\n', '<br>') + "</span>")
        return wrap


    def create_subwin_parameter(self, parameters=[]):
        widget = QtGui.QWidget()
        self.subwin_abq_param = Ui_Form4(parameters)
        self.subwin_abq_param.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( 'Parameters' )
        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()
        return self.subwin_abq_param
        
        
    def plot(self):
        widget = QtGui.QWidget()
        self.subwin_abq_plot = Ui_Form2()
        self.subwin_abq_plot.setupUi(widget)

        active_sub_window = self.mdiArea.activeSubWindow()
        if not active_sub_window :
            QtGui.QMessageBox.about(None, "Alert", 'No window is selected.')
            return
        if not active_sub_window in self.windows:
            QtGui.QMessageBox.about(None, "Alert", 'The selected window is not correct.')
            return
        string = self.windows[active_sub_window][0].edit.toPlainText()
        string = str(string)

        self.statusbar.showMessage("Processing...")
        if not DEBUG:
            sys.stdout = stdProxy(self.create_subwin_stdout())
        if self.subwin_parameters:
            self.context = self.subwin_parameters.getParam()            
        result = RunScript(string, self.context).run()
        result.globals['Parameter'].list = list()
        sys.stdout = __stdout__
        self.statusbar.showMessage("Done.")
        if result.traceback<>'':
            stderr = self.create_subwin_stderr()
            stderr(result.traceback)
            return

        series = []
        for k, v in result.globals.iteritems():
            if isinstance(v, Primitives) and k[0]<>'_':
                for primitive in v:
                    series.append( [k, 'fill', [pt.x for pt in primitive], [pt.y for pt in primitive]] )
            if isinstance(v, Paths) and k[0]<>'_':
                for path in v:
                    series.append( [k, 'line', [pt.x for pt in path], [pt.y for pt in path]] )
        self.subwin_abq_plot.plot(series)

        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( 'Plot' )

        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()

    def clear_cache(self):
        QtGui.QMessageBox.about(None, "Alert", 'Not yet implemented')

        

    def extract_gds(self):
        widget = QtGui.QWidget()
        self.subwin_abq_plot = Ui_Form2()
        self.subwin_abq_plot.setupUi(widget)

        active_sub_window = self.mdiArea.activeSubWindow()
        if not active_sub_window :
            QtGui.QMessageBox.about(None, "Alert", 'No window is selected.')
            return
        if not active_sub_window in self.windows:
            QtGui.QMessageBox.about(None, "Alert", 'The selected window is not correct.')
            return
        string = self.windows[active_sub_window][0].edit.toPlainText()
        string = str(string)

        self.statusbar.showMessage("Processing...")
        if not DEBUG:
            sys.stdout = stdProxy(self.create_subwin_stdout())
        if self.subwin_parameters:
            self.context = self.subwin_parameters.getParam()            
        result = RunScript(string, self.context).run()
        sys.stdout = __stdout__
        self.statusbar.showMessage("Done.")
        if result.traceback<>'':
            stderr = self.create_subwin_stderr()
            stderr(result.traceback)
            return

        series = []
        for k, v in result.globals.iteritems():
            if isinstance(v, Primitives) and k[0]<>'_':
                series.append((k, v))


        self.filename = QtGui.QFileDialog.getSaveFileName(None,
            'Save a data file', '.', 'python files (*.gds);;All Files (*.*)')
        if self.filename:
            self.statusbar.showMessage("Saved " + self.filename)
        else:
            self.statusbar.showMessage("Failed to save to a file...")
            return

        for name, layer in series:
            cell = gdsii.Cell(str(name))
            cell.append( layer, layer=1 )
        gdsii.export(str(self.filename))
         
    def new_cst(self):
        widget = QtGui.QWidget()
        self.subwin_abq = Ui_Form()
        self.subwin_abq.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( 'Untitled' )
        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()
        self.windows_cst[self.subwindow] = self.subwin_abq.plainTextEdit

    def open_cst(self):
        self.filename = QtGui.QFileDialog.getOpenFileName(None,
            'Open a data file', '.', 'python files (*.py);;All Files (*.*)')
        if self.filename:
            self.statusbar.showMessage("Loaded " + self.filename)
        else:
            self.statusbar.showMessage("Failed to load a file...")
            return
        self.filename = str(self.filename)

        # add path for import module
        sys.path.append(os.path.dirname(self.filename))

        widget = QtGui.QWidget()
        self.subwin_abq = Ui_Form()
        self.subwin_abq.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( self.filename.split('/')[-1].split('\\')[-1] + ' - ' + self.filename )
        self.subwindow.resize(700, 500)

        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()

        f = open(self.filename)
        string = f.read()
        self.subwin_abq.plainTextEdit.setText(string)
        f.close()

        self.windows_cst[self.subwindow] = self.subwin_abq.plainTextEdit


    def saveas_cst(self):
        active_sub_window = self.mdiArea.activeSubWindow()
        if not active_sub_window :
            QtGui.QMessageBox.about(None, "Alert", 'No window is selected.')
            return
        if not active_sub_window in self.windows_cst:
            QtGui.QMessageBox.about(None, "Alert", 'The selected window is not correct.')
            return

        self.filename = QtGui.QFileDialog.getSaveFileName(None,
            'Save a data file', '.', 'python files (*.py);;All Files (*.*)')
        if self.filename:
            self.statusbar.showMessage("Saved " + self.filename)
        else:
            self.statusbar.showMessage("Failed to save to a file...")
            return

        f = open(self.filename, 'w')
        string = self.windows_cst[active_sub_window].edit.toPlainText()
        string = str(string)
        f.write(string)
        f.close()

        active_sub_window.setWindowTitle( self.filename.split('/')[-1].split('\\')[-1] + ' - ' + self.filename )

    def extract_cst(self):
        active_sub_window = self.mdiArea.activeSubWindow()
        if not active_sub_window :
            QtGui.QMessageBox.about(None, "Alert", 'No window is selected.')
            return
        if not active_sub_window in self.windows_cst:
            QtGui.QMessageBox.about(None, "Alert", 'The selected window is not correct.')
            return
        string = self.windows_cst[active_sub_window].edit.toPlainText()
        string = str(string)

        self.statusbar.showMessage("Processing...")
        if not DEBUG:
            sys.stdout = stdProxy(self.create_subwin_stdout())
        if hasattr(self, 'subwin_parameters') and self.subwin_parameters:
            self.context = self.subwin_parameters.getParam()            
        obj = RunScript(string, self.context).run()
        obj.globals['Parameter'].list = list()
        sys.stdout = __stdout__

        series = []
        for k, v in obj.globals.iteritems():
            if isinstance(v, CST) and k[0]<>'_':
                break
        widget = QtGui.QWidget()
        self.subwin_abq = Ui_Form3()
        self.subwin_abq.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.mdiArea)
        self.subwindow.setWindowTitle( 'CST VBA SCript' )
        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)
        self.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()
        self.subwin_abq.plainTextEdit.append(str(v))

        self.subwin_abq.plainTextEdit.selectAll()


        self.statusbar.showMessage("Done.")
        


         
class QTLayout(QtGui.QMainWindow):
    def __init__(self, parent=None):
        super(QTLayout, self).__init__(parent)
        self.ui = Ui_MainWindow()
        self.setWindowIcon(QtGui.QIcon('LayoutCreator.ico'))        
        self.ui.setupUi(self)

    def Add_Subwindow(self):
        widget = QtGui.QWidget()
        self.subwin_abq = Ui_Form()
        self.subwin_abq.setupUi(widget)
        self.subwindow = QtGui.QMdiSubWindow(self.ui.mdiArea) 
        widget.setParent(self.subwindow)
        self.subwindow.setWidget(widget)  
        self.subwindow.setWindowTitle("QTLayout")
        self.ui.mdiArea.addSubWindow(self.subwindow)
        widget.show()
        self.subwindow.show()
        self.subwindow.widget().show()



def main():
    global app
    app = QtGui.QApplication(sys.argv)
    form = QTLayout()
    form.show()
    app.exec_()


if __name__ == "__main__":
    main()
    
