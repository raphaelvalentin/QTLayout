from PyQt4 import QtCore, QtGui

import matplotlib
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg as NavigationToolbar
from matplotlib.figure import Figure


class linestyles(list):
    def __init__(self):
        list.__init__(self)
        self.append( {'color':'green', 'hatch':"//////"} )
        self.append( {'color':'red', 'hatch':"\\\\\\\\\\\\"} )
        self.append( {'color':'blue', 'hatch':"||||"} )
        self.append( {'color':'cyan', 'hatch':"---"} )
        self.append( {'color':'black', 'hatch':"...."} )
    def __getitem__(self, indx):
        indx = indx % len(self)
        return list.__getitem__(self,indx)
linestyles = linestyles()

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName( ("Form"))
        Form.resize(500, 1000)
        self.main_frame = QtGui.QWidget()
        self.dpi = 100
        self.fig = Figure((6.0, 6.0), dpi=self.dpi, facecolor='white', edgecolor='white')
        self.canvas = FigureCanvas(self.fig)

        self.canvas.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        self.canvas.updateGeometry() 
        
        self.canvas.setParent(self.main_frame)
        self.axes = self.fig.add_subplot(111)
        self.mpl_toolbar = NavigationToolbar(self.canvas, self.main_frame)
        self.main_frame.setObjectName( ("matplotlib"))

        self.series_list_model = QtGui.QStandardItemModel()
        self.series_list_view = QtGui.QListView()
        self.series_list_view.setFixedWidth(100)
        self.series_list_view.setModel(self.series_list_model)

        self.series_list_view.clicked[QtCore.QModelIndex].connect(self.slot)
        

        #self.Form = Form

        layout1 = QtGui.QHBoxLayout()
        layout1.setMargin(0)
        layout1.setSpacing(0)
        layout1.addWidget(self.canvas)
        layout1.addWidget(self.series_list_view)

        layout2 = QtGui.QHBoxLayout()
        layout2.setMargin(0)
        layout2.setSpacing(0)
        layout2.addWidget(self.mpl_toolbar)

        
        layout = QtGui.QVBoxLayout(Form)
        layout.setMargin(0)
        layout.setSpacing(0)
        layout.addLayout(layout2)
        layout.addLayout(layout1)
        
        Form.setLayout(layout)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def plot(self, series):
        
        self.series = series
        self.axes.clear()        
        self.axes.grid(True)
        self.axes.set_aspect('equal', 'datalim')
        self.axes.set_xlabel('x-coordinate (um)')
        self.axes.set_ylabel('y-coordinate (um)')
        names = []
        i=0
        xmin, xmax, ymin, ymax = 1e300, -1e300, 1e300, -1e300
        self.lines = []
        for k, t, x, y in series:
            xmin, xmax, ymin, ymax = min(min(x), xmin), max(max(x), xmax), min(min(y), ymin), max(max(y), ymax)
            isplotted = False
            for ki, linei in self.lines:
                if k==ki:
                    isplotted = True
                    break
            if isplotted:
                if t=='fill':
                    polygon = self.axes.fill(x, y, color='None', hatch=linei[0].get_hatch(), edgecolor=linei[0].get_edgecolor())
                elif t=='line':
                    polygon = self.axes.plot(x, y, color='b')
            else:
                if t=='fill':
                    polygon = self.axes.fill(x, y, color='None', hatch=linestyles[i]['hatch'], edgecolor=linestyles[i]['color'])
                elif t=='line':
                    polygon = self.axes.plot(x, y, color='b')
                i+=1
            self.lines.append((k, polygon))
            #self.axes.plot(x, y,'+')
            names.append(k)
            #i+=1
        self.fill_series_list(names)
        self.axes.set_xlim(xmin-(xmax-xmin)*0.1, xmax+(xmax-xmin)*0.1)
        self.axes.set_ylim(ymin-(ymax-ymin)*0.1, ymax+(ymax-ymin)*0.1)
        
    def fill_series_list(self, names):
        self.series_list_model.clear()

        _a = []
        for name in names:
            if name in _a:
                continue
            item = QtGui.QStandardItem(name)
            item.setCheckState(QtCore.Qt.Checked)
            item.setCheckable(True)
            self.series_list_model.appendRow(item)
            _a.append(name)

    def slot(self, index):
        item = self.series_list_model.itemFromIndex(index)
       
        for row in range(self.series_list_model.rowCount()):
            model_index = self.series_list_model.index(row, 0)
            checked = self.series_list_model.data(model_index,
                QtCore.Qt.CheckStateRole) == QtCore.QVariant(QtCore.Qt.Checked)
            name = str(self.series_list_model.data(model_index).toString())

            if checked:
                for line in self.lines:
                    if line[0]==name:
                        line[1][0].set_visible(True)
            else:
                for line in self.lines:
                    if line[0]==name:
                        line[1][0].set_visible(False)
            self.canvas.draw()

                
