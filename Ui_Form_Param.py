from PyQt4 import QtCore, QtGui

        
class Ui_Form(object):
    def __init__(self, parameters):
        self._parameters = parameters
    def setupUi(self, Form):
        Form.setObjectName( ("Form"))
        Form.resize(900, 300)
        self.layout = QtGui.QGridLayout()
        self.layout.setColumnStretch(0, 1)
        self.layout.setColumnStretch(3, 1)
        self.layout.setRowStretch(0, 1)
        self.layout.setRowStretch(3, 1)

        
        font = QtGui.QFont()
        font.setFamily(u"DejaVu Sans Mono") # police de Qt4
        font.setStyleHint(QtGui.QFont.Courier) # si la police est indisponible
        font.setPointSize(10)

        self.parameters = []
        for i, p in enumerate(self._parameters):
            if len(p)==3:
                label, value, comment = p
            elif len(p)==2:
                label, value = p
            label = QtGui.QLabel(label)
            edit = QtGui.QLineEdit()
            edit.setText(str(value))
            label.setFont(font)
            edit.setFont(font)
            self.parameters.append((label, edit))
            self.layout.addWidget(label, i, 1)
            self.layout.addWidget(edit, i, 2)
            
        Form.setLayout(self.layout)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def getParam(self):
        o = {}
        for l, e in self.parameters:
            o[str(l.text())] = eval(str(e.displayText()))
        return o
 
        
