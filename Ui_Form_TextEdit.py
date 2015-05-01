from PyQt4 import QtCore, QtGui
import QTextEditCS 

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName( ("Form"))
        Form.resize(900, 900)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setMargin(0)
        self.gridLayout.setSpacing(0)
        self.gridLayout.setObjectName( ("gridLayout"))
        self.plainTextEdit = QTextEditCS.QTextEditCS("", Form)
        self.plainTextEdit.setObjectName( ("plainTextEdit"))
        self.gridLayout.addWidget(self.plainTextEdit, 0, 0, 0, 0)
        QtCore.QMetaObject.connectSlotsByName(Form)
