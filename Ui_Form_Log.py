import sys
from PyQt4 import QtCore, QtGui

class stdProxy:
    def __init__(self, write_func):
        self.write_func = write_func
    def write(self, text):
        stripped_text = text.strip('\n')
        if len(stripped_text):
            self.write_func(stripped_text)
            QtCore.QCoreApplication.processEvents()


class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName( ("Form"))
        Form.resize(800, 640)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setMargin(0)
        self.gridLayout.setSpacing(0)
        self.gridLayout.setObjectName( ("gridLayout"))
        self.plainTextEdit = QtGui.QTextEdit("", Form)
        self.plainTextEdit.setObjectName( ("plainTextEdit"))
        self.gridLayout.addWidget(self.plainTextEdit, 0, 0, 0, 0)
        QtCore.QMetaObject.connectSlotsByName(Form)
        


