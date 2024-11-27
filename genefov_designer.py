# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'genefov_designer_v4.ui'
#
# Created by: PyQt5 UI code generator 5.15.11
# irradiancia kWh/m2 (por el tamaño . cambio Nico)
#generacion  
# Pst unidades de potencia 
# Revisar unidades de Coef pot/T
#advertencias en rango de ingreso de datos 
# WARNING: Any manual changes made to this file will be lost when pyuic5 is
# run again.  Do not edit this file unless you know what you are doing.


from PyQt5 import QtCore, QtGui, QtWidgets
import sys
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import pandas as pd
from matplotlib.figure import Figure
import os
from PyQt5.QtWidgets import QMessageBox
from Genefov import Genefov

def resource_path(relative_path):
    """Obtiene la ruta absoluta a un recurso, funciona tanto en desarrollo como en el ejecutable."""
    try:
        # PyInstaller crea una carpeta temporal y almacena el nombre del ejecutable en _MEIPASS
        base_path = sys._MEIPASS
    except AttributeError:
        base_path = os.path.abspath(".")
        print(f"base:{base_path}")

    return os.path.join(base_path, relative_path)

def set_custom_font():
    font_path = resource_path("fuentes/Roboto-Regular.ttf")
    font_id = QtGui.QFontDatabase.addApplicationFont(font_path)
    if font_id != -1:
        font_family = QtGui.QFontDatabase.applicationFontFamilies(font_id)[0]
        custom_font = QtGui.QFont(font_family)
        QtWidgets.QApplication.setFont(custom_font)
    else:
        print("No se pudo cargar la fuente.")



class Ui_MainWindow(object):
        def setupUi(self, MainWindow):
                MainWindow.setObjectName("MainWindow")
                MainWindow.resize(1000, 750)
                sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Expanding)
                sizePolicy.setHorizontalStretch(0)
                sizePolicy.setVerticalStretch(0)
                sizePolicy.setHeightForWidth(MainWindow.sizePolicy().hasHeightForWidth())
                MainWindow.setSizePolicy(sizePolicy)
                self.centralwidget = QtWidgets.QWidget(MainWindow)
                self.centralwidget.setObjectName("centralwidget")
                self.horizontalLayout = QtWidgets.QHBoxLayout(self.centralwidget)
                self.horizontalLayout.setContentsMargins(2, 2, 2, 2)
                self.horizontalLayout.setSpacing(2)
                self.horizontalLayout.setObjectName("horizontalLayout")
                self.frame = QtWidgets.QFrame(self.centralwidget)
                self.frame.setMinimumSize(QtCore.QSize(500, 200))
                self.frame.setStyleSheet("background-color: rgb(150, 150, 150);")
                self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
                self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
                self.frame.setObjectName("frame")
                self.verticalLayout_2 = QtWidgets.QVBoxLayout(self.frame)
                self.verticalLayout_2.setContentsMargins(0, 0, 0, 0)
                self.verticalLayout_2.setSpacing(0)
                self.verticalLayout_2.setObjectName("verticalLayout_2")
                self.frame_ventana = QtWidgets.QFrame(self.frame)
                self.frame_ventana.setStyleSheet("background-color: rgb(167, 167, 167);\n"
                "border-color: rgb(255, 255, 255);")
                self.frame_ventana.setFrameShape(QtWidgets.QFrame.StyledPanel)
                self.frame_ventana.setFrameShadow(QtWidgets.QFrame.Raised)
                self.frame_ventana.setObjectName("frame_ventana")
                self.horizontalLayout_2 = QtWidgets.QHBoxLayout(self.frame_ventana)
                self.horizontalLayout_2.setContentsMargins(5, 5, 5, 5)
                self.horizontalLayout_2.setSpacing(5)
                self.horizontalLayout_2.setObjectName("horizontalLayout_2")
                self.frame_Configuracion = QtWidgets.QFrame(self.frame_ventana)
                self.frame_Configuracion.setStyleSheet("background-color: rgb(150, 150, 150);\n"
                "border-color: rgb(255, 255, 255);")
                self.frame_Configuracion.setFrameShape(QtWidgets.QFrame.StyledPanel)
                self.frame_Configuracion.setFrameShadow(QtWidgets.QFrame.Raised)
                self.frame_Configuracion.setObjectName("frame_Configuracion")
                self.verticalLayout = QtWidgets.QVBoxLayout(self.frame_Configuracion)
                self.verticalLayout.setObjectName("verticalLayout")
                self.groupBox_Sitios = QtWidgets.QGroupBox(self.frame_Configuracion)
                self.groupBox_Sitios.setMinimumSize(QtCore.QSize(221, 91))
                self.groupBox_Sitios.setMaximumSize(QtCore.QSize(16777215, 91))
                self.groupBox_Sitios.setBaseSize(QtCore.QSize(0, 0))
                font = QtGui.QFont()
                font.setBold(True)
                font.setWeight(75)
                self.groupBox_Sitios.setFont(font)
                self.groupBox_Sitios.setAlignment(QtCore.Qt.AlignLeading|QtCore.Qt.AlignLeft|QtCore.Qt.AlignVCenter)
                self.groupBox_Sitios.setObjectName("groupBox_Sitios")
                self.verticalLayout_4 = QtWidgets.QVBoxLayout(self.groupBox_Sitios)
                self.verticalLayout_4.setObjectName("verticalLayout_4")
                self.comboBox_Sitios = QtWidgets.QComboBox(self.groupBox_Sitios)
                self.comboBox_Sitios.setMaximumSize(QtCore.QSize(200, 30))
                self.comboBox_Sitios.setObjectName("comboBox_Sitios")
                self.comboBox_Sitios.addItem("")
                self.comboBox_Sitios.addItem("")
                self.comboBox_Sitios.addItem("")
                self.verticalLayout_4.addWidget(self.comboBox_Sitios)
                self.verticalLayout.addWidget(self.groupBox_Sitios)
                self.groupBox_orientacion = QtWidgets.QGroupBox(self.frame_Configuracion)
                self.groupBox_orientacion.setMaximumSize(QtCore.QSize(16777215, 200))
                font = QtGui.QFont()
                font.setBold(True)
                font.setWeight(75)
                self.groupBox_orientacion.setFont(font)
                self.groupBox_orientacion.setObjectName("groupBox_orientacion")
                self.gridLayout = QtWidgets.QGridLayout(self.groupBox_orientacion)
                self.gridLayout.setObjectName("gridLayout")
                self.label_2 = QtWidgets.QLabel(self.groupBox_orientacion)
                self.label_2.setObjectName("label_2")
                self.gridLayout.addWidget(self.label_2, 0, 0, 1, 1)
                self.lineEdit = QtWidgets.QLineEdit(self.groupBox_orientacion)
                self.lineEdit.setMaximumSize(QtCore.QSize(80, 30))
                self.lineEdit.setObjectName("lineEdit")
                self.gridLayout.addWidget(self.lineEdit, 0, 1, 1, 1)
                self.label_4 = QtWidgets.QLabel(self.groupBox_orientacion)
                self.label_4.setObjectName("label_4")
                self.gridLayout.addWidget(self.label_4, 0, 2, 1, 1)
                self.label_3 = QtWidgets.QLabel(self.groupBox_orientacion)
                self.label_3.setObjectName("label_3")
                self.gridLayout.addWidget(self.label_3, 1, 0, 1, 1)
                self.lineEdit_2 = QtWidgets.QLineEdit(self.groupBox_orientacion)
                self.lineEdit_2.setMaximumSize(QtCore.QSize(80, 30))
                self.lineEdit_2.setObjectName("lineEdit_2")
                self.gridLayout.addWidget(self.lineEdit_2, 1, 1, 1, 1)
                self.label_5 = QtWidgets.QLabel(self.groupBox_orientacion)
                self.label_5.setObjectName("label_5")
                self.gridLayout.addWidget(self.label_5, 1, 2, 1, 1)
                self.verticalLayout.addWidget(self.groupBox_orientacion)
                self.groupBox_caracteristicas = QtWidgets.QGroupBox(self.frame_Configuracion)
                font = QtGui.QFont()
                font.setBold(True)
                font.setWeight(75)
                self.groupBox_caracteristicas.setFont(font)
                self.groupBox_caracteristicas.setObjectName("groupBox_caracteristicas")
                self.gridLayout_2 = QtWidgets.QGridLayout(self.groupBox_caracteristicas)
                self.gridLayout_2.setObjectName("gridLayout_2")
                self.label_6 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_6.setObjectName("label_6")
                self.gridLayout_2.addWidget(self.label_6, 0, 0, 1, 1)
                self.lineEdit_3 = QtWidgets.QLineEdit(self.groupBox_caracteristicas)
                sizePolicy = QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
                sizePolicy.setHorizontalStretch(0)
                sizePolicy.setVerticalStretch(0)
                sizePolicy.setHeightForWidth(self.lineEdit_3.sizePolicy().hasHeightForWidth())
                self.lineEdit_3.setSizePolicy(sizePolicy)
                self.lineEdit_3.setMaximumSize(QtCore.QSize(80, 30))
                self.lineEdit_3.setAutoFillBackground(True)
                self.lineEdit_3.setObjectName("lineEdit_3")
                self.gridLayout_2.addWidget(self.lineEdit_3, 0, 1, 1, 1)
                self.label_7 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_7.setLayoutDirection(QtCore.Qt.LeftToRight)
                self.label_7.setObjectName("label_7")
                self.gridLayout_2.addWidget(self.label_7, 0, 2, 1, 1)
                self.label_8 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_8.setObjectName("label_8")
                self.gridLayout_2.addWidget(self.label_8, 1, 0, 1, 1)
                self.lineEdit_4 = QtWidgets.QLineEdit(self.groupBox_caracteristicas)
                self.lineEdit_4.setMaximumSize(QtCore.QSize(80, 30))
                self.lineEdit_4.setObjectName("lineEdit_4")
                self.gridLayout_2.addWidget(self.lineEdit_4, 1, 1, 1, 1)
                self.label_9 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_9.setObjectName("label_9")
                self.gridLayout_2.addWidget(self.label_9, 1, 2, 1, 1)
                self.label_10 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_10.setObjectName("label_10")
                self.gridLayout_2.addWidget(self.label_10, 2, 0, 1, 1)
                self.lineEdit_5 = QtWidgets.QLineEdit(self.groupBox_caracteristicas)
                self.lineEdit_5.setMaximumSize(QtCore.QSize(80, 30))
                self.lineEdit_5.setObjectName("lineEdit_5")
                self.gridLayout_2.addWidget(self.lineEdit_5, 2, 1, 1, 1)
                self.label_11 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_11.setObjectName("label_11")
                self.gridLayout_2.addWidget(self.label_11, 2, 2, 1, 1)
                self.label_16 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_16.setObjectName("label_16")
                self.gridLayout_2.addWidget(self.label_16, 3, 0, 1, 1)
                self.lineEdit_8 = QtWidgets.QLineEdit(self.groupBox_caracteristicas)
                self.lineEdit_8.setMaximumSize(QtCore.QSize(80, 30))
                self.lineEdit_8.setObjectName("lineEdit_8")
                self.gridLayout_2.addWidget(self.lineEdit_8, 3, 1, 1, 1)
                self.label_17 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_17.setObjectName("label_17")
                self.gridLayout_2.addWidget(self.label_17, 3, 2, 1, 1)
                self.label_12 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_12.setObjectName("label_12")
                self.gridLayout_2.addWidget(self.label_12, 4, 0, 1, 1)
                self.lineEdit_6 = QtWidgets.QLineEdit(self.groupBox_caracteristicas)
                self.lineEdit_6.setMaximumSize(QtCore.QSize(80, 30))
                self.lineEdit_6.setObjectName("lineEdit_6")
                self.gridLayout_2.addWidget(self.lineEdit_6, 4, 1, 1, 1)
                self.label_14 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_14.setObjectName("label_14")
                self.gridLayout_2.addWidget(self.label_14, 4, 2, 1, 1)
                self.label_13 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_13.setObjectName("label_13")
                self.gridLayout_2.addWidget(self.label_13, 5, 0, 1, 1)
                self.lineEdit_7 = QtWidgets.QLineEdit(self.groupBox_caracteristicas)
                self.lineEdit_7.setMaximumSize(QtCore.QSize(80, 30))
                self.lineEdit_7.setObjectName("lineEdit_7")
                self.gridLayout_2.addWidget(self.lineEdit_7, 5, 1, 1, 1)
                self.label_15 = QtWidgets.QLabel(self.groupBox_caracteristicas)
                self.label_15.setObjectName("label_15")
                self.gridLayout_2.addWidget(self.label_15, 5, 2, 1, 1)
                self.verticalLayout.addWidget(self.groupBox_caracteristicas)
                self.pushButton_calcular = QtWidgets.QPushButton(self.frame_Configuracion)
                self.pushButton_calcular.setMaximumSize(QtCore.QSize(16777215, 16777215))
                self.pushButton_calcular.setObjectName("pushButton_calcular")
                self.pushButton_calcular.clicked.connect(self.calcular_generacion)
                self.verticalLayout.addWidget(self.pushButton_calcular)
                self.horizontalLayout_2.addWidget(self.frame_Configuracion)
                self.frame_resultados = QtWidgets.QFrame(self.frame_ventana)
                self.frame_resultados.setStyleSheet("background-color: rgb(167, 167, 167);\n"
                "border-color: rgb(255, 255, 255);")
                self.frame_resultados.setFrameShape(QtWidgets.QFrame.StyledPanel)
                self.frame_resultados.setFrameShadow(QtWidgets.QFrame.Raised)
                self.frame_resultados.setObjectName("frame_resultados")
                self.verticalLayout_3 = QtWidgets.QVBoxLayout(self.frame_resultados)
                self.verticalLayout_3.setContentsMargins(5, 5, 5, 5)
                self.verticalLayout_3.setSpacing(5)
                self.verticalLayout_3.setObjectName("verticalLayout_3")
                self.frame_grafica = QtWidgets.QFrame(self.frame_resultados)
                self.frame_grafica.setStyleSheet("background-color: rgb(150, 150, 150);\n"
                "border-color: rgb(255, 255, 255);")
                self.frame_grafica.setFrameShape(QtWidgets.QFrame.StyledPanel)
                self.frame_grafica.setFrameShadow(QtWidgets.QFrame.Raised)
                self.frame_grafica.setObjectName("frame_grafica")

                #crear un layout horizontal
                self.horizontalLayout_5 = QtWidgets.QHBoxLayout(self.frame_grafica)
                self.horizontalLayout_5.setObjectName("horizontalLayout_5")
                ## Cavas aqui
                self.figure = plt.figure()
                self.canvas = FigureCanvas(self.figure)
                ## fin de canvas
                ## Agregar canvas
                self.horizontalLayout_5.addWidget(self.canvas)
                ##fin de horizontal layout 

                self.verticalLayout_3.addWidget(self.frame_grafica)
                self.frame_resultaods = QtWidgets.QFrame(self.frame_resultados)
                self.frame_resultaods.setStyleSheet("background-color: rgb(150, 150, 150);\n"
                "border-color: rgb(255, 255, 255);")
                self.frame_resultaods.setFrameShape(QtWidgets.QFrame.StyledPanel)
                self.frame_resultaods.setFrameShadow(QtWidgets.QFrame.Raised)
                self.frame_resultaods.setObjectName("frame_resultaods")
                self.horizontalLayout_3 = QtWidgets.QHBoxLayout(self.frame_resultaods)
                self.horizontalLayout_3.setObjectName("horizontalLayout_3")
                self.groupBox_resultados_irr = QtWidgets.QGroupBox(self.frame_resultaods)
                font = QtGui.QFont()
                font.setBold(True)
                font.setWeight(75)
                self.groupBox_resultados_irr.setFont(font)
                self.groupBox_resultados_irr.setObjectName("groupBox_resultados_irr")
                self.gridLayout_3 = QtWidgets.QGridLayout(self.groupBox_resultados_irr)
                self.gridLayout_3.setObjectName("gridLayout_3")
                self.label_18 = QtWidgets.QLabel(self.groupBox_resultados_irr)
                self.label_18.setObjectName("label_18")
                self.gridLayout_3.addWidget(self.label_18, 0, 0, 1, 1)
                self.textEdit = QtWidgets.QTextEdit(self.groupBox_resultados_irr)
                self.textEdit.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit.setObjectName("textEdit")
                self.gridLayout_3.addWidget(self.textEdit, 0, 1, 1, 1)
                self.label_22 = QtWidgets.QLabel(self.groupBox_resultados_irr)
                self.label_22.setObjectName("label_22")
                self.gridLayout_3.addWidget(self.label_22, 0, 2, 1, 1)
                self.label_19 = QtWidgets.QLabel(self.groupBox_resultados_irr)
                self.label_19.setObjectName("label_19")
                self.gridLayout_3.addWidget(self.label_19, 1, 0, 1, 1)
                self.textEdit_2 = QtWidgets.QTextEdit(self.groupBox_resultados_irr)
                self.textEdit_2.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit_2.setObjectName("textEdit_2")
                self.gridLayout_3.addWidget(self.textEdit_2, 1, 1, 1, 1)
                self.label_23 = QtWidgets.QLabel(self.groupBox_resultados_irr)
                self.label_23.setObjectName("label_23")
                self.gridLayout_3.addWidget(self.label_23, 1, 2, 1, 1)
                self.label_20 = QtWidgets.QLabel(self.groupBox_resultados_irr)
                self.label_20.setObjectName("label_20")
                self.gridLayout_3.addWidget(self.label_20, 2, 0, 1, 1)
                self.textEdit_3 = QtWidgets.QTextEdit(self.groupBox_resultados_irr)
                self.textEdit_3.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit_3.setObjectName("textEdit_3")
                self.gridLayout_3.addWidget(self.textEdit_3, 2, 1, 1, 1)
                self.label_24 = QtWidgets.QLabel(self.groupBox_resultados_irr)
                self.label_24.setObjectName("label_24")
                self.gridLayout_3.addWidget(self.label_24, 2, 2, 1, 1)
                self.label_21 = QtWidgets.QLabel(self.groupBox_resultados_irr)
                self.label_21.setObjectName("label_21")
                self.gridLayout_3.addWidget(self.label_21, 3, 0, 1, 1)
                self.textEdit_4 = QtWidgets.QTextEdit(self.groupBox_resultados_irr)
                self.textEdit_4.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit_4.setObjectName("textEdit_4")
                self.gridLayout_3.addWidget(self.textEdit_4, 3, 1, 1, 1)
                self.label_25 = QtWidgets.QLabel(self.groupBox_resultados_irr)
                self.label_25.setObjectName("label_25")
                self.gridLayout_3.addWidget(self.label_25, 3, 2, 1, 1)
                self.horizontalLayout_3.addWidget(self.groupBox_resultados_irr)
                self.groupBox_ResultadosGeneracion = QtWidgets.QGroupBox(self.frame_resultaods)
                font = QtGui.QFont()
                font.setBold(True)
                font.setWeight(75)
                self.groupBox_ResultadosGeneracion.setFont(font)
                self.groupBox_ResultadosGeneracion.setObjectName("groupBox_ResultadosGeneracion")
                self.gridLayout_4 = QtWidgets.QGridLayout(self.groupBox_ResultadosGeneracion)
                self.gridLayout_4.setObjectName("gridLayout_4")
                self.label_27 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_27.setObjectName("label_27")
                self.gridLayout_4.addWidget(self.label_27, 1, 0, 1, 1)
                self.label_30 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_30.setObjectName("label_30")
                self.gridLayout_4.addWidget(self.label_30, 4, 0, 1, 1)
                self.label_32 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_32.setObjectName("label_32")
                self.gridLayout_4.addWidget(self.label_32, 1, 2, 1, 1)
                self.label_33 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_33.setObjectName("label_33")
                self.gridLayout_4.addWidget(self.label_33, 2, 2, 1, 1)
                self.label_29 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_29.setObjectName("label_29")
                self.gridLayout_4.addWidget(self.label_29, 3, 0, 1, 1)
                self.label_26 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_26.setObjectName("label_26")
                self.gridLayout_4.addWidget(self.label_26, 0, 0, 1, 1)
                self.label_31 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_31.setObjectName("label_31")
                self.gridLayout_4.addWidget(self.label_31, 0, 2, 1, 1)
                self.label_28 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_28.setObjectName("label_28")
                self.gridLayout_4.addWidget(self.label_28, 2, 0, 1, 1)
                self.label_35 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_35.setObjectName("label_35")
                self.gridLayout_4.addWidget(self.label_35, 4, 2, 1, 1)
                self.label_34 = QtWidgets.QLabel(self.groupBox_ResultadosGeneracion)
                self.label_34.setObjectName("label_34")
                self.gridLayout_4.addWidget(self.label_34, 3, 2, 1, 1)
                self.textEdit_5 = QtWidgets.QTextEdit(self.groupBox_ResultadosGeneracion)
                self.textEdit_5.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit_5.setObjectName("textEdit_5")
                self.gridLayout_4.addWidget(self.textEdit_5, 0, 1, 1, 1)
                self.textEdit_6 = QtWidgets.QTextEdit(self.groupBox_ResultadosGeneracion)
                self.textEdit_6.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit_6.setObjectName("textEdit_6")
                self.gridLayout_4.addWidget(self.textEdit_6, 1, 1, 1, 1)
                self.textEdit_7 = QtWidgets.QTextEdit(self.groupBox_ResultadosGeneracion)
                self.textEdit_7.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit_7.setObjectName("textEdit_7")
                self.gridLayout_4.addWidget(self.textEdit_7, 2, 1, 1, 1)
                self.textEdit_8 = QtWidgets.QTextEdit(self.groupBox_ResultadosGeneracion)
                self.textEdit_8.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit_8.setObjectName("textEdit_8")
                self.gridLayout_4.addWidget(self.textEdit_8, 3, 1, 1, 1)
                self.textEdit_9 = QtWidgets.QTextEdit(self.groupBox_ResultadosGeneracion)
                self.textEdit_9.setMaximumSize(QtCore.QSize(110, 40))
                self.textEdit_9.setObjectName("textEdit_9")
                self.gridLayout_4.addWidget(self.textEdit_9, 4, 1, 1, 1)
                self.horizontalLayout_3.addWidget(self.groupBox_ResultadosGeneracion)
                self.pushButton_exportar = QtWidgets.QPushButton('Exportar Resultados a CSV',self.frame_resultaods)
                self.pushButton_exportar.setObjectName("pushButton_exportar")

                #self.exportar_button.setFont(font)
                self.pushButton_exportar.clicked.connect(self.exportar_resultados)

                self.horizontalLayout_3.addWidget(self.pushButton_exportar)
                self.verticalLayout_3.addWidget(self.frame_resultaods)
                self.verticalLayout_3.setStretch(0, 3)
                self.horizontalLayout_2.addWidget(self.frame_resultados)
                self.horizontalLayout_2.setStretch(0, 1)
                self.horizontalLayout_2.setStretch(1, 3)
                self.verticalLayout_2.addWidget(self.frame_ventana)
                self.frame_footer = QtWidgets.QFrame(self.frame)
                self.frame_footer.setStyleSheet("background-color: rgb(167, 167, 167);\n"
                "border-color: rgb(255, 255, 255);")
                self.frame_footer.setFrameShape(QtWidgets.QFrame.StyledPanel)
                self.frame_footer.setFrameShadow(QtWidgets.QFrame.Raised)
                self.frame_footer.setObjectName("frame_footer")
                self.horizontalLayout_4 = QtWidgets.QHBoxLayout(self.frame_footer)
                self.horizontalLayout_4.setContentsMargins(0, 0, 0, 0)
                self.horizontalLayout_4.setSpacing(0)
                self.horizontalLayout_4.setObjectName("horizontalLayout_4")
                self.label_footer = QtWidgets.QLabel(self.frame_footer)
                self.label_footer.setStyleSheet("font: italic 8pt \"MS Sans Serif\";")
                self.label_footer.setAlignment(QtCore.Qt.AlignCenter)
                self.label_footer.setObjectName("label_footer")
                self.horizontalLayout_4.addWidget(self.label_footer)
                self.verticalLayout_2.addWidget(self.frame_footer)
                self.verticalLayout_2.setStretch(0, 30)
                self.verticalLayout_2.setStretch(1, 1)
                self.horizontalLayout.addWidget(self.frame)
                MainWindow.setCentralWidget(self.centralwidget)

                self.retranslateUi(MainWindow)
                QtCore.QMetaObject.connectSlotsByName(MainWindow)

        def retranslateUi(self, MainWindow):
                _translate = QtCore.QCoreApplication.translate
                MainWindow.setWindowTitle(_translate("MainWindow", "GENEFOV (0.3.0)"))
                self.groupBox_Sitios.setTitle(_translate("MainWindow", "Sitios Disponibles"))
                self.comboBox_Sitios.setItemText(0, _translate("MainWindow", "Salta"))
                self.comboBox_Sitios.setItemText(1, _translate("MainWindow", "ElRosal60"))
                self.comboBox_Sitios.setItemText(2, _translate("MainWindow", "Yuto"))
                self.groupBox_orientacion.setTitle(_translate("MainWindow", "Orientacion del Panel"))
                self.label_2.setText(_translate("MainWindow", "Orientacion"))
                self.lineEdit.setText(_translate("MainWindow", "180"))
                self.label_4.setText(_translate("MainWindow", "[ ° ]"))
                self.label_3.setText(_translate("MainWindow", "Inclinacion"))
                self.lineEdit_2.setText(_translate("MainWindow", "24"))
                self.label_5.setText(_translate("MainWindow", "[ ° ]"))
                self.groupBox_caracteristicas.setTitle(_translate("MainWindow", "Caracteristicas del Panel"))
                self.label_6.setText(_translate("MainWindow", "K_temp"))
                self.lineEdit_3.setText(_translate("MainWindow", "0.029"))
                self.label_7.setText(_translate("MainWindow", "[ °C/Wm2  ]"))
                self.label_8.setText(_translate("MainWindow", "Pfv-stc"))
                self.lineEdit_4.setText(_translate("MainWindow", "450"))
                self.label_9.setText(_translate("MainWindow", "[ W ]"))
                self.label_10.setText(_translate("MainWindow", "Coef. Δpot/T"))
                self.lineEdit_5.setText(_translate("MainWindow", "-0.36"))
                self.label_11.setText(_translate("MainWindow", "[ %/°C ]"))
                self.label_16.setText(_translate("MainWindow", "Area Panel"))
                self.lineEdit_8.setText(_translate("MainWindow", "2.21"))
                self.label_17.setText(_translate("MainWindow", "[ m2 ]"))
                self.label_12.setText(_translate("MainWindow", "Albedo"))
                self.lineEdit_6.setText(_translate("MainWindow", "0.2"))
                self.label_14.setText(_translate("MainWindow", "[ - ]"))
                self.label_13.setText(_translate("MainWindow", "N° Paneles"))
                self.lineEdit_7.setText(_translate("MainWindow", "1"))
                self.label_15.setText(_translate("MainWindow", "[ - ]"))
                self.pushButton_calcular.setText(_translate("MainWindow", "Calcular"))
                self.groupBox_resultados_irr.setTitle(_translate("MainWindow", "Irradiancia Solar"))
                self.label_18.setText(_translate("MainWindow", "GHI_total"))
                self.label_22.setText(_translate("MainWindow", "[ kWh/m2 ]"))
                self.label_19.setText(_translate("MainWindow", "DHI_total"))
                self.label_23.setText(_translate("MainWindow", "[ kWh/m2 ]"))
                self.label_20.setText(_translate("MainWindow", "DNI_total"))
                self.label_24.setText(_translate("MainWindow", "[ kWh/m2 ]"))
                self.label_21.setText(_translate("MainWindow", "T_promedio"))
                self.label_25.setText(_translate("MainWindow", " [ ° C ]"))
                self.groupBox_ResultadosGeneracion.setTitle(_translate("MainWindow", "Reultados de la Generacion "))
                self.label_27.setText(_translate("MainWindow", "Otoño"))
                self.label_30.setText(_translate("MainWindow", "Primavera"))
                self.label_32.setText(_translate("MainWindow", "[ kWh ]"))
                self.label_33.setText(_translate("MainWindow", "[ kWh ]"))
                self.label_29.setText(_translate("MainWindow", "Invierno"))
                self.label_26.setText(_translate("MainWindow", "Total"))
                self.label_31.setText(_translate("MainWindow", "[ kWh ]"))
                self.label_28.setText(_translate("MainWindow", "Verano"))
                self.label_35.setText(_translate("MainWindow", "[ kWh ]"))
                self.label_34.setText(_translate("MainWindow", "[ kWh ]"))
                self.pushButton_exportar.setText(_translate("MainWindow", "Exportar"))
                self.label_footer.setText(_translate("MainWindow", "Grupo de Estudio y Evaluacion del Recurso Solar - INENCO (UNSA-CONICET)"))            

        def calcular_generacion(self):
                self.sitio_seleccionado = self.comboBox_Sitios.currentText()
                print(f"sitio: {self.sitio_seleccionado}")
                valores_por_sitio = {
                        'Salta': {'latitud': -24.7289, 'longitud': -65.4098, 'altura': 1233, 'gtm': 0},
                        'ElRosal60': {'latitud': -24.3928, 'longitud': -65.7683, 'altura': 3348, 'gtm': 0},
                        'yu': {'latitud': -20, 'longitud': -60, 'altura': 1100, 'gtm': 0},
                }
                
                parametros_sitio = valores_por_sitio.get(self.sitio_seleccionado, {})
                print(parametros_sitio)
                # Asignar los valores a los atributos de la clase
                self.latitud = parametros_sitio.get('latitud', 0)
        
                self.longitud = parametros_sitio.get('longitud', 0)
                self.altura = parametros_sitio.get('altura', 0)
                self.gtm = parametros_sitio.get('gtm', 0)

                # Verificar si el archivo CSV existe
                archivo_csv = resource_path(f"sitios\\{self.sitio_seleccionado}.csv")  # Usar resource_path para obtener la ruta correcta
                print(f" archivo csv: {archivo_csv}")
                if not os.path.exists(archivo_csv):
                        # Mostrar un mensaje de advertencia si el archivo no se encuentra
                        QMessageBox.warning(self.centralwidget, "Archivo no encontrado", "No se encontró el archivo necesario")
                        return  # Salir de la función si el archivo no existe

                # Obtener los valores de los QLineEdit y aplicar valores por defecto si el campo está vacío o no es válido

                    # Validación de parámetros con rangos
                def validar_parametro(valor_str, default, rango, nombre_parametro):
                        try:
                                valor = float(valor_str) if valor_str else default
                                if not (rango[0] <= valor <= rango[1]):
                                        raise ValueError
                        except ValueError:
                                valor = default
                                QMessageBox.warning(
                                        self.centralwidget,
                                        f"Valor fuera de rango: {nombre_parametro}",
                                        f"El valor ingresado para '{nombre_parametro}' está fuera del rango válido "
                                        f"({rango[0]} a {rango[1]}). Se usará el valor por defecto: {default}."
                                )
                        return valor

                # Validar y asignar parámetros
                self.orientacion = validar_parametro(self.lineEdit.text(), 180, (-180, 180), "Orientación")
                self.inclinacion = validar_parametro(self.lineEdit_2.text(), 24, (0, 90), "Inclinación")
                self.K_temp = validar_parametro(self.lineEdit_3.text(), 0.029, (0.01, 0.1), "Coeficiente Térmico (K_temp)")
                self.albedo = validar_parametro(self.lineEdit_6.text(), 0.2, (0, 1), "Albedo")
                self.area_paneles = validar_parametro(self.lineEdit_8.text(), 2.21, (0.1, 10), "Área de Paneles")
                self.PST = validar_parametro(self.lineEdit_4.text(), 450, (100, 1000), "PST")
                self.cvptm = validar_parametro(self.lineEdit_5.text(), -0.36, (-1, 0), "Coeficiente CVPTM")
                self.cantidad = validar_parametro(self.lineEdit_7.text(), 1, (1, 100), "Cantidad")



                # Crear una instancia de Genefov
                genefov = Genefov(
                        self.latitud, self.longitud, self.altura, self.gtm, self.orientacion, self.inclinacion, 
                        self.K_temp, self.PST, self.cvptm, self.albedo, self.area_paneles, self.cantidad, 
                        self.sitio_seleccionado
                ) # cantidad, pstc, coef. 

                # Obtener los resultados
                self.dias, self.vec_suma, total = genefov.diaria3
                GHI, DNI, DNIh ,DHI, tempprom = genefov.clausura
                otono, invierno, primavera, verano = genefov.Iestacional
                # Mostrar los resultados
                self.textEdit_5.setText(str(round(total, 1)))
                self.textEdit.setText(str(round(GHI, 1)))
                self.textEdit_2.setText(str(round(DHI, 1)))
                self.textEdit_3.setText(str(round(DNI, 1)))
                self.textEdit_4.setText(str(round(tempprom, 1)))


                self.textEdit_6.setText(str(round(otono, 1)))
                self.textEdit_8.setText(str(round(invierno, 1)))
                self.textEdit_9.setText(str(round(primavera, 1)))
                self.textEdit_7.setText(str(round(verano, 1)))
                
                self.plotOnCanvas(self.dias, self.vec_suma)
        
        def plotOnCanvas(self, dias, vec_suma):
                # Límites de las estaciones en el hemisferio sur con colores correspondientes
                estaciones_límites = [
                        (264, 354, 'Primavera', 'green'),
                        (355, 79, 'Verano', 'yellow'),
                        (80, 171, 'Otoño', 'orange'),
                        (172, 263, 'Invierno', 'skyblue')
                ]

                # Limpiar el canvas antes de graficar
                self.figure.clear()
                ax = self.canvas.figure.add_subplot(111)
                ax.plot(dias, vec_suma, label='Generación diaria')

                # Agregar áreas sombreadas para las estaciones y texto

                for inicio, fin, estacion, color in estaciones_límites: 
                        if inicio < fin: 
                                ax.axvspan(inicio, fin, color=color, alpha=0.3, label=estacion) 
                                ax.text((inicio + fin) / 2, max(vec_suma) * 1, estacion, ha='center', va='center') 
                        else: 
                                ax.axvspan(inicio, 365, color=color, alpha=0.3, label=estacion) 
                                ax.axvspan(0, fin, color=color, alpha=0.3, label=estacion) 
                                if estacion == 'Verano': ax.text((0 + fin) / 2, max(vec_suma) * 1, estacion, ha='center', va='center')

                ax.set_title('Generación Fotovoltaica')
                ax.set_xlabel('Días')
                ax.set_ylabel('Generación (kWh)')
                ax.grid()

                self.canvas.draw()




        
        def exportar_resultados(self):
                # Verificar si los arrays están vacíos
                if self.dias is None or self.vec_suma is None:
                        print("No hay datos para exportar.")
                        return

                # Asegurarse de que no estén vacíos
                if len(self.dias) == 0 or len(self.vec_suma) == 0:
                        print("Los arrays están vacíos. No hay datos para exportar.")
                        return

                # Asegurarse de que ambos arrays tienen la misma longitud
                if len(self.dias) != len(self.vec_suma):
                        print("Los arrays 'dias' y 'vec_suma' no tienen la misma longitud.")
                        return

                # Crear un diccionario con los datos para exportar
                data = {
                        'Días': self.dias,
                        'Generación (kWh)': self.vec_suma,
                        #'Sitio': self.sitio_seleccionado,  # Añadir información del sitio
                        #'Latitud': self.latitud,            # Añadir latitud
                        #'Longitud': self.longitud,          # Añadir longitud
                        #'Altura': self.altura,              # Añadir altura
                        #'Inclinación': self.inclinacion,    # Añadir inclinación
                        #'Orientación': self.orientacion,    # Añadir orientación
                        #'K_temp': self.K_temp,     # Añadir temperatura
                        #'Albedo': self.albedo,               # Añadir albedo
                        #'Área Paneles': self.area_paneles    # Añadir área de paneles
                }

                df = pd.DataFrame(data)
                filename = 'resultados_generacion.csv'

                with open(filename, 'w') as f:
                        # Agregar comentarios con la información del sitio y los parámetros
                        f.write(f"# Sitio: {self.sitio_seleccionado}")
                        f.write(f" # Latitud: {self.latitud}")
                        f.write(f" # Longitud: {self.longitud}")
                        f.write(f" # Altura: {self.altura}\n")
                        f.write(f"# Inclinacion: {self.inclinacion}")
                        f.write(f" # Orientacion: {self.orientacion}")
                        f.write(f" # K_temp: {self.K_temp}")
                        f.write(f" # Albedo: {self.albedo}")
                        f.write(f" # Area Paneles: {self.area_paneles}\n")
                        f.write("\n")  # Línea en blanco antes de los datos

                df.to_csv(filename, mode='a', index=False)
                # Mostrar un mensaje de éxito
                QMessageBox.information(self.centralwidget, "Exportación exitosa", f"Los resultados han sido exportados a: {filename}")






if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    set_custom_font()
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    ui.label_footer.setStyleSheet("font-size: 14px; font-style: italic;")
    MainWindow.show()
    sys.exit(app.exec_())

