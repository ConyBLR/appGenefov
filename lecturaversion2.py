#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:11:09 2024

@author: nicolas
"""
import math
import numpy as np
from Genefov import Genefov
import pylab as plt

g = Genefov(-24.7288, -65.4095, 1233,0, 180, 0, 0.029, 450, -0.99, 0.2, 2.21,1, 'datos') #Lat, Long, Alt, gmt, Orientacion, Beta, k_temp, Pf_STC, Coef var pot temp mod, albedo, area, cant

df = g.df # DataFrame
dias, diaria, total, k = g.diaria3 #numero de dias;Energiadiaria kwh ;total producido; k: plotea
GHI, DNI, DNIh ,DHI, tempprom = g.clausura #GHI global ; DNI; DNI sobre el plano horizontal; DHI; temperatura promedio
otono, invierno, primavera, verano = g.Iestacional #Totales para cada estación
