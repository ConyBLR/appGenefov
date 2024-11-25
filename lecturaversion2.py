#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 15:11:09 2024

@author: Nicolas Rivera
"""
#
from Genefov import Genefov


g = Genefov(-24.7288, -65.4095, 1233,0, 180, 24, 0.029, 450, -0.36, 0.2, 2.21,1, 'sitios\Salta') #Lat, Long, Alt (m), UTC, Orientacion (°), Beta, k_temp, Pf_STC (W), Coef var pot temp mod, albedo, area (m2), cantidad de paneles       

df = g.df # DataFrame
dias, diaria, total, k = g.diaria3 #Numero de dias;Energiadiaria Wh/m2 ;total producido kWh/m2; k: plotea
GHI, DNI, DNIh ,DHI, tempprom = g.clausura #GHI horaria total (Wh/m2) global ; DNI horaria total (Wh/m2); DNI sobre el plano horizontal horaria total (Wh/m2); DHI horaria total (Wh/m2); temperatura promedio (°c)
otono, invierno, primavera, verano = g.Iestacional #Totales para cada estación en Wh/m2