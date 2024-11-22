#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 09:08:33 2024

@author: nicolas
"""

import math
import pandas as pd
import numpy as np

class Genefov:
    
    def __init__(self, lat, lon, alt, gamma, beta, temp, panelfv, Coef, albedo, area, estacion):
        """
        Constructor de la clase Genefov.
        
        :param lat: Latitud
        :param lon: Longitud
        :param alt: Altitud
        :param gamma: Ángulo de orientación
        :param beta: Ángulo de inclinación
        :param temp: Temperatura
        :param panelfv: Tipo de panel fotovoltaico
        :param Coef: Coeficientes
        :param albedo: Albedo
        :param area: Área del panel
        :param estacion: Nombre del archivo CSV con datos meteorológicos
        """
        self.lat = lat
        self.lon = lon
        self.alt = alt
        self.gamma = gamma
        self.beta = beta
        self.temp = temp
        self.panelfv = panelfv
        self.Coef = Coef
        self.albedo = albedo
        self.area = area
        self.estacion = estacion
        
        # Inicializa el DataFrame
        self.df = pd.DataFrame()
        
        # Intenta cargar el archivo CSV
        try:
            self.df2 = pd.read_csv(f"{self.estacion}.csv", header=2)
            self.df["J"] = self.df2["Day"]
            self.df["H"] = self.df2["Hour"]
            self.df["Hdec"] = self.df2["Hour"] + self.df2["Minute"] / 60
            self.df["GHI"] = self.df2["GHI"]
            self.df["DNI"] = self.df2["DNI"]
            self.df["DHI"] = self.df2["DHI"]
            self.df["Temp"] = self.df2["Temperature"]
        except FileNotFoundError:
            print(f"Error: El archivo {self.estacion}.csv no se encuentra.")
        except KeyError as e:
            print(f"Error: Falta una columna en el archivo CSV - {e}")
        except Exception as e:
            print(f"Error al leer el archivo CSV: {e}")
