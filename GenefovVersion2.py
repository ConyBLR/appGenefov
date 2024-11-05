#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 14:16:30 2024

@author: nicolas
"""

import math
import pandas as pd
import numpy as np
import pylab as plt
class Genefov:
    
    def __init__(self, lat, lon, alt, gmt, orientacion, beta, ktemp, PST, cvptm, albedo, area, cant,estacion):
        """
        df2 es un dataframe que tiene los datos de GHI DNI DHI TMY
        """
        self.lat = lat
        self.gmt = gmt
        self.lon = lon
        self.alt = alt
        self.orientacion = orientacion
        self.beta = beta
        self.ktemp = ktemp
        self.cvptm = cvptm
        self.PST = PST
        self.albedo = albedo
        self.area = area
        #self.df2 = df2 
        self.cant = cant
        self.df = pd.DataFrame()
        self.estacion = estacion
        self.df2 = pd.read_csv(f"{self.estacion}.csv", header=2) #separado por coma
        self.df["Y"] = self.df2["Year"]
        self.df["M"] = self.df2["Month"]
        self.df["D"] = self.df2["Day"]
        
        
        self.df['fecha'] = pd.to_datetime(self.df[["Y","M","D"]].rename(columns={'Y':'year',"M":"month","D":"day"}))
        self.df["J"] = self.df['fecha'].dt.dayofyear
        
        

        T = self.df2["Minute"]
        self.dt = T[2] - T[1]
        self.df["GHI"] = self.df2["GHI"]
        self.df["DNI"] = self.df2["DNI"]
        self.df["DHI"] = self.df2["DHI"]
        self.df["Temp"] = self.df2["Temperature"]
        self.df["E"] = list(map(self.getE, self.df.J, self.df.Y))
        self.df["H"] = self.df2['Hour']
        self.df["Minutes"] = self.df2["Minute"]
        self.df["Hdec"] = self.df2["Hour"] + (self.df2["Minute"])/60
        self.df['HS'] = self.getHS()
        
        self.df["Fn"] = list(map(self.Fn, self.df["J"], self.df["Y"]))
        self.df["Declinacion"] = list(map(self.delta, self.df.J, self.df.Y))
        self.df["DecGrados"] = self.df.Declinacion.apply(math.degrees)
        
        self.df['w'] = 15 * (12 - self.df['HS'])    
        self.df['wRad'] = self.df['w'].apply(math.radians)

        
        self.df["Ws"] = list(map(self.getWs, self.df["Declinacion"]))    
        self.df["CSZ"] = list(map(self.getCTZ, self.df["Declinacion"], self.df["wRad"]))
        self.df["CSTita"] = list(map(self.getCTita, self.df["Declinacion"], self.df["wRad"]))
        self.df["Icoll"] = list(map(self.getIcoll, self.df["CSTita"],self.df["GHI"], self.df["DNI"],self.df["DHI"] ))
        self.df["Tpanel"] = list(map(self.gettpanel, self.df["Icoll"],self.df["Temp"]))
       # self.diario = list(map(self.diario, np.array(self.df["Tpanel"].apply(np.array))))
        self.diaria3 = self.diaria(np.array(self.df["Tpanel"]),self.cant)
        self.clausura = self.clausuratotal(self.df["GHI"], self.df["DNI"], self.df["DHI"],self.df["Temp"],self.df["CSTita"])
        self.dias, self.diaria, total, k = self.diaria3
        self.Iestacional = self.estacional(self.dias, self.diaria,self.cant)

    #Calculo de CSZ
    
    def getHS(self):
        A = 1
        if self.gmt<=0:
            A = -1
        return self.df['Hdec'] + (4 * ((A * 15 * self.gmt)- (A*self.lon))+ self.df['E'])/60
    
            
    def getKtrp(self):
        if(self.alt>1000):
            return 0.7 + 1.6391 * 10**-3 * self.alt ** 0.5500 
        else:
            return 0.7570 + 1.0112 * 10**-5 * self.alt ** 1.1067
    
    def getE(self,n,y):
        
        
       
        if(y%4) == 0:
            gamma = 2 * np.pi * (n-1)/366
        else:
            gamma = 2 * np.pi * (n-1)/365
        
        
        cosg = np.cos(gamma) 
        sing = np.sin(gamma)
        
        cos2g = np.cos(2*gamma)
        sin2g = np.sin(2*gamma)
        
        
        E = 229.18 * (0.000075+ 0.001868 * cosg - 0.032077*sing - 0.014615*cos2g - 0.04089* sin2g )
        return E
    
    def Fn(self, n, y):
            
        if(y%4) == 0:
            gamma = 2 * math.pi * (n-1)/366
        else:
            gamma = 2 * math.pi * (n-1)/365
            
        fn = 1.000110 + 0.034221*math.cos(gamma) + 0.001280*math.sin(gamma) + 0.000719*math.cos(2*gamma) + 0.000077*math.sin(2*gamma)
        return fn
    
    #Hora solar 
    #devuelve valor hora solar
    #recorre cada fila del df
    def getHS(self):
        A = 1
        if self.gmt<=0:
            A = -1
        return self.df['Hdec'] + (4 * ((A * 15 * self.gmt)- (A*self.lon))+ self.df['E'])/60
    
  
    
    #Declinación solar
    #dado dia ordinal
    #devuelve declinacion en radianes
    def delta(self, n, y):
        
        if(y%4) == 0:
            gamma = 2 * math.pi * (n-1)/366
        else:
            gamma = 2 * math.pi * (n-1)/365
            
        delta = 0.006918 - 0.399912 * math.cos(gamma) + 0.070257 * math.sin(gamma) - 0.006758 * math.cos(2*gamma) + 0.000907*math.sin(2*gamma) - 0.002697*math.cos(3*gamma)+ 0.00148*math.sin(3*gamma)
        return delta
            
    def getWs(self, delta):
        return math.acos(-math.tan(delta) * (math.tan(math.radians(self.lat))))
    
    
    
    
    
    
    
    
    #CTZ
    #dado decliancio, lat y angulo horario
    #devuelve cos tita z en radianes
    def getCTZ(self, delta, omega):
        latR = math.radians(self.lat)    
        return (math.cos(latR) * math.cos(delta)* math.cos(omega)) + (math.sin(latR)*math.sin(delta))
        #return math.sin(delta) * math.sin(math.radians(lat)) + math.cos(delta) * math.cos(math.radians(lat))* math.cos(omega)

    def getCTita(self, delta, omega):
        latR = np.radians(self.lat)
        inclinacion = np.radians(self.beta)
        orien = np.radians(self.orientacion)
        return np.sin(delta) * np.sin(latR)* np.cos(inclinacion) - np.sin(delta) * np.cos(latR) * np.sin(inclinacion) * np.cos(orien) + np.cos(delta) * np.cos(latR) * np.cos(inclinacion) * np.cos(omega) + np.cos(delta) * np.sin(latR) * np.sin(inclinacion) * np.cos(orien) * np.cos(omega) + np.cos(delta) * np.sin(inclinacion) * np.sin(orien) * np.sin(omega)                       

    def getIcoll(self, CTita, ghi, dni, dhi):
        #ghi = self.df2["GHI"]
        #dni = self.df["DNI"]
        #dhi = self.df["DHI"]
        albedo = self.albedo
        inclinacion = np.radians(self.beta)
        #if CTita > 0 and CTZ > 0:
        Icoll = np.where(CTita >0, CTita * dni + dhi * ((1 + np.cos(inclinacion))/2) + albedo * ghi * ((1 - np.cos(inclinacion))/2), 0)
               
        return Icoll
        
    def gettpanel(self, Icoll, temperatura):
        #temperatura = self.df["Temp"]
        ktemp = self.ktemp
        PST = self.PST
        cvptm = self.cvptm
        area = self.area
        dt = self.dt
        tpanel = Icoll * ktemp + temperatura
        pele = PST * (Icoll/1000) * (1 + (cvptm/100) * (tpanel -25))
        if dt != 0:
            Energia = pele * (dt/60) 
        else: 
            Energia = pele
        
        return Energia/area
    
    def diaria(self, Tpanel,cant):
        dias = []
        cant = cant
       # s = Tpanel
        dt = self.dt
    #    h = len(s)
        if dt != 0:
            desp = 60/dt
        else:
            desp = 1

        vec = Tpanel.reshape(-1, 24)
        vec_suma = vec.sum(axis = 1)

        dias = np.linspace(1,365,365)

       

        total = cant*vec_suma.sum()/1000
     
        return dias, cant*vec_suma, total, plt.plot(dias, cant*vec_suma)
        
    def grafica(self, horario, dias):
        return plt.plot(dias, horario)
    
    def clausuratotal(self,ghi,dni,dhi,temp,csz):
        GHI=ghi
        DNI=dni
        DHI=dhi
        #cant= cant
        CSZ = csz
        temp = temp
        GHItotal = GHI.sum()
        DNItotal = DNI.sum()
        DHItotal = DHI.sum()
        Temprom = temp.mean()
        DNIhorizon = DNI * CSZ
        return GHItotal, DNItotal, DNIhorizon, DHItotal, Temprom
    
    def estacional(self,dias,total,cant):
        total = np.array(total)
        cant = cant
        dias = np.array(dias)
        df = pd.DataFrame({"dias":dias,"ghi":total})
        #df.index.name = "dias"
        #df = df.rename(columns= {0:"GHI"})
        #total = df.total
        dias = df.dias
        mask1 = ((dias < 172) & (80 <= dias))
        mask2 = ((172 <= dias) & (dias < 263))
        mask3 = ((263 <= dias) & (dias < 365))
        mask4 = ((dias < 80))
        mask5 = (355 <= dias)
        #df = df.set_index("dias")
        ghi = df.ghi
        otono = ghi[mask1]
        invierno = ghi[mask2]
        primavera = ghi[mask3]
        verano = ghi[mask4 | mask5]
        return cant*otono.sum(), cant*invierno.sum(), cant*primavera.sum(), cant*verano.sum()
       
        
    
        