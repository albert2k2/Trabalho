# -*- coding: utf-8 -*-
"""
Created on Mon Sep 23 15:27:12 2024

@author: ALBERTO


# Modified Peng-Robinson equation of state for CO2/hydrocarbon systems within nanopores Gang Yang,Xiaoli Li
"""



import numpy as np
from scipy.optimize import fsolve


#variáveis zeradas apenas para começar e não dar bug
P,V,T=[0,0,0]




# vamos tentar uma PR hoje

w= 0.2249 #co2 pr (fator acêntrico)

temperatura_celsius= 40
T=temperatura_celsius+ 273.15
P= 50e5                          # 50 bar = 50*10**5  #Pa

mm_co2=44.009 #g/mol


R=8.314462            #J/K mol


Tc=304.1288           #K
Pc=7.3773*(10**6)     #Pa

Tt=216.593            #K
Pt=0.51795*(10**6)    #Pa

Tr=T/Tc

k=0.37464 + 1.54226*w - 0.26992*(w**2)

alpha=(1+k*(1-Tr**(1/2)))**2

a=0.45724*((R**2)*(Tc**2)/(Pc))*alpha

b=0.00778*R*Tc/Pc

def peng_robinson_volume(V,P,T):
    return -P+ (R*T)/(V-b)-(a/((V**2)+(2*b*V)-(b**2)))

def peng_robinson_pressao(P,V,T):
    return -P+ (R*T)/(V-b)-(a/((V**2)+(2*b*V)-(b**2)))

def peng_robinson_temperatura(T,P,V):
    return -P+ (R*T)/(V-b)-(a/((V**2)+(2*b*V)-(b**2)))

#var = P , V ou T
variavel='V'
# Chutes iniciais (usando PV=RT)
if variavel == 'V':
    V_chute=R*T/P
elif variavel == 'P':
    P_chute=R*T/V
elif variavel == 'T':
    T_chute=P*V/R
else:
    print('variável não detectada  Apenas V , P e T')

    
#aplicando a solução


if variavel == 'V':
    V_sol=fsolve(peng_robinson_volume,V_chute,args=(P,T) )
    print(f"O volume molar é {V_sol[0]} m³/mol")
    rho=(V_sol[0]**-1)*mm_co2/1000
    print(f"A massa específica é {rho} kg/m³")
    conc= 1/V_sol # que será mol/m³
    conc_L=conc/(1000)
    print(f"A concentraão no sistema puro será : {conc} mol/m³ ou {conc_L}mol/L")
elif variavel == 'P':
    P_sol = fsolve(peng_robinson_pressao, P_chute, args=(V, T))
    print(f"A pressão é {P_sol[0]} Pa")
elif variavel == 'T':
    T_sol = fsolve(peng_robinson_temperatura, T_chute, args=(P, V))
    print(f"A temperatura é {V_sol[0]} kelvin")
else:
    print('variável não detectada  Apenas V , P e T')