# -*- coding: utf-8 -*-
"""
Created on Tue Jun 18 11:46:41 2024

@author: Windows
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize as opt
from tabulate import tabulate
from pathlib import Path #para organizar arquivos
import math

plt.style.use('seaborn-v0_8-dark-palette')

#aqui eu coloco todas as fontEs dos gráficos em Times new roman
plt.rcParams["font.family"] = "Times New Roman"

# plt.style.use('default')
file=Path('..')/'30C_CO2_LTA.xlsx'


#foram 3 medidas por ex
medidas=5
primeira_temperatura=30
#delta de temperatura em cada medida
passo=10

#definindo uma lista de erros para cada método 
s_freundlich=[]
s_lang=[]
s_sips=[]
s_dsl=[]
s_toth=[]

#definindo a coluna para tabela
column_sips=[]
column_dsl=[]
column_lang=[]
column_toth=[]
column_freund=[]

models=['Langmuir','Freundlich','DSL','Toth','Sips']

for model in range(len(models)):

#%%organizando os dados
    for i in range(medidas):
        if i==0:
            CF=file
            j=primeira_temperatura+i*passo
        else:
            j=primeira_temperatura+i*passo
            j_ant= primeira_temperatura+(i-1)*passo
            CF = str(CF).replace(str(j_ant), str(j))
            CF = Path(CF) 
        
        #cada medida está separada na planilha 1, 2 ,3 de cada carquivo
        sheets=['1','2','3']
    
        #aqui ele cria um dict contendo 3 dfs , as planilhas do aquivo
        df = {sheet: pd.read_excel(CF, sheet_name=sheet) for sheet in sheets}
        
        P = [df[sheet]['Pressão'] for sheet in sheets]
        conc = [df[sheet]['qab(mmol/g)'] for sheet in sheets]
    
        P_mean = pd.concat(P, axis=1).mean(axis=1)
        P_std = pd.concat(P, axis=1).std(axis=1)
    
        conc_mean = pd.concat(conc, axis=1).mean(axis=1)
        conc_std = pd.concat(conc, axis=1).std(axis=1)
    
        # tirando o zero de cadda um
        P_mean1 =P_mean[1:]
        conc_mean1 =conc_mean[1:]
        conc_std1 = conc_std [1:]
        P_std1 = P_std[1:]
        
        # para scipy optimize é importante transformar os objetos de pandas em listas
        P_mean_list= P_mean1.tolist()
        conc_mean_list= conc_mean1.tolist()
        conc_std_list=conc_std1.tolist()
        
        #%% Praparando parâmetros dos gráficos
        
        #aqui vamos organizar a cor para poder escolher sempre a mesma cor para o modelo e dados experimentais
        colour=['b','g','r','c','m','y','k']
        
        markerlist=['.', '*', 'o', 'v', 's', '^','p', '<', '>', '1', '2', '3', '4',  'h', 'H', '+', 'x', 'D', 'd', '|', '_']
        
        custom_linestyles = [
        (0, (1, 1)),      # dotted
        (5, (10, 3)),     # long dash with offset
        (0, (5, 10)),     # loosely dashed
        (0, (5, 5)),      # dashed
        (0, (5, 1)),      # densely dashed
        (0, (3, 10, 1, 10)),  # loosely dashdotted
        (0, (3, 5, 1, 5)),    # dashdotted
        (0, (3, 1, 1, 1)),    # densely dashdotted
        (0, (3, 5, 1, 5, 1, 5)),  # dashdotdotted
        (0, (3, 10, 1, 10, 1, 10)),  # loosely dashdotdotted
        (0, (3, 1, 1, 1, 1, 1))   # densely dashdotdotted
    ]
        
        legenda='Dados experimentais em '+str(j)+' C°'
        
        
        #vamos usar o P_geral que é apenas um grupo de dados contínuos
        P_geral = np.linspace(0, P_mean_list[-1],1000)
        
        #%%finalização dos parâmetros do gráfico
        plt.xlabel('Pressão(Bar)',fontsize=18)
        plt.ylabel('Concentração (mmol/g)',fontsize=18)
        plt.xlim(-0.5,52)
        plt.ylim(0, 4.5)
        plt.legend(loc='lower right') # localização da legenda
        plt.subplots_adjust(top=0.98) #aqui vamos colocar a borda superior perto do fim da imagem
        plt.grid(linestyle='--')

        #este aqui coloca pequenos ticks antes dos ticks principais(eg 10 20 30) de 5 em 5
        plt.gca().xaxis.set_minor_locator(plt.MultipleLocator(5)) 
    
        #%% Daqui em diante estou apenas adicionando a parte de análise em cima de cada curva
        
        # para scipy optimize é importante transformar os objetos de pandas em listas
        P_mean_list= P_mean1.tolist()
        conc_mean_list= conc_mean1.tolist()
        conc_std_list =conc_std1.tolist()
        '''
        calculando o erro do z-score , onde para 95% de confiança se * 1,96 
        se=std/(sqr(n)) - nocaso nossos experimentos foram feitos em triplicatas
        '''
        conc_error= (conc_std1/(3**(1/2)))*1.96
        #%% DSl - dual site langmuir
        
        if models[model] =='DSL':
        
            def dsl(p,d_b1,d_b2,qs_1,qs_2):
                q = (qs_1*d_b1*p)/(1+d_b1*p)+(qs_2*d_b2*p)/(1+d_b2*p)
                return q
            
            
            def fobj_dsl(param,P_mean_list,conc_mean_list):
                (d_b1,d_b2, qs_1,qs_2) = param
                # print('testeando',l_b,s_n,qs_lang)
                NEXP=len(P_mean_list)
                s=0
                for a in range(NEXP):
                    ycalc_i=dsl(P_mean_list[a],d_b1,d_b2, qs_1,qs_2) # aqui tem mudar em !!!
                    desvio_i=conc_mean_list[a] -ycalc_i
                    if math.isnan(conc_std_list[a]): #testando se é NaN
                        distancia_i=desvio_i**2
                        s+=distancia_i
                    else:
                        distancia_i=(desvio_i/conc_std_list[a])**2
                        s+=distancia_i
                # print(s)    
                return s
        
            resultado_dsl = opt.minimize(fobj_dsl,x0=[0,1,4,3],args=(P_mean_list,conc_mean_list),method='Powell')
            #aqui o methodo utilizado MUDA muito o resultado não sei se seria um viés de descoberta se eu usar só o méthodo que gere melhor resultado
            b_prov = resultado_dsl['x']
            (d_b1,d_b2, qs_1,qs_2) = b_prov
            s_dslII =fobj_dsl(b_prov,P_mean_list,conc_mean_list)
            
            #criando a tabela em pandas
            head_dsl= ['Temp','s^2','b1','b2','qs_1','qs_2']
            column_dsl.append([str(j)+' C°',s_dslII,d_b1, d_b2, qs_1, qs_2])
            df_dsl = pd.DataFrame(column_dsl, columns=head_dsl)
            df_dsl = df_dsl.set_index(['Temp'])
            
        
            #vamos gravar o 's' desse erro junto com os outros onde s_freundlich[i]  será uma lista de i objetos
            s_dsl.append(s_dslII)
            
            conc_dsl=[]
            for k in P_geral:
                conc_dsl.append(dsl(k,d_b1, d_b2, qs_1, qs_2))
                
            legenda_modelo_sips='DSL para ' +str(j)+' C°'  
            plt.figure(1)
            plt.plot(P_geral,conc_dsl,marker='',label=legenda_modelo_sips,color=colour[i],linestyle=custom_linestyles[i])
            plt.errorbar(P_mean1, conc_mean1,marker=markerlist[i],yerr=conc_error , xerr=None ,fmt='.',capsize= 3,label=legenda,color=colour[i],markersize=5)
            
        
        #%% Langmuir   
        elif models[model] =='Langmuir':
            def lang(p,l_b,qs_lang):
                return qs_lang*(l_b*p)/(1+l_b*p)
            
            
            def fobj_lang(param,P_mean_list,conc_mean_list):
                (l_b , qs_lang) = param
                NEXP=len(P_mean_list)
                s=0
                for a in range(NEXP):
                    ycalc_i=lang(P_mean_list[a],l_b,qs_lang) # aqui tem mudar em !!!
                    desvio_i=conc_mean_list[a] -ycalc_i
                    if math.isnan(conc_std_list[a]): #testando se é NaN
                        distancia_i=desvio_i**2
                        s+=distancia_i
                    else:
                        distancia_i=(desvio_i/conc_std_list[a])**2
                        s+=distancia_i
                # print(s)    
                return s

            resultado_lang = opt.minimize(fobj_lang,x0=[1,5],args=(P_mean_list,conc_mean_list),method='Powell')
            b_prov = resultado_lang['x']
            (l_b,qs_lang) = b_prov
            
            
            s_langII = fobj_lang(b_prov,P_mean_list,conc_mean_list)
            #Redefinindo a função com a correção
            def lang(p,l_b,qs_lang):
                return qs_lang*(l_b*p)/(1+l_b*p)
            
            #vamos gravar o 's' desse erro junto com os outros onde s_freundlich[i]  será uma lista de i objetos
            s_lang.append(fobj_lang(b_prov,P_mean_list,conc_mean_list))
            
            head_lang= ['Temp','s^2','b','qs']
            column_lang.append([str(j)+' C°',s_langII,l_b,qs_lang])
            df_lang = pd.DataFrame(column_lang, columns=head_lang)
            df_lang = df_lang.set_index(['Temp'])
            
            conc_lang=[]
            for k in P_geral:
                conc_lang.append(lang(k,l_b,qs_lang))
                
            legenda_modelo_lang='Langmuir para ' +str(j)+' C°'  
            plt.figure(2)
            plt.plot(P_geral,conc_lang,marker='',label=legenda_modelo_lang,color=colour[i],linestyle=custom_linestyles[i])
            plt.errorbar(P_mean1, conc_mean1,marker=markerlist[i],yerr=conc_error , xerr=None ,fmt='.',capsize= 3,label=legenda,color=colour[i],markersize=5)
            

        #%%freundlich
        elif models[model] =='Freundlich':
            def freund(p,b,n):
                return b*p**(1/n)

            def fobj_freund(parametros,P_mean_list,conc_mean_list):
                (b, n) = parametros
                NEXP=len(P_mean_list)
                s=0
                for a in range(NEXP):
                    ycalc_i=freund(P_mean_list[a],b,n) # aqui tem mudar em !!!
                    desvio_i=conc_mean_list[a] -ycalc_i
                    if math.isnan(conc_std_list[a]): #testando se é NaN
                        distancia_i=desvio_i**2
                        s+=distancia_i
                    else:
                        distancia_i=(desvio_i/conc_std_list[a])**2
                        s+=distancia_i
                # print(s)    
                return s

            resultado_freundlich = opt.minimize(fobj_freund,x0=[1,5],args=(P_mean_list,conc_mean_list),method='Powell')
            (b, n) = resultado_freundlich['x']
            resultado=(b,n)
            s_freund =fobj_freund(resultado,P_mean_list,conc_mean_list)

            #Redefinindo a função com a correção
            def freund(p,b,n):
                return b*p**(1/n)
            
            #vamos gravar o 's' desse erro junto com os outros onde s_freundlich[i]  será uma lista de i objetos
            s_freundlich.append(fobj_freund((b,n),P_mean_list,conc_mean_list))
            
            conc_freund=[]
            for k in P_geral:
                conc_freund.append(freund(k,b,n))
                
            legenda_modelo_freund='Freundlich para ' +str(j)+' C°'  
            #gráfico do fit de modelo
            plt.figure(3)
            plt.plot(P_geral,conc_freund,marker='',label=legenda_modelo_freund,color=colour[i],linestyle=custom_linestyles[i])
            #dados
            plt.errorbar(P_mean1, conc_mean1,marker=markerlist[i],yerr=conc_error , xerr=None ,fmt='.',capsize= 3,label=legenda,color=colour[i],markersize=5)
            
            
            #criando a tabela em pandas
            head_freund= ['Temp','s^2','b','n']
            column_freund.append([str(j)+' C°',s_freund,b,n])
            df_freund = pd.DataFrame(column_freund, columns=head_freund)
            df_freund = df_freund.set_index(['Temp'])
        #%%toth
        elif models[model] =='Toth':
            def toth(p,t_b,t_n,qs_toth):
                q = qs_toth*(t_b*p)/((1+(t_b*p)**(t_n))**(1/t_n))
                return q
            

            
            def fobj_toth(param,P_mean_list,conc_mean_list):
                (t_b ,t_n, qs_toth) = param
                # print('testeando',l_b,s_n,qs_lang)
                NEXP=len(P_mean_list)
                s=0
                for a in range(NEXP):
                    ycalc_i=toth(P_mean_list[a],t_b,t_n,qs_toth) # aqui tem mudar em !!!
                    desvio_i=conc_mean_list[a] -ycalc_i
                    if math.isnan(conc_std_list[a]): #testando se é NaN
                        distancia_i=desvio_i**2
                        s+=distancia_i
                    else:
                        distancia_i=(desvio_i/conc_std_list[a])**2
                        s+=distancia_i
                # print(s)    
                return s

            resultado_toth = opt.minimize(fobj_toth,x0=[1000,0.2,5],args=(P_mean_list,conc_mean_list))
            b_prov = resultado_toth['x']
            (t_b,t_n,qs_toth) = b_prov
            s_tothII =fobj_toth(b_prov,P_mean_list,conc_mean_list)
            
            head_toth= ['Temp','s^2','b','n','qs']
            column_toth.append([str(j)+' C°',s_tothII,t_b,t_n,qs_toth])
            df_toth = pd.DataFrame(column_toth, columns=head_toth)
            df_toth = df_toth.set_index(['Temp'])
            
            
            #vamos gravar o 's' desse erro junto com os outros onde s_freundlich[i]  será uma lista de i objetos
            s_toth.append(s_tothII)
            
            conc_sips=[]
            for k in P_geral:
                conc_sips.append(toth(k,t_b,t_n,qs_toth))
                
            legenda_modelo_lang='toth para ' +str(j)+' C°'  
            plt.figure(4)
            plt.plot(P_geral,conc_sips,marker='',label=legenda_modelo_lang,color=colour[i],linestyle=custom_linestyles[i])
            plt.errorbar(P_mean1, conc_mean1,yerr=conc_error ,xerr=None ,fmt='.',capsize= 3,label=legenda,color=colour[i],marker=markerlist[i])
            

            
            
        #%%Sips
        elif models[model] =='Sips':
            
            def sips(p,s_b,s_n,qs_sips):
                q = qs_sips*((s_b*p)**(1/s_n))/(1+s_b*p**(1/s_n))
                return q
            
            
            def fobj_sips(param,P_mean_list,conc_mean_list):
                (s_b ,s_n, qs_sips) = param
                # print('testeando',l_b,s_n,qs_lang)
                NEXP=len(P_mean_list)
                s=0
                for a in range(NEXP):
                    ycalc_i=sips(P_mean_list[a],s_b,s_n,qs_sips) # aqui tem mudar em !!!
                    desvio_i=conc_mean_list[a] -ycalc_i
                    if math.isnan(conc_std_list[a]): #testando se é NaN
                        distancia_i=desvio_i**2
                        s+=distancia_i
                    else:
                        distancia_i=(desvio_i/conc_std_list[a])**2
                        s+=distancia_i
                # print(s)    
                return s

            resultado_sips = opt.minimize(fobj_sips,x0=[0.5,8,5],args=(P_mean_list,conc_mean_list),method='Nelder-Mead')
            b_prov = resultado_sips['x']
            (s_b,s_n,qs_sips) = b_prov
            s_sipsII =fobj_sips(b_prov,P_mean_list,conc_mean_list)
            
            #criando a tabela em pandas
            head_sips= ['Temp','s^2','b','n','qs']
            column_sips.append([str(j)+' C°',s_sipsII,s_b,s_n,qs_sips])
            df_sips = pd.DataFrame(column_sips, columns=head_sips)
            df_sips = df_sips.set_index(['Temp'])
            

            #vamos gravar o 's' desse erro junto com os outros onde s_freundlich[i]  será uma lista de i objetos
            s_sips.append(s_sipsII)
            
            
            conc_sips=[]
            for k in P_geral:
                conc_sips.append(sips(k,s_b,s_n,qs_sips))
                
            legenda_modelo_sips='SIPS para ' +str(j)+' C°'  
            plt.figure(5)
            plt.plot(P_geral,conc_sips,marker='',label=legenda_modelo_sips,color=colour[i],linestyle=custom_linestyles[i])
            plt.errorbar(P_mean1, conc_mean1,marker=markerlist[i],yerr=conc_error , xerr=None ,fmt='.',capsize= 3,label=legenda,color=colour[i],markersize=5)
            
        
        #%%erro

        
        
    
    
    


#%%Preparando a tabela para salvar

# vamos arredondar os dados da tabela
'''
df_dsl= df_dsl.round(4)



df_dsl.loc['Total'] = pd.Series(df_dsl['s^2'].sum(), index=['s^2'])


print(df_dsl)
df_dsl.to_csv('dados_dsl.csv')
'''
