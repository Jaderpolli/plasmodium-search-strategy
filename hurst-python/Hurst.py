import os
import csv
import numpy as np
from sklearn.linear_model import LinearRegression
import math 
import matplotlib.pyplot as plt
from scipy import stats

#THIS PROGRAM CALCULATES THE HURST COEFICIENT OF A TEMPORAL SEREIES OF ORGANISM MOVEMENT.

tau_max=25            #<--- maximum time_window range.
tau=0                 #<--- the variable that contains the time_window range in analyses (will be iterated to tau_max)
treshhold=0           #<--- velocity treshhold. (if not necessary put 0)
resolucao_temporal=4  #<--- temporal resolution
msd_medio=0.1         #<--- average mean square displacement (will be used on the calculation, initial value doesn't matter)
janela_temporal=[]    #<--- array that stores the time data for plotting
serie_msd_bicho=[]    #<--- array that temporarily stores the MSD of a given organism
matriz_resultados=[]  #<--- matrix that stores the Hurst coefficient of every organism
resultado=[]          #<--- stores the Hurst data into a table, that we will save into a csv file


#first, we open RESULTADOS.csv: a file to organize the displacement data time series for every organism
with open('RESULTS.csv', 'w', encoding='UTF8', newline='') as result:
    writer = csv.writer(result)
    header=[]
    header.append("tempo (s):") #<--- here we create the header (the time labels)
    for i in range(tau_max):
        header.append((i+1)*resolucao_temporal) #<--- and here we populate it with the time instants of appropriate time resolution.
    writer.writerow(header)
    #The movement data of every organism are organized into separeted files.
    #the "For" loop bellow will sweep across all those files.
    for filename in os.listdir("raw_data"):                      #<--- for every file...
        with open(os.path.join("raw_data", filename), 'r') as f: #<--- ...it opens it and do the analysis
            serie_msd_bicho=[]
            janela_temporal=[]
            csvreader = csv.reader(f, delimiter=",")
            header = []
            header = next(csvreader)
            contador=0           #<--- this is a general variable to count things
            hurst_medio=[]       #<--- this stores average hurst coefficient for the organism paths, at every time window
            dados =[]            #<--- this array will stores the organism data from the file...
            for row in csvreader:#<--- ...here its storing it
                dados.append(row)
            n_passos=len(dados)
            for tau in range(tau_max): #<--- ...now, for every time window, it proceds with the hurst calculations
                soma_msd=0
                contador=0
                if n_passos>=(4*(tau+1)): #<--- ...checks if the organism path is suficiently big to generate hurst data at the current time window under analysis.
                                          #The organism trajectory must be at least 4 times longer than the current time window, otherwise the statistics wont be too reliable.
                    for i in range(len(dados)-(tau+1)):  #Here it computes the mean square displacement.
                        msd=pow(float(dados[i+(tau+1)][0])-float(dados[i][0]),2)+pow(float(dados[i+(tau+1)][1])-float(dados[i][1]),2)
                        velocidade=pow(msd,0.5)/((tau+1)*resolucao_temporal)
                        if velocidade>treshhold:         #If it is faster than the treshold, it stores de MSD...
                            soma_msd=soma_msd+msd
                            contador=contador+1
                    msd_medio=soma_msd/contador          #...and here we take the MSD average
                    serie_msd_bicho.append(msd_medio)    #<--- storing the average MSD into a list...
                    janela_temporal.append((tau+1))      #...and the corresponding time intervals into another list
            
            #We now proceed to the organism hurst calculations. First, we find the curves in log scale:
            for i in range(len(serie_msd_bicho)):
                serie_msd_bicho[i]=math.log10(serie_msd_bicho[i])
                janela_temporal[i]=math.log10(janela_temporal[i])
            hurst=[]
            t=[]

            for i in range(len(serie_msd_bicho)-1):
                hurst_dummy=0.5*(serie_msd_bicho[i+1]-serie_msd_bicho[i])/(janela_temporal[i+1]-janela_temporal[i]) #<--- here we calculate the Hurst coefficient of a path
                if hurst_dummy>0:
                    hurst.append(hurst_dummy)
                if hurst_dummy<=0:
                    hurst.append("---")
            resultado=[]
            resultado.append(f.name)       #<--- on the results table, we first store the file (and hence the organism) name...
            for l in range(len(hurst)):
                resultado.append(hurst[l]) #<--- and then the respective hurst coefficient calculated from the average MSDs
            writer.writerow(resultado)
            matriz_resultados.append(hurst)#<--- The hurst of this single organism will then be saved on this matrix (more precisely, in one of its lines).
    hurst_medio=[0]*(tau_max-1)            #<--- This variable will be used to average the Hurst coefficients of all the organisms.
    contador=0
    t=[]
    for j in range(len(hurst_medio)):
        contador=0
        for i in range(len(matriz_resultados)):
            if j<=(len(matriz_resultados[i])-1):#<--- since not every organism has paths the same size, not every organism contributes to every time window averages. Here we check it.
                if type(matriz_resultados[i][j])== float:
                    hurst_medio[j]=hurst_medio[j]+matriz_resultados[i][j]
                    contador=contador+1
        hurst_medio[j]=hurst_medio[j]/contador #<--- averaging the Hurst coefficients.
        t.append((j+1)*4)
    print(hurst_medio)
    resultado=[]
    resultado.append("MÉDIO")  #<--- now we'll store the average into the results file aswell.
    for l in range(len(hurst_medio)):
        resultado.append(hurst_medio[l])
    writer.writerow(resultado)

plt.plot(t,hurst_medio,lw=0.2, alpha=0.8, label="treshhold=0\u03BCm/s", marker = 'o') #and finally, plot the graph.
plt.xlabel('time window')
plt.ylabel('Hurst coefficient')
plt.ylim(0.2,0.8)
plt.legend(loc='upper center')
plt.show()