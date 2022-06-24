% plot_resultados_general.m
% ================================
% DESCRIPCIÓN: generaliza el proceso y visualiza las salidas dado un cierto
% porcentaje de filtraje para la calibración
% SALIDA: Por_Ef_Est Stadistic_Results T_entrena ...
% Data_Stations label_e  IRISIDS Dist_respectivas Nombre_caso Data_Stations
% =========================================================================
clear all
%% carpeta madre
cd 'C:\Tesis\Programas\Programas_Misael_22_Oct_2020'
%PARA 7, MUESTRA 0%
%% ENTRADAS
label_e ='c';
Data_Stations =8;
Por_Ef_Est=25;
%% PROCESO
label = {label_e};
Nombre_general = 'STADISTIC_Results_';
Nombre_caso = {'Japon','DOU','USA','Chile','DUR','PEG','PET','KAK'};
Directorio = dir([Nombre_general char(Nombre_caso(Data_Stations)) char(label) '.mat']); 
clc
display(['Nombre archivo cargado: ' Directorio .name ])
display(['Zona: ' char(Nombre_caso(Data_Stations)) '. Label: ' char(label) '.'])
%% ENTRENA, PRODUCE Y PRECURSOR
% ENTRENA =================================================================
band =0;
load(Directorio .name,'Stadistic_Results')
Stadistic_Entrena=Stadistic_Results(:,:,1);
Optimal_Time_Interval_Entrena
T_entrena = T;
% PRODUCE =================================================================
clearvars  -except Directorio  Por_Ef_Est Stadistic_Results T_entrena ...
 Data_Stations label_e Nombre_caso Data_Stations 
band = 0;
Stadistic_Produce=Stadistic_Results(:,:,2);
Optimal_Time_Interval_Produce
% PRECURSOR ===============================================================
clearvars  -except Directorio  Por_Ef_Est Stadistic_Results T_entrena ...
 Data_Stations label_e  IRISIDS Dist_respectivas Nombre_caso Data_Stations
Precursor_Activation_Gamma_Dependence_new
% =========================================================================