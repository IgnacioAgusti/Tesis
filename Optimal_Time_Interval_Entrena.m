% Optimal_Time_Interval_Entrena.m
% ================================
% DESCRIPCIÓN: Programa permite leer la primera parte de la estructura
% STADISTIC_Results y visualizar resultados
% SALIDA: INDICES_ENTRENA
%==========================================================================
% PARAMETROS DE INICIACION
% -----------------------------------------
sigma_so = 1;            % UMBRAL DE LA DESVIACION
Por_Ef_Est = 20;         % PORCENTAJE DE EFICIENCIA ESTIMADO (EMPIRICO)

% DATA GEOMAGNETICA DE LAS ESTACIONES A USAR
% ----------------------------------------------
%      Data_Stations = 3; % (1) Japan (2) Europa (3) otras ....
[Data,Label_Station] = Data_Geomagnetica(Data_Stations);
% ----------------------------------------------------
%  COORDENADAS Y ETIQUETAS DE LOS EVENTOS SISMICOS
% ----------------------------------------------------
switch Data_Stations
 case 1
  load ('Entrena_Etiqueta_Eventos_Japan_II_6M_2000_2019.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont');
  disp(' Entrenamiento JAPON')
 case 2
  load ('Etiquetas_DOUEntrenamiento2_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont');
  disp(' Entrenamiento Europa silencio, sismos simulados')
  
 case 3
  load ('Etiquetas_USA_Entrenamiento_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont');
  disp('Entrenamiento USA')
 case 4
  load ('Etiquetas_CHILE_Entrenamiento_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont')  
  disp(' Entrenamiento CHILE')
 case 5
  load ('Etiquetas_EUROPA_DUREntrenamiento_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont')  
  disp(' Entrenamiento DUR')
  
 case 6
  load ('Etiquetas_EUROPA_PEGEntrenamiento_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont')
  disp(' Entrenamiento PEG')
 case 7
  load ('Etiquetas_PETEntrenamiento_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont')
  disp(' Entrenamiento PET')
 case 8
  load ('Etiquetas_KAK_Entrenamiento_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont')
  disp(' Entrenamiento KAK')
  
end

Lxx = Lx;
Lyy = Ly;

%%
% -------------------------------------------------------------------------
% CAMPOS
% -------
INDICES_ACIERTOS_ALL = [];
INDICES_FALLOS_ALL = [];
Indices_Aciertos_Produce = {};
Indices_Entrena = {};
PA = [];
Por_Ef = [];
Mag_max = 10;
% -------------------------------------------------------------------------
%

% --------------------------------------------------------
% VALOR ACTUAL DE GAMMA
% --------------------------------------------------------
S_Ent = Stadistic_Entrena;
s_ent = S_Ent{4,1}; 
Station_Used_Now =  s_ent.Station_Used_Now;
Actual_Gamma = s_ent.Gamma_2;
% -------------------------------------------------------------------------
% OBTENCION DEL PORCENTAJE DE EFICIENICIA PROMEDIO (MINIMA) 
% -------------------------------------------------------------------------
Total_Stations_Combinations = 4; 
for Actual_Station = 1:Total_Stations_Combinations;
 s_ent = S_Ent{Actual_Station,1};
 Aciertos = s_ent.Aciertos;
 iMc = s_ent.iMc;
 % --------------------------------------------------------
 % INDICES CON EFICIENCIA PORCENTUAL >= Por_Ef
 % -----------------------------------------------------------------------
 % MAXIMO DE PROBABILIDAD EN LOS ACIERTOS
end
% -------------------------------------------------------------------------
POR_EF = min(Por_Ef); % PORCENTAJE DE EFICIENCIA MINIMA
% -------------------------------------------------------------------------
% LAZO PARA VARIOS VALORES DE GAMMA
% ------------------------------------

for j_gamma = 1:length(Actual_Gamma);
 
 % -----------------------------------
 if band == 1
  disp('----------------------------------------------------------')
  disp(['Indice de gamma: ' num2str([j_gamma]) ' Valor da gamma ' ...
   num2str(Actual_Gamma(j_gamma))])
  disp('----------------------------------------------------------')
 end
 % 
 Por_Ef = [];
 for Actual_Station = 1:Total_Stations_Combinations;
  s_ent = S_Ent{Actual_Station,j_gamma};
  Aciertos = s_ent.Aciertos;
  Fallos = s_ent.Fallos;
  iMc = s_ent.iMc;
  M_dh = s_ent.M_dh;
  Dh = s_ent.Dh;
  Num_h = s_ent.Num_h;
  Int_H = s_ent.Int_H;
  Total_a = s_ent.Total_a;
  Pa = s_ent.Pa;
  t = s_ent.t;
  T = t - t(1);
  % INDICES CON EFICIENCIA PORCENTUAL >= Por_Ef
  % -----------------------------------------------------------------------
  %CASO USA Y JAPON:=======================================================
  Por_Ef(Actual_Station) = Por_Ef_Est;
  % =======================================================================
  %PRUEBA Chile:
  %Por_Ef(Actual_Station) = round(Aciertos(iMc) - sigma_so*std(Aciertos)); 
  % =======================================================================
  ind_aciertos = find(Aciertos >= Por_Ef(Actual_Station));
  ind_fallos = find(Aciertos < Por_Ef(Actual_Station));
  
  % CREA LOS INDICES DEL ENTRENMIENTO POR ESTACION
  % -----------------------------------------------------------------------
  Indices_Aciertos_Produce{Actual_Station} = struct('ind_aciertos',...
   ind_aciertos,'ind_fallos',ind_fallos);
  INDICES_ACIERTOS_ALL = [INDICES_ACIERTOS_ALL;ind_aciertos];
  INDICES_FALLOS_ALL = [INDICES_FALLOS_ALL;ind_fallos];
  Pa(:,ind_aciertos) = Pa(:,ind_aciertos);
  Pa(:,ind_fallos) = 0.125 ; % SOLO ES EL BACKGROUND DE LA IMAGEN
  PA(:,:,Actual_Station) = Pa;
 end
 % -----------------------------
 Por_Ef =  min(Por_Ef);
 % -----------------------------
 % ELIMINA INDICES REPETIDOS (aciertos)
 % ------------------------------------------------------------------------
 INDICES_ACIERTOS_ALL_sort = sort(INDICES_ACIERTOS_ALL);
 IND_ACIERTOS = [];
 
 while length(INDICES_ACIERTOS_ALL_sort) > 0
  ind = find(INDICES_ACIERTOS_ALL_sort==INDICES_ACIERTOS_ALL_sort(1));
  IND_ACIERTOS = [IND_ACIERTOS;INDICES_ACIERTOS_ALL_sort(1)];
  INDICES_ACIERTOS_ALL_sort(ind) = [];
 end
 % ------------------------------------------------------------------------
 % ELIMINA INDICES REPETIDOS (fallos)
 % ------------------------------------------------------------------------
 INDICES_FALLOS_ALL_sort = sort(INDICES_FALLOS_ALL);
 IND_FALLOS = [];
 
 while length(INDICES_FALLOS_ALL_sort) > 0
  ind = find(INDICES_FALLOS_ALL_sort==INDICES_FALLOS_ALL_sort(1));
  IND_FALLOS = [IND_FALLOS;INDICES_FALLOS_ALL_sort(1)];
  INDICES_FALLOS_ALL_sort(ind) = [];
 end
 % ------------------------------------------------------------------------
 Indices_Entrena{j_gamma} = struct('IND_ACIERTOS',IND_ACIERTOS,...
  'IND_FALLOS',IND_FALLOS,'Indices_Aciertos_Produce',...
  Indices_Aciertos_Produce,'Por_Ef',Por_Ef,'Actual_Gamma',...
  Actual_Gamma(j_gamma));
 % ========================================================================
 INT_H = Int_H(IND_ACIERTOS);
 lx_ac = Int_H;
 ly_ac = Por_Ef*ones(1,length(Int_H));
 % ========================================================================
 % VISUALIZACION
 if band == 1
  close all
  set(figure(1),'Position',[ 50 490 1245 206],'Color','W')
  plot(Int_H,Aciertos,'.-r',Int_H,Aciertos,'or',Int_H,Fallos,'.-b',...
  Int_H,Fallos,'ob',Int_H,Total_a,'.-m',Int_H,Total_a,'om','LineWidth',[2])
  hold on
  plot(Int_H(iMc),Aciertos(iMc),'+k','MarkerSize',[12],'LineWidth',[2])
  plot(Int_H(ind_aciertos),Aciertos(ind_aciertos),'og',...
   Int_H,Aciertos,'+k','MarkerSize',[14],'LineWidth',[2])
  xlabel('Time interval from 0 to 24 h')
  ylabel('Eficiency in (%)')
  title(['Efmax: [' num2str(Aciertos(iMc)) '(%)] --- Dhm: ['...
   num2str(round(Int_H(iMc))) '] hours' ' (Training phase) '])
  plot(lx_ac,ly_ac,'--g','LineWidth',[2])
  hold off
  xlim([Dh(1) Dh(end)/Num_h])
  ylim([0 100])
  g = get(1);
  g1 = g.Children;
  set(g1,'Position',[0.0415    0.1250    0.9514    0.8103],'Box','on')
  legend('Thrue','','False','','Total','Location','Best')
  ax1 = gca;                    % current axes
  ax1.XTick = [Dh(1):1:Dh(end)/Num_h];
  grid
  Ejes_Visibles(1)
  % =======================================================================
  set(figure(2),'Position',[51 44 248 365],'Color','W')
  for i=1:Total_Stations_Combinations
   Pm = PA(:,:,i)';
   subplot(4,2,i)
   % ------------
   % SURF
   % ------------
   [xa,ya] = meshgrid(linspace(T(1),T(end),size(Pm,2)),...
    linspace(1,M_dh(end),size(Pm,1)));
   surf(xa,ya/Num_h,Pm),shading interp,view([0 90])
   colorbar,colormap jet
   grid,axis xy,axis on
   Pm = mean(PA(:,:,i)');
   Pm = (Pm-max(Pm))/(max(Pm)-min(Pm))*Mag_max;
   hold on
   plot((Lx'-Lxx(1,1)),Ly'-Mag_max,'-b','LineWidth',[2])
   plot((Lx_Cont'-Lxx(1,1)),Ly_Cont'-Mag_max,'--g','LineWidth',[3])
   plot(T,Pm,'-r','LineWidth',[3]),grid
   plot(T,Pm,'-k','LineWidth',[1]),grid
   hold off
   if i < 4 title(cell2mat(['Station: ' Label_Station(i)]))
   else
    title('All Station')
   end
   ax1 = gca;                   % current axes
   ax1.YTick = [-10:2:M_dh(end)];
   ax1.GridLineStyle = '--';
   ax1.GridColor = [1 0 0];
   xlim([0 T(end)])
   ylim([-10 24])
  end
  g = get(figure(2));
  g1 = g.Children;
  set(g1(1),'Position',[0.9674    0.0807    0.0141    0.4025],'Box','on');
  set(g1(2),'Position',[0.5363    0.0806    0.4239    0.4050],'Box','on');
  if not(Total_Stations_Combinations == 1)
   set(g1(3),'Position',[0.4751    0.0773    0.0141    0.4095],'Box','on');
   set(g1(4),'Position',[0.0454    0.0773    0.4225    0.4095],'Box','on');
   set(g1(5),'Position',[0.9685    0.5594    0.0141    0.4060],'Box','on');
   set(g1(6),'Position',[0.5363    0.5596    0.4251    0.4058],'Box','on');
   set(g1(7),'Position',[0.4760    0.5563    0.0141    0.4029],'Box','on');
   set(g1(8),'Position',[0.0449    0.5552    0.4241    0.4040],'Box','on');
  end
  grid
  Ejes_Visibles(2)
  pause
 end
end
% -------------------------------------------------------------------------
fprintf('-------------------------------------------------------------------------\n')
disp('   La base datos de Indices de Entrenamiento fue creada para las estaciones:  ')
disp(Label_Station')
fprintf('-------------------------------------------------------------------------\n')
% -------------------------------------------------------------------------
% SALIDAS
% ----------------------------------------
save INDICES_ENTRENA Indices_Entrena
% ----------------------------------------