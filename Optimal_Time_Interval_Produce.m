% Optimal_Time_Interval_Produce.m
% ================================
% DESCRIPCIÓN: Programa permite leer la segunda parte de la estructura
% STADISTIC_Results y visualizar resultados
% SALIDA: Precursor_Activation_Parameters
% =========================================================================
close all
% CAMPOS
%--------
PA = [];
IA_MEAN_Y = [];
IA_MEAN_X = [];
Mag_max = 10;
ka = 2e-1;
n1 = 4;      % FILTRO EN X DE LA IMAGEN
n2 = 4;      % FILTRO EN Y DE LA IMAGEN
% ----------------------------------------------
% DATA GEOMAGNETICA DE LAS ESTACIONES A USAR
% ----------------------------------------------
%        Data_Stations = 3;          % (1) Japan (2) Europa (3) otras ....
[Data,Label_Station] = Data_Geomagnetica(Data_Stations);

switch Data_Stations
 case 1
  load('Produce_Etiqueta_Eventos_Japan_II_6M_2000_2019.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont','IRISIDS','Dist_respectivas')
  disp(' Produccion JAPON')
 case 2
  load('Etiquetas_DOUProduccion3_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont');
  disp('Produccion  Europa silencio, sismos simulados')
  
 case 3
  load('Etiquetas_USA_Produccion_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont','IRISIDS','T','Dist_respectivas');
  disp('Produccion USA')
 case 4
  load('Etiquetas_CHILE_Produccion_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont','IRISIDS','T','Dist_respectivas');
  disp(' Produccion CHILE')
 case 5
  load('Etiquetas_EUROPA_DURProduccion_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont','IRISID','T','Dist_respectivas');
  disp(' Produccion DUR')
  
 case 6
  load('Etiquetas_EUROPA_PEGProduccion_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont','IRISIDS','T','Dist_respectivas');
  disp(' Produccion PEG')
 case 7
   load('Etiquetas_PETProduccion2_Ready.mat','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont','IRISIDS');
  disp(' Produccion PET')
 case 8
  load('Etiquetas_KAK_Produccion_Ready','Lx_c','Ly','Lx'...
   ,'Ly_c','Mag_S','Lon_S','Lat_S','Depthkm_S','Lx_Cont','Ly_Cont','IRISIDS');
  disp(' Produccion KAK')
  
end
Lxx = Lx;
Lyy = Ly;
%%
% -------------------------------------------------------------------------

% CARGA LOS INDICES DEL ENTRENAMIENTO
% ---------------------------------------
load INDICES_ENTRENA
% OBTIENE EL VALOR DE Gamma_2 de la estacion 1 y gamma 1
% ----------------------------------------------------------
S_Prod = Stadistic_Produce;
s_prod = S_Prod{1,1};
Gamma_2 = s_prod.Gamma_2;
% ----------------------------------------------------------

% LAZO ARA GAMMA
% ---------------------------------------------------------------
for j_gamma = 1:length(Gamma_2);
 Ind_Ent = cell2mat(Indices_Entrena(j_gamma));
 IND_ACIERTOS = Ind_Ent.IND_ACIERTOS;
 IND_FALLOS = Ind_Ent.IND_FALLOS;
 Indices_Aciertos_Produce = Ind_Ent.Indices_Aciertos_Produce;
 Por_Ef = Ind_Ent.Por_Ef;
 Actual_Gamma = Ind_Ent.Actual_Gamma;
 % ------------------------------------------------------------------------
 if band == 1
  disp('----------------------------------------------------------')
  disp(['Valor de gamma: ' num2str(Actual_Gamma)] )
  disp('----------------------------------------------------------')
 end
 % ------------------------------------------------------------------------
 Total_Stations_Combinations = 4;      % Station_Used_Now = [1 2 3 4];
 S_Prod = Stadistic_Produce;
 
 for Actual_Station = 1:Total_Stations_Combinations;
  s_prod = S_Prod{Actual_Station,j_gamma};
  f1 = s_prod.f1;
  f2 = s_prod.f2;
  Md = s_prod.Md;
  M_dh = s_prod.M_dh;
  Dh = s_prod.Dh;
  Num_h = s_prod.Num_h;
  Int_H = s_prod.Int_H;
  Nv = s_prod.Nv;
  Pa = s_prod.Pa;
  Pt = s_prod.Pt;
  t = s_prod.t;
  T = t-t(1);
  Station_Used_Now = s_prod.Station_Used_Now;
  Ind_Ac = Indices_Aciertos_Produce;
  ind_aciertos = Ind_Ac.ind_aciertos;
  ind_fallos = Ind_Ac.ind_fallos;
  Pa(:,ind_aciertos) = Pa(:,ind_aciertos) + ka;
  Pa(:,ind_fallos) = 0;
  PA(:,:,Actual_Station) = Pa;
 end
 % ========================================================================
 if band == 1
  close all
  
  set(figure(1),'Position',[232    45   951   622],'Color','W')
  for i=1:Total_Stations_Combinations
   Pm = PA(:,:,i)';
   subplot(2,2,i)
   
   % SURF
   % ----
   [xa,ya] = meshgrid(linspace(T(1),T(end),size(Pm,2)),linspace(1,M_dh(end),size(Pm,1)));
   surf(xa,ya/Num_h,Pm),shading interp,view([0 90])
   colorbar,colormap jet,alpha .7
   grid,axis xy,axis on
   Pm = mean(PA(:,:,i)');
   Pm = (Pm-max(Pm))/(max(Pm)-min(Pm))*Mag_max;
   hold on
   plot((Lx'-Lxx(1,1)),Ly'-Mag_max,'-b','LineWidth',[2])
   plot((Lx_Cont'-Lxx(1,1)),Ly_Cont'-Mag_max,'--g','LineWidth',[3])
   plot(T,Pm,'-r','LineWidth',[3]),grid
   plot(T,Pm,'-k','LineWidth',[1]),grid
   hold off
   grid
   ax1 = gca;                   % current axes
   ax1.YTick = [-10:2:M_dh(end)];
   ax1.GridLineStyle = '--';
   ax1.GridColor = [0 0 0];
   ax1.GridAlpha = 1;
   ax1.LineWidth = 0.5;
   ax1.YMinorGrid = 'on';
   ax1.XMinorGrid = 'on';
   if i < 4 title(cell2mat(['Estación: ' Label_Station(i)]))
   else
    title('Todas las estaciones')
    %title('All Station')
   end
   xlim([0 T(end)])
   ylim([-10 24])
  end
  g = get(figure(1));
  g1 = g.Children;
  set(g1(1),'Position',[0.9674    0.0807    0.0141    0.4025],'Box','on');
  set(g1(2),'Position',[0.5363    0.0806    0.4239    0.4050],'Box','on');
  set(g1(3),'Position',[0.4751    0.0773    0.0141    0.4095],'Box','on');
  set(g1(4),'Position',[0.0454    0.0773    0.4225    0.4095],'Box','on');
  set(g1(5),'Position',[0.9685    0.5594    0.0141    0.4060],'Box','on');
  set(g1(6),'Position',[0.5363    0.5596    0.4251    0.4058],'Box','on');
  set(g1(7),'Position',[0.4760    0.5563    0.0141    0.4029],'Box','on');
  set(g1(8),'Position',[0.0449    0.5552    0.4241    0.4040],'Box','on');
  Ejes_Visibles(1)
  % =======================================================================
  set(figure(2),'Position',[ 84   110   900   588],'Color','W')
  ko = 3e1;
  ym = max(PA(:));
  for i=1:4,subplot(2,2,i)
   hold on
   plot((Lx'-Lx(1,1)),Ly','-b','LineWidth',[3])
   plot((Lx_Cont'-Lx(1,1)),Ly_Cont','--g','LineWidth',[2])
   Mpa = mean(PA(:,:,i)')*ko;
   Mpa = (Mpa - min(Mpa))*ko/15;
   plot(T,Mpa,'-r','LineWidth',[6])
   plot(T,Mpa,'-k','LineWidth',[3])
   hold off
   if i < 4 title(cell2mat(['Estacion: ' Label_Station(i)])),axis xy,axis on
   else
    title('Todas las estaciones'),axis xy,axis on
   end
   xlim([T(1) T(end)])
   % ylim([0 max(PA(:))*1])
   ax1 = gca;                   % current axes
   ax1.GridLineStyle = '--';
   ax1.GridColor = [0 0 0];
   ax1.GridAlpha = 1;
   ax1.LineWidth = 0.5;
   ax1.YMinorGrid = 'on';
   ax1.XMinorGrid = 'on';
   grid
  end
  Ejes_Visibles(2)
  g = get(figure(2));
  g1 = g.Children;
  set(g1(1),'Position',[0.5385 0.0861 0.4180 0.3651],'Box','on');
  set(g1(2),'Position',[0.0414 0.0795 0.4477 0.3657],'Box','on');
  set(g1(3),'Position',[0.5347 0.5497 0.4286 0.3852],'Box','on');
  set(g1(4),'Position',[0.0419 0.5453 0.4467 0.3874],'Box','on');
  
 end
 % IMAGEN MEDIA DE LAS TRES ESTACIONES
 % ------------------------------------------------------------------------
 % PARAMETROS DE VISUALIZACION
 % -------------------------------
 h1 = 1;
 ko = 2e1;
 k1 = 200;
 k2 = 5;
 MA = [];
 % -------------------------------
 for i = 1:1
  MA = [MA; mean(PA(:,:,i)')];
 end
 % ------------------------------------------------------------------------
 Pro = prod(MA)/max(prod(MA))*ko;
 Sum = sum(MA)/max(sum(MA))*ko;
 Pp = permute(PA,[1 3 2]);Pp = permute(Pp,[2 1 3]);
 Ia = squeeze(mean(Pp))';
 Ia_mean_x = medfilt2(Ia,[n1 n2]);
 Ia_mean_y = mean(Ia);
 Ia_mean_y = smooth(Ia_mean_y - min(Ia_mean_y),h1)*k1;
 yv = linspace(0,24,length(Ia_mean_y));
 % ------------------------------------------------------------------------
 if band == 1
  set(figure(3),'Position',[ 47 59 1307 616],'Color','W')
  subplot(2,1,1)
  % SURF
  % -----
  [xa,ya] = meshgrid(linspace(T(1),T(end),size(Ia,2)),linspace(1,M_dh(end),size(Ia,1)));
  surf(xa,ya/Num_h,Ia_mean_x),shading interp,view([0 90])
  colorbar,colormap jet,alpha .7
  hold on
  plot((Ia_mean_y*k2-k1),yv,'-k','LineWidth',[5]),
  plot((Ia_mean_y*k2-k1),yv,'-g','LineWidth',[1],'MarkerSize',[4])
  hold off
  title(['Probabilidad de activación en estacion geomagnetica para gamma: [' num2str(Gamma_2(j_gamma)) ']'])
  ylabel('Intervalo horario (0-24 h)')
  ax1 = gca;                   % current axes
  ax1.YTick = [0:2:M_dh(end)];
  ax1.GridLineStyle = '--';
  ax1.GridColor = [0 0 0];
  ax1.GridAlpha = 1;
  ax1.LineWidth = 0.25;
  xlim([-k1 T(end)])
  ylim([0 24])
  subplot(2,1,2)
  plot((Lx'-Lxx(1,1)),Ly','-b','LineWidth',[2]), hold on
  plot((Lx_Cont'-Lxx(1,1)),Ly_Cont','--g','LineWidth',[3])
  plot(T,Pro,'-k','LineWidth',[6])
  plot(T,Pro,'-r','LineWidth',[2],'MarkerSize',[4])
  plot(T,Ia_mean_y/k2,'-k','LineWidth',[6])
  plot(T,Ia_mean_y/k2,'-m','LineWidth',[2])
  hold off
  xlabel(['Fecha inicial y final: ' datestr(round(f1)) ' --- hasta --- '...
   datestr(round(f2)) ' --- Tiempo total: [\approx' ...
   num2str(round((T(end)/365.25))) '] años --- [eje x en días]'])
  ylabel('Eventos (azul y verde) y precursor (rojo)')
  xlim([0 T(end)])
  g = get(figure(3));
  g1 = g.Children;
  set(g1(1),'Position',[0.0828 0.0727 0.8688 0.4611],'Box','on');
  set(g1(2),'Position',[0.9569 0.5836 0.0141 0.3679],'Box','on');
  set(g1(3),'Position',[0.0377 0.5838 0.9121 0.3643],'Box','on');
  set(g1(3),'XTick',[-200:200:T(end)])
  set(g1(1),'XTick',[0:200:T(end)])
  ax1 = gca; 
  ax1.YTick = [0:5:M_dh(end)];
  ax1.GridLineStyle = '--';
  ax1.GridColor = [0 0 0];
  ax1.GridAlpha = 1;
  ax1.LineWidth = 0.25;
  ax1.XAxisLocation = 'bottom';
  grid
  Ejes_Visibles(3)
  % =======================================================================
  INT_H = Int_H(IND_ACIERTOS);
  lx_ac = Int_H;
  ly_ac = Por_Ef*ones(1,length(Int_H));
  % =======================================================================
  set(figure(4),'Position',[49 136 1306 540],'Color','W')
  plot(yv,(Ia_mean_y),'-k','LineWidth',[8])
  hold on
  plot(yv,(Ia_mean_y),'*g','LineWidth',[4],'MarkerSize',[4])
  plot(yv,(Ia_mean_y),'.r','LineWidth',[2],'MarkerSize',[1])
  plot(lx_ac,ly_ac,'-r','LineWidth',[4])
  plot(lx_ac,ly_ac,'--g','LineWidth',[2])
  hold off
  title(['Probabilidad de activación en estacion geomagnetica para gamma: ['num2str(Gamma_2(j_gamma)) ']'])
  ylabel('Probabilidad de activación')
  xlabel('Intervalo horario (0-24 h)')
  xlim([-1 Dh(end)/Num_h+1])
  %ylim([0 1.1*max(Ia_mean_y)])
  g = get(4);
  g1 = g.Children;
  set(g1,'Position',[0.0415    0.1250    0.9514    0.8103],'Box','on')
  ax1 = gca;                    % current axes
  ax1.XTick = [-1:1:Dh(end)/Num_h];
  ax1.GridLineStyle = '--';
  ax1.GridColor = [0 0 0];
  ax1.GridAlpha = 1;
  ax1.LineWidth = 0.25;
  grid
  Ejes_Visibles(4)
  pause
 end
 % ========================================================================
 IA_MEAN_Y = [IA_MEAN_Y;Ia_mean_y'];
 IA_MEAN_X = [IA_MEAN_X;mean(Ia_mean_x')];
end
% SALIDAS
% ---------------------------------------------------------------------------------------------------------
save Precursor_Activation_Parameters IA_MEAN_Y IA_MEAN_X Lx Ly Lx_Cont Ly_Cont T k2 f1 f2 Gamma_2 Md
% ---------------------------------------------------------------------------------------------------------
% ---------------------------------------------------------------------------
fprintf('-------------------------------------------------------------------------\n')
disp('   La base datos de la fase de Produccion fue creada para las estaciones:  ')
disp(Label_Station')
fprintf('-------------------------------------------------------------------------\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%