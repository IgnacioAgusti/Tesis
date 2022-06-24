% Programa_General_Precursores_Relativistas.m
% ============================================
% DESCRIPCIÓN: prepara las componentes Bx, By, Bz en ventanas horarias que
% puedan ser calculadas en la función:
% Relativistic_Flux_as_Seismic_Precursors_Tensor
% Esta función tambien necesitará la posición geografica de las estaciones
% geomagneticas. 
% SALIDA: 
% CONTEO UMBRALIZADO Y FILTRADO PARA CADA GAMMA Y PARA CADA INTERVALO DE HORAS DEL DIA
% -------------------------------------------------------------------------
%         Count_Por % USANDO EFICIENCIA PORCENTUAL PARA CADA (Md) DIAS
%         P         % USANDO AMPLITUD DEL PRECURSOR
%         Fp % DENSIDAD DE MOMENTUM (VECTOR DE POINTING)
% =========================================================================
% DEFINICION DE PARAMETROS DE ENTRADA
% ======================================
% DIMENSION DE LOS DATOS
% ------------------------
[bn,bm,bk] = size(Bxn(:,1:mf,:));
% REDIMENSIONA LOS DATOS HASTA CONSEGUIR UN MULTIPLO DEL NUMERO DE VENTANAS X
% NUMERO DE DIAS POR VENTANA APROX A LA DIMENSION DE LOS DATOS (NUMERO TOTAL DE DIAS)
% -------------------------------------------------------------------------
n = 0;
while rem(bm-n,Md) > 0
 n = n + 1;
end
Nv = (bm - n)/Md;                 % NUMERO DE VENTANAS
% UMBRAL
% ---------
u1 = 0;                                 % UMBRAL INFERIOR (EN VALORES DE SIGMA)
u2 = 70;                                % UMBRAL SUPERIOR (EN VALORES DE SIGMA
n_u = 20;                                % NUMERO TOTAL DE UMBRALES
Uo = linspace(u1,u2,n_u);               % VECTOR DE UMBRALES (EN VALORES DE SIGMA)
% HORAS DEL DIA
% -----------------
nh = 1;                                 % NUMERO DE HORAS CORRESPONDIENTE AL ANCHO DE LA VENTANA DE DATOS
ng = 2;                                 % NUMERO DE GAMMAS USADOS PARA PROMEDIAR
St = [1 2 3];                             % ESTACIONES A UTILIZAR
[Gamma,E] = Energia_Vs_Gamma(gamma,0);         % ENERGIA EN GeV
% ------------------------
grados_km = 111.11;                            % CONVESRION DE GRADOS A KM
ko = 1;                                 % PARAMETRO DE AMPLIFICACION PARA VISUALIZAR
hp = 6;                                 % FILTRO DE SUAVIZACION PARA EL PRECURSOR (usando amplitud)
hc = 6;                                 % FILTRO DE SUAVIZACION PARA EL PRECURSOR (usando conteo (eficiencia porcentual %))
h2 = 14;                                % FILTRO DE SUAVIZACION PARA EL GRAFICO
h3 = 8;                                 % NUMERO DE BINES DEL HISTOGRAMA DE FRECUENCIA SISMICA

na = length(Fo)/365.5; % NUMERO TOTAL DE A??OS CONSIDERADO
% FRECUENCIA SISMICA
% -----------------------------
[ny,nx] = hist(Lxx(:,1),h3);
fs = round((na*365.25)/h3);  % FRECUECIA SISMICA EN SISMOS/TIEMPO
% =========================================================================
% ESTACIONES INTERMAGNET
% =========================================================================
load Estaciones_intermagnet_IEBQuakesExport.mat
latitud = 90 - Colatitute  ;        % 90 - COLATITUD = LATITUD
Long_Intermagnet = EastLongitude;
ind = find(Long_Intermagnet >= 180);
Long_Intermagnet(ind) = Long_Intermagnet(ind) - 360;
idx = [];
Ind_M = [];
for j=1:length(Label_Station);
 for i=1:length(IAGA)
  idx = strfind(IAGA(i),Label_Station(j));
  if (cell2mat(idx)== 1)
   Ind_M = [Ind_M,i];
   i = length(IAGA);
  end
 end
end
% ESTACIONES UTILIZADAS
% ====================================================================
% REALES
% -------------
x_station = latitud(Ind_M)'*grados_km;
y_station = Long_Intermagnet(Ind_M)'*grados_km;
z_station = [0.042 0.033 0.107]*grados_km;
% SIMULADAS
% -------------
%   x_station = x_station_Sim;
%   y_station = y_station_Sim;
%   z_station = z_station_Sim;
% -------------------------------
% LOCALIZACION DEL EVENTO SISMICO DE JAPON
% ====================================================================
switch Data_Stations
 case 1 %JAPON
  HLJ = 1.157401129603386e-05*60*60*9;% HORA LOCAL Japon
  Mag = 9.1;                     % MAGNITUD DEL EVENTO
  ts = '11-Mar-2011 05:46:23';   % GRAN TERREMOTO DE JAPON HORA UTC
  Lat_sismo = 38.19;
  Long_sismo = 142.22;
  x_event = Lat_sismo*grados_km;
  y_event = Long_sismo*grados_km;
  z_event = 29;
 case 2 %Europa silencio sismico sismos simulados
  HLJ = 1.157401129603386e-05*60*60*(1.3393);      % HORA LOCAL EUROPA
  Mag = 7.2; % MAGNITUD DEL EVENTO
  ts = datenum([2005 3 2]);
  Lat_sismo =  49;
  Long_sismo = 4.56;
  x_event = Lat_sismo*grados_km;
  y_event = Long_sismo*grados_km;
  z_event = 11.8;
  
 case 3 %USA
  %T=timezones('America');
  %      mean(T.UTCOffset)
  HLJ = 1.157401129603386e-05*60*60*(-4.9660);      % HORA LOCAL USA
  Mag = 6.8; % MAGNITUD DEL EVENTO
  ts = datenum([2001 02 28]);
  Lat_sismo =  47.15;
  Long_sismo = -122.63;
  x_event = Lat_sismo*grados_km;
  y_event = Long_sismo*grados_km;
  z_event = 56;
 case 4%CHILE
  %T=timezones('America');
  %      mean(T.UTCOffset)
  HLJ = 1.157401129603386e-05*60*60*(-4.9660);      % HORA LOCAL CHILE
  Mag = 8; % MAGNITUD DEL EVENTO
  ts = datenum([2007 08 15]);
  
  Lat_sismo =  -13.38;
  Long_sismo = -76.55;
  x_event = Lat_sismo*grados_km;
  y_event = Long_sismo*grados_km;
  z_event = 40;
 case 5%EUROPA DUR
  %T=timezones('America');
  %      mean(T.UTCOffset)
  HLJ = 1.157401129603386e-05*60*60*(1.3393);      % HORA LOCAL EUROPA
  Mag = 6.2; % MAGNITUD DEL EVENTO
  ts = datenum([2016 08 24]);
  
  Lat_sismo =  42.72;
  Long_sismo = 13.19;
  x_event = Lat_sismo*grados_km;
  y_event = Long_sismo*grados_km;
  z_event = 4.4;
 case 6%EUROPA  PAN IZN PEN
  %T=timezones('America');
  %      mean(T.UTCOffset)
  HLJ = 1.157401129603386e-05*60*60*(1.3393);      % HORA LOCAL EUROPA
  Mag = 6.6; % MAGNITUD DEL EVENTO
  ts = datenum([2017 07 20]);
  
  Lat_sismo =  36.92;
  Long_sismo = 27.41;
  x_event = Lat_sismo*grados_km;
  y_event = Long_sismo*grados_km;
  z_event = 7;
 case 7%JAPON PET PET PET
  %T=timezones('America');
  %      mean(T.UTCOffset)
  HLJ = 1.157401129603386e-05*60*60*9;      % HORA LOCAL JAPON
  Mag = 7.2; % MAGNITUD DEL EVENTO
  ts = datenum([2016 1 30]);
  
  Lat_sismo =  54.00;
  Long_sismo = 158.51;
  x_event = Lat_sismo*grados_km;
  y_event = Long_sismo*grados_km;
  z_event = 163.2;
  
end
% ====================================================================
% ETIQUETAS DE LOS PRECURSORES DEL TENSOR ELECTROMAGNETICO
% ---------------------------------------------------------
Precursor_Usado = {'Energia total','Presion de radiacion','Flujo de momentum en  x ==> plano yz)','Flujo de momentum en y ==> plano xz)',...
 'Flujo de momentum en z ==> plano xy)','Traza del tensor','Densidad de momentum (Vector de Pointing)'};
% ================================================================================================
if band == 1
 if band_display == 1
  fprintf('=====================status==========================\n')
  disp('Selected stations')
  fprintf('--------------------------\n')
  disp(Label_Station)
  fprintf('--------------------------\n')
  disp('Selected data')
  fprintf('--------------------------\n')
  disp(Data)
  fprintf('--------------------------\n')
  fprintf('Numero de dias por ventana Md: %1.0f\n',Md)
  fprintf('Numero de ventanas Nv: %1.0f\n',Nv)
  fprintf('Umbral inferior: %1.0f\n',u1)
  fprintf('Umbral superior: %1.0f\n',u2)
  fprintf('Numero de umbrales: %1.0f\n',n_u)
  Ih = ['Intervalo horario [0 to 24h] => [' repmat(' %d ',1,length(dh)) ']\n'];
  fprintf(Ih,dh)
  fprintf('Gamma inferior: %1.0f\n',gamma_1)
  fprintf('Gamma superior: %1.0f\n',gamma_2)
  fprintf('Valor medio de Gamma: %1.2f\n',mean(G))
  fprintf('Gammas usados: %1.0f\n',Ng)
  fprintf('Numero de gammas usados para promediar: %1.0f\n',ng)
  fprintf('Numero de estaciones usadas: %d\n',length(St))
  for i=1:length(St)
   fprintf('Estacion: [%1.1f] %s\n',St(i),cell2mat(IAGA(Ind_M(St(i))))),end
  fprintf('Parametro de amplificacion: %1.0f\n',ko)
  fprintf('Filtro de suavizacion precursor: %1.0f\n',hp)
  fprintf('Filtro de suavizacion del grafico: %1.0f\n',h2)
  fprintf('Precursor_Usado: %s',cell2mat(Precursor_Usado(Componente_Tensor_Maxwell)))
  fprintf('\n=====================status==========================\n')
 end
end
% GENERA PRECURSORES
% ==========================================================================================================================================
if band == 0
 % ==========================================================================================================================================
 M = [];
 Bxt = []; Byt = []; Bzt = [];
 
 % DETERMINA LAS DIMENSIONES MINIMAS DE LAS ESTACIONES
 % ====================================================
 for i=1:length(St)
  S = whos('-file',Data(St(i),:));
  M = [M;str2num(mat2str(S(2).size))];
 end
 % ----------------------------------------------------
 moo = min(M(:,1));
 noo = min(M(:,2));
 vm = 1:moo;
 vn = 1:noo;
 v = 1:moo*noo;
 [Bxn,Byn,Bzn,F,Fo] = Proceso_Normaliza_media0_desv1(Label_Station,Data,2,1);
 % ----------------------------------------------------
end
if band == 0
 % VISUALIZACION DE LA DATA DE LAS TRES ESTACIONES
 % =========================================================
 % COMPONENTE X
 % =============
 set(figure(1),'Position',[1 28 1920 950],'Color','W')
 for st=1:length(St)
  subplot(3,1,st)
  imagesc(Bxt(:,:,st)),colormap jet,colorbar
  title(['X component from Station number: ' num2str(St(st))]),grid, axis on
  ylabel('[0 to 24 hours]')
 end
 g = get(figure(1));
 g1 = g.Children;
 set(g1(1),'Position',[0.9607    0.1095    0.0139    0.2158]);
 set(g1(2),'Position',[0.0510    0.1100    0.8938    0.2157]);
 set(g1(3),'Position',[0.9629    0.4095    0.0139    0.2158]);
 set(g1(4),'Position',[0.0510    0.4096    0.8954    0.2157]);
 set(g1(5),'Position',[0.9633    0.7095    0.0139    0.2158]);
 set(g1(6),'Position',[0.0505    0.7093    0.8964    0.2157]);
 xlabel(['Time in days from: ' datestr(Fo(1)) '  to  ' datestr(Fo(end))])
 Ejes_Visibles(1)
 % COMPONENTE Y
 % =============
 set(figure(2),'Position',[1 28 1920 950],'Color','W')
 for st=1:length(St)
  subplot(3,1,st)
  imagesc(Byt(:,:,st)),colormap jet,colorbar
  title(['Y component from Station number: ' num2str(St(st))]),grid, axis on
  ylabel('[0 to 24 hours]')
 end
 g = get(figure(2));
 g1 = g.Children;
 set(g1(1),'Position',[0.9607    0.1095    0.0139    0.2158]);
 set(g1(2),'Position',[0.0510    0.1100    0.8938    0.2157]);
 set(g1(3),'Position',[0.9629    0.4095    0.0139    0.2158]);
 set(g1(4),'Position',[0.0510    0.4096    0.8954    0.2157]);
 set(g1(5),'Position',[0.9633    0.7095    0.0139    0.2158]);
 set(g1(6),'Position',[0.0505    0.7093    0.8964    0.2157]);
 xlabel(['Time in days from: ' datestr(Fo(1)) '  to  ' datestr(Fo(end))])
 Ejes_Visibles(2)
 % COMPONENTE X
 % =============
 set(figure(3),'Position',[1 28 1920 950],'Color','W')
 for st=1:length(St)
  subplot(3,1,st)
  imagesc(Bzt(:,:,st)),colormap jet,colorbar
  title(['Z component from Station number: ' num2str(St(st))]),grid, axis on
  ylabel('[0 to 24 hours]')
 end
 g = get(figure(3));
 g1 = g.Children;
 set(g1(1),'Position',[0.9607    0.1095    0.0139    0.2158]);
 set(g1(2),'Position',[0.0510    0.1100    0.8938    0.2157]);
 set(g1(3),'Position',[0.9629    0.4095    0.0139    0.2158]);
 set(g1(4),'Position',[0.0510    0.4096    0.8954    0.2157]);
 set(g1(5),'Position',[0.9633    0.7095    0.0139    0.2158]);
 set(g1(6),'Position',[0.0505    0.7093    0.8964    0.2157]);
 xlabel(['Time in days from: ' datestr(Fo(1)) '  to  ' datestr(Fo(end))])
 Ejes_Visibles(3)
 % ===========================================================================================
end
% GENERACION DE LA VENTANA DESLIZANTE
% =============================================
NH = Num_h*24*dho; % NUMERO TOTAL DE HORAS EN UN DIA
np = round(NH/nh); % NUMERO TOTAL DE PARTES IGUALES EN LA SE DESEA DIVIDIR LA TOTALIDAD DE LOS DATOS

% ENCUENTRA MULTIPLO DE LA MATRIZ DE ACUERDO AL NUMERO DE DIAS POR VENTANA ELEGIDOS
% ==================================================================================
[m,n,k] = size(Bxn);
dn = 0;
while rem((m+dn)/nh,np)
 dn = dn + 1;
end
% ---------------------------------
vn = m-dn+1:m;
Bmx = Bxn(vn,:,:);                % MATRIZ MASCARA PARA ASEGURAR UNA DIMENSION DE MATRIZ QUE SEA MULTIPLO DE nh
Bmy = Byn(vn,:,:);                % PARA ELLO SE TOMAN LOS ULTIMOS dn ELEMENTOS DE LA MATRIZ Bxn,Byn y Bzn
Bmz = Bzn(vn,:,:);
BxN = [Bxn;Bmx];                  % SE AGREGAN A LA MATRIX Bxn,Byn y Bzn CREANDO UNA NUEVA MATRIZ BxN,ByN y BzN
ByN = [Byn;Bmy];
BzN = [Bzn;Bmz];
% VISUALIZACION DE LA MATRIZ MASCARA
% =============================================================================================
if band == -3
 close all
 set(figure(1),'Position',[5 261 1912 717],'Color','W')
 BxNv = BxN;
 BxNv(end-dn+1:end,:,:) = BxNv(end-dn+1:end,:,:)+10;
 imagesc(BxNv(:,:,1)),axis xy
 title(['Original matrix [' num2str(size(Bxn,1)) 'x' num2str(size(Bxn,2)), '] plus the mask ' num2str(dn) ...
  ' rows at the end: New matriz: [' num2str(size(BxN,1)) 'x' num2str(size(BxN,2)), ']' ])
 ylim([1 size(BxN,1)+50])
 xlabel('Time in days')
 ylabel('Time interval (0 to 24 hours)')
 axis on,grid
 colormap jet,colorbar
 g = get(figure(1));
 g1 = g.Children;
 set(g1(1),'Position',[0.9582    0.0823    0.0139    0.8466]);
 set(g1(2),'Position',[0.0460    0.0823    0.8950    0.8466]);
 Ejes_Visibles(1)
 pause
end
% =============================================================================================
% ------------------------------------
mo = round((m+dn)/np);
Bxnr = reshape(BxN,np,mo,n,k);
Bynr = reshape(ByN,np,mo,n,k);
Bznr = reshape(BzN,np,mo,n,k);
% PARA VISUALIZAR EL FORMATO DE LOS DATOS
% PERMITE TESTAR EL FORMATO DE LOS DATOS
% ================================================================================================
if band == -3
 N_h = linspace(0,24,size(BxN,1));
 v1 = 1:np:size(BxN,1);v1(1)=1;
 v2 = np:np:size(BxN,1);
 Bv = BxN;
 nk = 10;     % NUMERO DE FILAS A TESTAR
 close all
 set(figure(1),'Position',[6 32 1910 946],'Color','W')
 for i=1:size(Bxnr,2)
  Bc = BxN;
  a = squeeze(Bv(1:nk,:,1))';
  b = squeeze(Bxnr(1:nk,i,:,1))';
  [Xp,Yp,p] = polinomio2(a(:),b(:),1,10,0);
  Test = sum(a(:)-b(:));
  if band_dysplay == 1
   disp(Test)
  end
  % ----------------------------------------
  subplot(2,2,1)
  plot(a(:),b(:),'.'),hold on
  plot(Xp,Yp,'-r')
  grid
  title(['Mean interval: ' num2str(mean(N_h(v1(i)):N_h(v2(i)))) ' hours and [mo : bo] ==> [' num2str(p(1)) ' : ' num2str(round(p(2))) ']' ])
  xlabel('Old format'),ylabel('New format')
  xlim([-max(abs(a(:)))-5 max(abs(a(:)))+5])
  ylim([-max(abs(a(:)))-5 max(abs(a(:)))+5])
  subplot(2,2,2)
  plot(squeeze(Bv(1:nk,:,1))','-b')
  hold on,plot(squeeze(Bxnr(1:nk,i,:,1))','--r'),hold off,grid
  title(['Old format (-b) and New format (-r) for first: ' num2str(nk) ' rows'])
  % -------------------------------------------------
  Bv(1:size(Bxnr,1),:,:) = [];
  Bc(v1(i):v2(i),:,1) = Bc(v1(i):v2(i),:,1) + 10;
  % -------------------------------------------------
  xlim([0 size(Bxnr,3)])
  ylim([-max(abs(a(:)))-5 max(abs(a(:)))+5])
  subplot(4,1,3)
  imagesc(1:size(Bxnr,3),N_h,Bc(:,:,1)),axis xy,axis on,grid,colormap jet,colorbar
  title(['Horary interval: ' num2str(N_h(v1(i))) ' to ' num2str(N_h(v2(i)))])
  ylabel('0 to 24h')
  subplot(4,1,4)
  plot(mean(Bc(v1(i):v2(i),:,1))),grid
  xlim([0 size(Bxnr,3)])
  ylim([0 20])
  xlabel('Time in days')
  ylabel('Amplitude')
  g = get(figure(1));
  g1 = g.Children;
  set(g1(1),'Position',[0.0393 0.0698 0.9382 0.1871]);
  set(g1(2),'Position',[0.9644 0.3277 0.0140 0.1808]);
  set(g1(3),'Position',[0.0403 0.3277 0.9047 0.1808]);
  set(g1(4),'Position',[0.5241 0.5740 0.4529 0.3795]);
  set(g1(5),'Position',[0.0393 0.5719 0.4550 0.3784]);
  Ejes_Visibles(1)
  pause(.01)
 end
 pause
end
% =========================================================================
% ------------------------------------
% VALOR MEDIO POR VENTANA DESLIZANTE (nh) EN LA DIRECCION VERTICAL
% ================================================================
Bxm = squeeze(mean(permute(Bxnr,[2 1 3 4])));
Bym = squeeze(mean(permute(Bynr,[2 1 3 4])));
Bzm = squeeze(mean(permute(Bznr,[2 1 3 4])));
Bxm = Bxm(:,:,St);
Bym = Bym(:,:,St);
Bzm = Bzm(:,:,St);
% VISUALIZA LOS PROMEDIO DE CADA VENTANA
% ----------------------------------------
if band == -3
 close all
 for i=1:size(Bxm,3)
  set(figure(i),'Position',[5 32 1912 946],'Color','W')
  subplot(3,1,1)
  plot(Bxm(:,:,i)')
  ylim([-10 10])
  xlim([0 size(Bxm,2)])
  ylabel('Amplitude')
  title(['Bxm and Station: ' Label_Station(i)])
  grid
  subplot(3,1,2)
  plot(Bym(:,:,i)')
  ylim([-10 10])
  xlim([0 size(Bxm,2)])
  ylabel('Amplitude')
  title(['Bym and Station: ' Label_Station(i)])
  grid
  subplot(3,1,3)
  plot(Bzm(:,:,i)')
  ylim([-10 10])
  xlim([0 size(Bxm,2)])
  ylabel('Amplitude')
  xlabel('Time in days')
  title(['Bzm and Station: ' Label_Station(i)])
  grid
  g = get(figure(i));
  g1 = g.Children;
  set(g1(1),'Position',[0.0387 0.1100 0.9419 0.2157]);
  set(g1(2),'Position',[0.0345 0.4096 0.9440 0.2157]);
  set(g1(3),'Position',[0.0314 0.7093 0.9461 0.2157]);
  Ejes_Visibles(i)
 end
 pause
end
% ==================================================================================================
% GENERACION DEL PRECURSOR (CICLO EN GAMMA Y EN CADA UMBRAL Y PARA LAS HORAS DEL DIA ELEGIDAS)
% ==================================================================================================
GM = {};
GS = {};
% CICLOS PARA CADA INTERVALO HORARIO Y PARA CADA VALOR DE GAMMA
% -------------------------------------------------------------------
Gm = [];
Gs = [];
P = [];
fp = [];
for k=1:length(dh)-1 % ojo (-1)         % PARA CADA INTERVALO HORARIO
 for g = 1:length(gamma)-ng          % PARA CADA VALOR DE GAMMA
  Gamma = gamma(g:g+ng);         % SE CONSIDERAN DOS VALORES CONSECUTIVOS DE GAMMA (LUEGO SE PROMEDIA)
  Energy = E(g:g+ng);
  Pm = []; pm = [];
  Ps = []; ps = [];
  for j=1:length(Gamma)
   M_Gamma = [M_Gamma;mean(Gamma)];  % VALORS MEDIOS DE GAMMA POR INTERVALO PROMEDIADO DEL FLUJO
   % PRECURSOR USADO:TENSOR ELECTROMAGNETICO
   % -----------------------------------------------------------------------------------------------------------------------------------------------------
   % TENSOR DE MAXWELL
   % -------------------
   if Componente_Tensor_Maxwell ~= 7
    [W,Rp,Sxy,Sxz,Syz,Tz] = Relativistic_Flux_as_Seismic_Precursors_Tensor(Bxm,Bym,Bzm,...   % COMPONENTES PROMEDIO DEL CAMPO POR INTERVALO
     Gamma(j),...                                                     % VALOR DE GAMMA
     x_station,y_station,z_station,...                                % COORDENADAS DE POSICION DE LA ESTACION
     x_event,y_event,z_event,...                                      % COORDENADAS DE EVENTO SISMICO
     Mag,0);%pause(0.5)                                                       % MAGNITUD CONSIDERADA
   end
   k_w = 1;%1e11;
   k_p = 1;%1e7;
   k_syz = 1;%1e3;
   k_sxz = 1;%1e3;
   k_sxy = 1;%1e1;
   k_tz = 1;%1e5;
   k_fp = 1;%0.9882;%1e1;
   
   % ETIQUETAS DE LOS PRECURSORES DEL TENSOR ELECTROMAGNETICO
   % ---------------------------------------------------------
   switch  Componente_Tensor_Maxwell
    case 1                             % ENERGIA TOTAL
     Fp = W/k_w;
    case 2                             % PRESION DE RADIACION
     Fp = Rp/k_p;
    case 3
     Fp = abs(Syz)/k_syz;              % FLUJO DE MOMENTO (EN X ==> plano yz)
    case 4
     Fp = abs(Sxz)/k_sxz;              % FLUJO DE MOMENTO (EN Y ==> plano xz)
    case 5
     Fp = abs(Sxy)/k_sxy;              % FLUJO DE MOMENTO (EN Z ==> plano xy)
    case 6                             % TRAZA DEL TENSOR
     Fp = abs(Tz)/k_tz;
    case 7                             % DENSIDAD DE MOMENTUM (VECTOR DE POINTING)
     % --------------------------------------------------------------------------------------------------------------------------------------------
     % DENSIDAD DE MOMENTUM (VECTOR DE POINTING)
     % ----------------------------------------------
     Fp = Relativistic_Flux_as_Seismic_Precursors(Bxm,Bym,Bzm,...                   % COMPONENTES PROMEDIO DEL CAMPO POR INTERVALO
      Gamma(j),...                                                             % VALOR DE GAMMA
      x_station,y_station,z_station,...                                        % COORDENADAS DE POSICION DE LA ESTACION
      x_event,y_event,z_event,...                                              % COORDENADAS DE EVENTO SISMICO
      Mag,0);%pause                                                                % MAGNITUD DEL EVENTO
     Fp = Fp/k_fp;
    otherwise
     disp(['Warning ... Componente del tensor NO seleccionada'])
   end
   % PARA VISUALIZAR EL PRECURSOR PARA CADA ESTACION
   % -------------------------------------------------
   if band == -3
    close all
    M_Fp = max(Fp(:));
    set(figure(1),'Position',[6 32 1910 946],'Color','W')
    subplot(3,1,1)
    plot(Fp(:,:,1)')
    ylim([0 M_Fp])
    xlim([0 size(Fp,2)])
    ylabel('Amplitude')
    title(['Fp Station: ' Label_Station(1) ])
    grid
    subplot(3,1,2)
    plot(Fp(:,:,2)')
    ylim([0 M_Fp])
    xlim([0 size(Fp,2)])
    ylabel('Amplitude')
    title(['Fp and Station: ' Label_Station(2)])
    grid
    subplot(3,1,3)
    plot(Fp(:,:,3)')
    ylim([0 M_Fp])
    xlim([0 size(Fp,2)])
    ylabel('Amplitude')
    xlabel('Time in days')
    title(['Fp and Station: ' Label_Station(3)])
    grid
    g = get(figure(1));
    g1 = g.Children;
    set(g1(1),'Position',[0.0419 0.1100 0.9361 0.2157]);
    set(g1(2),'Position',[0.0408 0.4096 0.9377 0.2157]);
    set(g1(3),'Position',[0.0414 0.7093 0.9361 0.2157]);
    Ejes_Visibles(1)
    pause(.1)
   end
   % --------------------------------------------------------------------------------------------------------------------------------------------
   % PROMEDIO DEL FLUJO POR ESTACION EN EL INTERVALO DE HORAS ELEGIDO
   % ------------------------------------------------------------------
   if length(Station_Used_Now) > 1
    for i=1:length(Station_Used_Now)
     pm = [pm;mean(Fp(dh(k):dh(k+1),:,St(Station_Used_Now(i))))];
     ps = [ps; std(Fp(dh(k):dh(k+1),:,St(Station_Used_Now(i))))];
     % VISUALIZA EL PROMEDIO DE Fp POR INTERVALO HORARIO Y POR ESTACION
     % -------------------------------------------------------------------
     if band == -3
      %close all
      set(figure(1),'Position',[5 121 1912 857],'Color','W')
      M_Fp = max(Fp(:));
      subplot(4,1,1)
      plot(Fp(dh(k):dh(k+1),:,St(Station_Used_Now(i)))')
      title(['Fp for interval: [ ' num2str(dh(k)) ' to ' num2str(dh(k+1)) ' ] and station ' num2str(St(Station_Used_Now(i))) ' and Gamma: ' num2str(Gamma(j))])
      ylabel('Amplitude')
      ylim([0 M_Fp])
      xlim([0 size(Fp,2)])
      grid
      subplot(4,1,2)
      plot(mean(Fp(dh(k):dh(k+1),:,St(Station_Used_Now(i)))),'r'),grid
      title(['Fp mean (r) value '])
      ylim([0 M_Fp])
      xlim([0 size(Fp,2)])
      subplot(4,1,3)
      plot(std(Fp(dh(k):dh(k+1),:,St(Station_Used_Now(i)))),'-b')
      title('Fp std (b) value ')
      ylabel('Amplitude')
      ylim([0 M_Fp])
      xlim([0 size(Fp,2)])
      grid
      subplot(4,1,4)
      g = get(figure(1));
      g1 = g.Children;
      set(g1(1),'Position',[0.0497 0.1100 0.9372 0.1544]);
      set(g1(2),'Position',[0.0486 0.3459 0.9362 0.1320]);
      set(g1(3),'Position',[0.0476 0.5650 0.9362 0.1320]);
      set(g1(4),'Position',[0.0460 0.7841 0.9367 0.1320]);
      Ejes_Visibles(1)
      pause(.1)
     end
     % -------------------------------------------------------------------
     
    end
   else
    pm = [pm;mean(Fp(dh(k):dh(k+1),:,St(Station_Used_Now)))];
    ps = [ps; std(Fp(dh(k):dh(k+1),:,St(Station_Used_Now)))];
   end
  end
  Pm = prod(pm)/std(prod(pm));            % PROBABILIDAD CONDICIONADA: 1) PARA CADA ESTACION,2) PARA CADA INTERVALO HORARIO Y
  Ps = prod(ps)/std(prod(ps));            % 3) PARA CADA VALOR DE GAMMA (SE TOMA EL PRODUCTO)
  % -----------------------------------------------------------------------------------------------------------------------------------------------
  if band == -3
   % VISUALIZA LA PROBABILIDAD CONDICIONADA
   % ------------------------------------------
   plot(Pm,'-r'),hold on,plot(Ps,'-b'),hold off
   title('Conditioned probability Pm (-r) and Ps (-b)')
   xlabel('Time in days')
   ylabel('Amplitude')
   ylim([0 M_Fp])
   xlim([0 size(Fp,2)])
   set(g1(4),'Position',[0.0460 0.7841 0.9367 0.1320]);
   grid
   Ejes_Visibles(1)
   pause
   % ------------------------------------------
  end
  % ------------------------------------------------------------------------------------------
  fp = [fp;Pm + Ps];      % PRECURSOR USADO (DESVIACION MAS LA MEDIA DEL FLUJO)
  %fp = [fp;Pm];           % PRECURSOR USADO (DESVIACION DEL FLUJO)
  % ------------------------------------------------------------------------------------------
  if band == -3
   close all
   set(figure(1),'Position',[5 660 1852 318],'Color','W')
   plot(Pm,'-r'),hold on,plot(Ps,'--g'),plot(Pm+Ps,'-b'),hold off
   title('Conditioned probability: Pm (-r), Ps (-g)  and Pm + Ps (-b)')
   xlabel('Time in days')
   ylabel('Amplitude')
   ylim([0 M_Fp])
   xlim([0 size(Fp,2)])
   grid
   g = get(figure(1));
   g1 = g.Children;
   set(g1(1),'Position',[0.0416 0.1447 0.9438 0.7358]);
   Ejes_Visibles(1)
   pause
  end
 end
 % -------------------------
 Gm = [Gm;mean(fp)*ko];   % VALOR MEDIO PARA LOS VALORES ACTUALES DE GAMMA
 Gs = [Gs;std(fp)*ko];    % DESVIACION PARA LOS VALORES ACTUALES DE GAMMA
 % -------------------------
end
%          % CONTADOR DE EVENTOS (PARA CADA VALOR DE GAMMA)
% size(GM(Station).GM)
% ======================================================================
P = [];
Count_Por = []; % CONTEO PORCENTUAL POR CADA VENTANA DE Md DIAS
for go = 1:length(Uo)
 uoo = Uo(go);                           % UMBRAL ACTUAL
 if length(Station_Used_Now)>1
  a = mean(Gm) + mean(Gs);
 else
  a = Gm + Gs;
 end
 % UMBRALIZAMOS EL PRECURSOR
 % LO DIVIDIMOS EN Nv VENTANAS DE Md DIAS CADA UNA
 % ------------------------------------------------
 a = a(1:Md*Nv)/std(a(1:Md*Nv));  % NORMALIZAMOS A SU DESVIACION STANDARD (LA ESCALA SERA EN STD)
 Max_a = [Max_a;max(a)];          % GUARDAMOS SU VALOR MAXIMO
 b = reshape(a,Nv,Md);            % DIVIDIMOS EN Nv VANTANAS DE Md DIAS CADA UNA
 v = 1:Md:Md*Nv;
 count = zeros(1,Nv);
 count_por = zeros(1,Nv);
 to = Fo(1,1:Md*Nv);
 so = std(a);
 Ud = uoo*so;
 % UMBRAL USADO (PROPORCIONAL A LA DESVIACION STANDARD)
 % -----------------------------
 % CONTEOS MAYORES AL UMBRAL
 % ----------------------------
 if band == -3
  Bb = [];
  set(figure(1),'Position',[5 545 1852 393],'Color','W'),end
 for q=1:Nv
  s = find(b(q,:)>= Ud);
  count(q) = sum(b(q,s));        % PROPORCIONAL A LA AMPLITUD (suma de las amplitudes de las activaciones)
  count_por(q) = (length(s)/Md)*100; % ACTIVACION PORCENTUAL POR CADA Md DIAS
  % -----------------------------------------------------------------
  if band == -3
   Bb = [Bb b(q,:)];
   ly_u = Ud*ones(1,length(Bb));
   lx_u = linspace(1,length(Bb),length(Bb));
   subplot(2,1,1)
   plot(Bb),hold on,plot(lx_u,ly_u,'-r'),hold off
   title('Conditioned probability Gm + Gs')
   ylabel('Amplitude')
   ylim([-2 1.2*Ud])
   xlim([0 size(Fp,2)])
   grid
   subplot(2,1,2)
   plot(count)
   title(['Count for Ud: ' num2str(Ud)])
   xlabel('Time in days')
   ylabel('Counts')
   g = get(figure(1));
   g1 = g.Children;
   set(g1(1),'Position',[0.0373 0.1323 0.9449 0.3189]);
   set(g1(2),'Position',[0.0356 0.6296 0.9460 0.2954]);
   grid
   Ejes_Visibles(1)
   pause(.01)
  end
  % -------------------------------------------------------------------------
 end
 % ------------------------------------------------------
 p = smooth(smooth(count,hp),hp/2)*ko;        % SE FILTRA EL PRECURSOR
 % ------------------------------------------------------
 count_por = smooth(smooth(count_por,hc),hc/2)*ko;    % SE FILTRA EL PRECURSOR
 % ------------------------------------------------------
 % CONTEO UMBRALIZADO Y FILTRADO PARA CADA GAMMA Y PARA CADA INTERVALO DE HORAS DEL DIA
 % ------------------------------------------------------------------------------------
 P = [P;p'];                        % USANDO AMPLITUD DEL PRECURSOR
 Count_Por = [Count_Por;count_por'];% USANDO EFICIENCIA PORCENTUAL PARA CADA (Md) DIAS
 
end
% PARA VISUALIZAR LAS
% -----------------------------------------------------------------
if band == -3
 close all
 M_a = round(max(abs(a(:)))+1);
 set(figure(1),'Position',[5 655 1854 323],'Color','W')
 jv = 0:size(b,2);
 tv = 1:Md*Nv;
 hold on
 for i=1:Md
  bx = i*Nv*ones(1,Nv);%tv((jv(i))*Nv+1:i*Nv);
  by = linspace(-M_a,M_a,Nv);
  plot(tv((jv(i))*Nv+1:i*Nv),b(:,i),'.-','LineWidth',[3])
  plot(bx,by,'-b','LineWidth',[1])
 end
 plot(tv,a,'--c')
 xlim([0 size(Fp,2)])
 ylim([-M_a M_a])
 xlabel('Time in days')
 ylabel('Amplitude')
 title(['Format verification for: [' num2str(Md) '] windows each one with [' num2str(Nv) '] total [' num2str(Nv*Md) '] days'])
 hold off
 g = get(figure(1));
 g1 = g.Children;
 set(g1(1),'Position',[0.0345 0.1455 0.9509 0.7214]);
 grid
 Ejes_Visibles(1)
 pause
end
% -----------------------------------------------------------------
A = a*max(Mag_S)*1.5;             % SE??AL AMPLIFICADA RELATIVA A LA MAGNITUD MAXIMA (SOLO PARA VISUALIZAR)
% -----------------------------------------------------------------------------------------------------------
% VISUALIZA EL STATUS ACTUAL DEL PRECURSOR Y LA POSICION DE LOS EVENTOS SISMICOS
% -----------------------------------------------------------------------------------------------------------
if band == 1
 set(figure(1),'Position',[5 509 1912 469],'Color','W')
 Ea = num2str(Energy(j));
 f1 = Fo(1,1);
 x = [datenum(ts)-f1 datenum(ts)-f1];
 y = [0 50];
 bar(nx-f1,ny/10,'FaceColor',[0 1 0.1],'EdgeColor','b','BarWidth',[0.1]),hold on  % FRECUENCIA SISMICA
 plot(to-f1,A,'-b','LineWidth',[2])
 plot(Lxx'-f1,Lyy','-r','LineWidth',[1])
 line(x,y,'Color','g','LineStyle','--','LineWidth',[2])
 title(['dh: ' num2str(dh(1)) ' to ' num2str(dh(end)) ' Hours and Gamma: ' num2str(round(Gamma(j))) ' Energy '  Ea(1:4) ...
  ' (GeV)  == Seismic frequency each [' num2str(fs) '] days == Station: [ ' cell2mat(Label_Station(Station_Used_Now)) ' ]'])
 hold off
 xlabel(['Time in days from: ' datestr(round(Fo(1,1))) ' ==to== ' datestr(round(Fo(1,end))) ' [' num2str(round(na)) '] years' ])
 ylabel('Seismic amplitude')
 xlim([Fo(1,1)-f1 Fo(1,end)-f1])
 grid
 g = get(1);
 g1 = g.Children;
 set(g1,'Position',[0.0335    0.1250    0.9514    0.8103])
 Ejes_Visibles(1)
end
% =================================================
% MOMENTO ESTADISTICO UTILIZADO
% --------------------------------
switch Momento_Usado
 case 1         % MEDIA
  As = mean(P);             % USANDO AMPLITUD DEL PRECURSOR
 case 2         % KURTOSIS
  As = kurtosis(P);
  %  As = kurtosis(Count_Por);
  for i=1:length(As)
   ind_k = findstr(num2str(As(i)),'NaN');
   if length(ind_k)>0
    As(i) = 0;
   end
  end
 case 3
  As = std(P); % DESVIACION
 otherwise
  disp(['Momento estadistico no seleccionado'])
end
% ==================
% VISUALIZACION
% ==================
[mt,nt] = size(P);
t1 = Fo(1,1);
t2 = Fo(1,end);
t = linspace(t1,t2,nt);
% -------------------------------------------------------------------------
if band == 1
 % PLOT
 % ====
 set(figure(2),'Position',[5 514 1912 464],'Color','W')
 bar(nx-f1,ny/10,'FaceColor',[0 1 0.1],'EdgeColor','b','BarWidth',[0.1])  % FRECUENCIA SISMICA
 hold on
 plot(Lxx'-f1,Lyy','-r','LineWidth',[1])
 plot(t-f1,As,'.-b','LineWidth',[2])
 plot(datenum(ts)-f1,Mag,'*m')
 x = [datenum(ts)-f1 datenum(ts)-f1];
 y = [0 Mag];
 line(x,y,'Color','k','LineStyle','--','LineWidth',[2])
 xlabel(['Time in days from: ' datestr(round(Fo(1,1))) ' ==to== ' datestr(round(Fo(1,end))) '  nyears:  [' num2str(round(na)) ']' ])
 ylabel('Seismic amplitude')
 title(['dh: ' num2str(dh(1)) ' to ' num2str(dh(end)) ' Hours === Gamma from: ' num2str(round(gamma(1))) ' to ' num2str(round(gamma(end))) ' == Energy ' num2str(round(E(1)))...
  ' to ' num2str(round(E(end))) ' (Gev)  == Seismic frequency each [' num2str(fs) '] days == Station: [ ' cell2mat(Label_Station(Station_Used_Now)) ' ]'])
 xlim([Fo(1,1)-f1 Fo(1,end)-f1])
 grid
 hold off
 g = get(2);
 g2 = g.Children;
 set(g2,'Position',[0.0335    0.1250    0.9514    0.8103])
 Ejes_Visibles(2)
 % SURF P
 % ========
 eo = 10;
 set(figure(3),'Position',[5 32 921 906],'Color','W')
 surf(t-f1,Uo+Uo(1),medfilt2(P,[1 1])),shading interp,colormap jet,hold on
 plot(Lxx'-f1,Lyy'-eo,'-r','LineWidth',[2])
 line(x,y-eo,'Color','k','LineStyle','-.','LineWidth',[2])
 hold off
 xlim([Fo(1,1)-f1 Fo(1,end)-f1])
 ylim([-eo Uo(end)+Uo(1)])
 xlabel(['Time in days from: ' datestr(round(Fo(1,1))) ' ==to== ' datestr(round(Fo(1,end))) '  nyears:  [' num2str(round(na)) ']' ])
 ylabel('Magnitudes (<0) and threshold in sigma values')
 zlabel('Precursors amplitude')
 view([0 90])
 title(['dh: ' num2str(dh(1)) ' to ' num2str(dh(end)) ' Hours and Precursor amplitude for ' cell2mat([Label_Station(Station_Used_Now')] )])
 Ejes_Visibles(3)
 colorbar
 % SURF Count_Por
 % ================
 eo = 10;
 set(figure(4),'Position',[934 32 956 904],'Color','W')
 surf(t-f1,Uo+Uo(1),medfilt2(Count_Por,[1 1])),shading interp,colormap jet,hold on
 plot(Lxx'-f1,Lyy'-eo,'-r','LineWidth',[2])
 line(x,y-eo,'Color','k','LineStyle','-.','LineWidth',[2])
 hold off
 xlim([Fo(1,1)-f1 Fo(1,end)-f1])
 ylim([-eo Uo(end)+Uo(1)])
 xlabel(['Time in days from: ' datestr(round(Fo(1,1))) ' ==to== ' datestr(round(Fo(1,end))) '  nyears:  [' num2str(round(na)) ']' ])
 ylabel('Magnitudes (<0) and threshold in sigma values')
 zlabel('Precursors amplitude')
 view([0 90])
 title(['dh: ' num2str(dh(1)) ' to ' num2str(dh(end)) ' Hours and Pocentual (%) for ' cell2mat([Label_Station(Station_Used_Now')] )])
 Ejes_Visibles(4)
 colorbar
 pause(.1)
end
% FIN
