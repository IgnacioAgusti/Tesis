% Verifica_Eficiencia_Para_Varios_Parametros.m
%
% VERIFICA EFICIENCIA DE LOS PRECURSORES RELATIVISTAS
% ===================================================
% DESCRIPCIÓN: Programa transforma componentes Bx, By y Bz a estructuras
% con salidas que contienen el precursor sismico.
% SALIDA: STADISTIC_Results
%==========================================================================
clear all
close all
% =========================================================================
Gamma_2 = linspace(1,100,40); % Gammas que se quieren calcular
Gamma_2(1) = 1;
% =========================================================================
stations= [7]; % Estaciones que se quieren procesar 1 = Japon, 2 = Europa, 3= USA, etc...
band_display = 1; % Muestra en pantalla el proceso

for Data_Stations= stations
 Stadistic_Entrena  = {};
 Stadistic_Produce = {};
 for CASO = 0:1 % (0) PARA ENTRENAMIENTO (1) PARA PRODUCCION
  Caso_Ent_Prod = CASO;
  [Data,Label_Station] = Data_Geomagnetica(Data_Stations);
  
  % CARGAR UN DETERMINADO CONJUNTO DE DATOS GEOMAGNETICOS
  % =================================================================
  if Caso_Ent_Prod == 0
   switch Data_Stations
    case 1
     load ('Entrena_Etiqueta_Eventos_Japan_II_6M_2000_2019','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     disp('Fase de Entrenamiento JAPON')
    case 12
     load ('Etiquetas_KNY_Entrenamiento_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     disp('Fase de Entrenamiento JAPON(SOLO KNY)')
    case 13
     load ('Etiquetas_MMB_Entrenamiento_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     disp('Fase de Entrenamiento JAPON(SOLO MMB)')
    case 2
     % GENERA ETIQUETAS ENTRENAMIENTO
     % [Parametros_entrada, Sismos_Simulados ,Tabla_generada] = Simula_Terremotos(label_station_central...
     %  ,numero_terremotos,mag_mi,mag_max,fo1,fo2, Prof1, Prof2,radio_inferior,radio_superior,band);
     % [Parametros_entrada, Sismos_Simulados ,Tabla_generada]= Simula_Terremotos({'DOU'} ...
     %   , 20,6,8,'2005-01-01','2019-01-01',5,50,10,200,1);
     %  % AGREGAR ETIQUETAS DE SIMULACION EUROPA SILENCIO
     % load ('Etiquetas_DOUEntrenamiento_Ready.mat','Lx_c','Lx','Ly_c'...
     % ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     % load ('Etiquetas_DOUEntrenamiento2_Ready.mat','Lx_c','Lx','Ly_c'...
     % ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     % load ('Etiquetas_DOUEntrenamiento3_Ready.mat','Lx_c','Lx','Ly_c'...
     %   ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     % load ('Etiquetas_DOUEntrenamiento4_Ready.mat','Lx_c','Lx','Ly_c'...
     % ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     %                     disp('Fase de Entrenamiento Europa silencio sismico, con sismos simulados en sitios aleatorios')
    case 3
     load ('Etiquetas_USA_Entrenamiento_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     disp('Fase de Entrenamiento USA')
    case 4
     load ('Etiquetas_CHILE_Entrenamiento_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     disp('Fase de Entrenamiento CHILE')
    case 5
     load ('Etiquetas_EUROPA_DUREntrenamiento_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     disp('Fase de Entrenamiento DUR')
    case 6
     load ('Etiquetas_EUROPA_PEGEntrenamiento_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     disp('Fase de Entrenamiento PEG')
    case 7
     load ('Etiquetas_PETEntrenamiento_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     % Solo del 2007 al 2011 con 10 sismos entrenando y todo
     % el resto produciendo
     % load ('Etiquetas_PETEntrenamiento2_Ready.mat','Lx_c','Lx','Ly_c'...
     %  ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE ENTRENAMIENTO
     disp('Fase de Entrenamiento PET')
     
   end
   disp('Fase de Entrenamiento')
  else
   switch Data_Stations
    case 1
     load ('Produce_Etiqueta_Eventos_Japan_II_6M_2000_2019','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion JAPON')
     
    case 12
     load ('Etiquetas_KNY_Produccion_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion JAPON(SOLO KNY)')
    case 13
     load ('Etiquetas_MMB_Produccion_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion JAPON(SOLO MMB)')
    case 2
     % ETIQUETAS DE SIMULACION EUROPA SILENCIO
     % load ('Etiquetas_DOUProduccion_Ready.mat','Lx_c','Lx','Ly_c'...
     %  ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     % load ('Etiquetas_DOUProduccion2_Ready.mat','Lx_c','Lx','Ly_c'...
     %  ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     % load ('Etiquetas_DOUProduccion3_Ready.mat','Lx_c','Lx','Ly_c'...
     %  ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     load ('Etiquetas_DOUProduccion4_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion Europa silencio sismico, con sismos simulados en sitios aleatorios')
    case 3
     load ('Etiquetas_USA_Produccion_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion USA')
    case 4
     load ('Etiquetas_CHILE_Produccion_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion CHILE')
    case 5
     load ('Etiquetas_EUROPA_DURProduccion_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion DUR')
    case 6
     load ('Etiquetas_EUROPA_PEGProduccion_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion PEG')
    case 7
     load ('Etiquetas_PETProduccion_Ready.mat','Lx_c','Lx','Ly_c'...
      ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     % load ('Etiquetas_PETProduccion2_Ready.mat','Lx_c','Lx','Ly_c'...
     %  ,'Lon_S','Lat_S','Mag_S','Depthkm_S')  % DATOS DE PRODUCCION
     disp('Fase de Produccion PET')
   end
   disp('Fase de Produccion')
  end
  
  % VISUALIZACION
  % -----------------------------------------------------------------
  band = -2;    %  (0)  LA PRIMERA VEZ CALCULA Y SALVA LA DATA NORMALIZADA
  %  (1)  CARGA LA DATA NORMALIZADA, DA UN REPORTE BREVE, CALCULA Y VISUALIZA FIGURAS
  % (-1)  CARGA LA DATA NORMALIZADA CALCULA PERO NO VISUALIZA FIGURAS
  % (-2)  CARGA LA DATA NORMALIZADA CALCULA Y VISUALIZA FIGURAS SOBRE LA EFICIENCIA
  % =================================================================
  % CARGA CAMPO GEOMAGNETICO YA NORMALIZADO
  
  switch Data_Stations
   case 1 %JAPON:
    % ------------------------------------------------------
    load Campo_Geomagnetico_Normalizado  %Media 0 std(:)
    %load Campo_Geomagnetico_Normalizado_Japon.mat %Media 0 std 1
    % ------------------------------------------------------
   case 2 %SILENCIO EUROPA DOU_MAB_WNG
    % ------------------------------------------------------
    %load Campo_Geomagnetico_Normalizado_EUROPA_DOU.mat %Media 0 std 1
    % ------------------------------------------------------
    load Campo_Geomagnetico_Normalizado_EUROPA_DOU2.mat %Media 0 std (:)
   case 3 %USA:
    % ------------------------------------------------------
    % load('Campo_Geomagnetico_Normalizado_USA.mat')  %Maximo en 1 y distribucion por debajo
    load Campo_Geomagnetico_Normalizado_USA2.mat %Media 0 std(:)
    % load Campo_Geomagnetico_Normalizado_USA3.mat
    %load Campo_Geomagnetico_Normalizado_USA4 %Media 0 std 1
    % ------------------------------------------------------
   case 4 %CHILE:
    % ------------------------------------------------------
    % load Campo_Geomagnetico_Normalizado_CHILE  %Maximo en 1 y distribucion por debajo
    % load Campo_Geomagnetico_Normalizado_CHILE2  %Maximo en 1 y distribucion por debajo
    load Campo_Geomagnetico_Normalizado_CHILE3.mat  %Media 0 std(:)
    %load Campo_Geomagnetico_Normalizado_CHILE4 %Media 0 std 1
    % ------------------------------------------------------
   case 5
    % ------------------------------------------------------
    % load Campo_Geomagnetico_Normalizado_EUROPA_DUR %Maximo en 1 y distribucion por debajo
    load Campo_Geomagnetico_Normalizado_EUROPA_DUR2 %Media 0 std(:)
    %load Campo_Geomagnetico_Normalizado_EUROPA_DUR3 %Media 0 std 1
    % ------------------------------------------------------
   case 6
    % ------------------------------------------------------
    % load Campo_Geomagnetico_Normalizado_EUROPA_PEG  %Maximo en 1 y distribucion por debajo
    load Campo_Geomagnetico_Normalizado_EUROPA_PEG2  %Media 0 std(:)
    %load Campo_Geomagnetico_Normalizado_EUROPA_PEG3 %Media 0 std 1
    % ------------------------------------------------------
   case 7
    % ------------------------------------------------------
    % load Campo_Geomagnetico_Normalizado_PET.mat  %Media 0 std 1
    % ------------------------------------------------------
    load Campo_Geomagnetico_Normalizado_PET2.mat   %Media 0 std(:)
    
  end
  
  % DESCOMENTAR SECCION SI SE DESEA: %NORMALIZA A MEDIA 0 Y DESVIACION 1
  % -----------------------------------------------------------------
  % NORMALIZACION TIPO 2 SIGNIFICA DIVIDIR CADA DIA POR LA STD DEL DIA
  %              TIPO 1 SIGNIFICA DIVIDIR CADA DIA POR LA STD DEL CONJUNTO COMPLETO
  % -----------------------------------------------------------------
  % [Bxn,Byn,Bzn,F,Fo] = Proceso_Normaliza_media0_desv1(Label_Station,Data,1,0);
  % -----------------------------------------------------------------
  % [Bxn,Byn,Bzn,F,Fo] = Proceso_Normaliza_media0_desv1(Label_Station,Data,2,1);
  % -----------------------------------------------------------------
  % si se desea normalizar por todas las estaciones:(pero con todas las
  % estaciones en un vector)
  % [Bxt,Byt,Bzt]=NormalizaData2(Bxt,Byt,Bzt,1,band);
  % -----------------------------------------------------------------
  % save Campo_Geomagnetico_Normalizado_KNY.mat  Bxn Byn Bzn Fo  %Media 0 std 1
  % -----------------------------------------------------------------
  % save Campo_Geomagnetico_Normalizado_PET2.mat  Bxn Byn Bzn Fo  %Media 0 std(:)
  % -----------------------------------------------------------------
  % =================================================================
  % DEFINICION DE PARAMETROS DE ANALISIS
  % -----------------------------------------
  % PRECURSOR USADO
  % -----------------------------------------
  % (1)  ==> ENERGIA TOTAL
  % (2)  ==> PRESION DE RADIACION
  % (3)  ==> FLUJO DE MOMENTO (EN X ==> plano yz)
  % (4)  ==> FLUJO DE MOMENTO (EN Y ==> plano xz)
  % (5)  ==> FLUJO DE MOMENTO (EN Z ==> plano xy)
  % (6)  ==> TRAZA DEL TENSOR
  % (7)  ==> DENSIDAD DE MOMENTUM (VECTOR DE POINTING)
  % ---------------------------------------------------
  Componente_Tensor_Maxwell = 7;
  % ---------------------------------------------------
  % MOMENTO ESTADISTICO UTILIZADO PARA EL PRECURSOR
  % (1) MEDIA (2) KURTOSIS (3) DESVIACION STANDARD
  % ---------------------------------------------------
  Momento_Usado = 1;
  % ------------------------------------
  % UMBRAL CONSIDERADO PARA EL PRECURSOR
  % ------------------------------------
  uo = 1; % Umbral usado
  if band_display ==0
   display(['Umbral USADO:' num2str(uo) ' ']);
  end
  % -----------------------------------------------------------------
  % ESTACIONES A UTILIZAR: (1) KKA (2) KNY (3) MMB y (4) USA LAS 3 ESTACIONES SIMULTANEAMENTE.
  %                        (OJO) SE PUDEN USAR COMBINACIONES DE LAS
  % ESTACIONES ==> [1 2], [1 3], [2 3] o [1 2 3]
  % -----------------------------------------------------------------
  Station_Used_Now = [1 2 3 4];
  % -----------------------------------------------------------------
  
  % CALCULO DE LA EFICIENCIA PARA CADA ESTACION
  % -------------------------------------------
  
  for Actual_Station = 1:length(Station_Used_Now)
   
   Stadistic_Ef = [];% ESTADISTICA DE LA EFICIENCIA (para usar en ciclos de un prametro)
   Max_a = [];
   FP = [];
   
   if Actual_Station <=3
    Station_Used_Now = Actual_Station;  % SOLO UNA ESTACION POR LAZO
   else
    Station_Used_Now = [1 2 3];         % LAS TRES ESTACIONES SIUMULTANEAS
   end
   
   disp(['Station_Used_Now: ',num2str(Station_Used_Now)])
   
   for jg = 1:length(Gamma_2)
    if  band_display == 1
     disp([CASO Actual_Station Gamma_2(jg)])
    end
    % -----------------
    % VALORES DE GAMMA
    % -----------------
    Ng = 4;                                  % NUMEROS DE GAMMAS A GENERAR
    gamma_1 = Gamma_2(jg);                   % VALOR INFERIOR DE GAMMA
    gamma_2 = Gamma_2(jg) + Gamma_2(jg)/1e6; % VALOR SUPERIOR DE GAMMA
    G = linspace(gamma_1,gamma_2,Ng);        % VALORES DE GAMMA
    dg = 3;                                  % TAMAÑO DE LA MUESTRA DE GAMMA (GAMMAS USADOS EN LA ESTADISTICA) (debe ser >= 3)
    % ---------------------------------------------------------------------
    % NUMERO DE DIAS POR VENTANA
    % ---------------------------
    Md = 14; % NUMERO DE VENTANAS (SE CALCULAN AUTOMATICAMENTE A PARTIR DE Md Y TIEMPO TOTAL DE LOS DATOS)
    % ---------------------------------------------------------------------
    % ACTIVACION DEL PRECURSOR
    % --------------------------
    nd_b = 60/Md; %90/Md; % EL SISMO OCURRE (nd_b) ANTES DE LA ACTIVACION DEL PRECURSOR
    Nd_a = 60/Md; %6=120 dias%650/Md;%(50:20*5:95*Md)/Md; % EL SISMO OCURRE (nd_a) DESPUES DE LA ACTIVACION DEL PRECURSOR  ( (Nd_a) ETIQUETA DE ANCHO DE VENTANA)
    % ---------------------------------------------------------------------
    % INTERVALO HORARIO TOTAL
    % ------------------------------
    num_h = 1; % NUMERO DE HORAS DEL INTERVALO EN ESTUDIO
    Num_h = 1/num_h;
    dho = 1;
    Dh = [1:dho:24*Num_h];
    Dho = 2; % INTERVALO INCREMENTAL HORARIO (debe ser >=2)
    M_dh = [];
    for i=1:length(Dh) - Dho,dh = Dh(i):Dh(i+Dho); M_dh(i) = mean(dh); end % (M_dh) ETIQUETA DE INTERVALO HORARIO
    % ---------------------------------------------------------------------
    % CAMPOS A LLENAR
    % -----------------------
    IM = [];
    ST = {};
    k_ev = 0;
    Pt = [];
    Pa = [];
    AS_IND_C = [];
    AS_IND_F = [];
    AS_IND_CF = [];
    M_Gamma = [];
    AMPLITUD = {};
    IMAGE_HISTORY = {};
    % --------------------------
    %  Dr_C ==> DISTRIBUCION DE TIEMPOS PARA LOS EVENTOS NO CONTIGUOS
    % Ind_C ==> INDICES PARA LOS EVENTOS NO CONTIGOUS
    % CICLO PARA CADA EVENTO NO CONTIGUO SELECCIONADO
    % ===================================================
    for ev=1:1;%length(Dr_C)
     
     % PARA ETIQUETAR LOS EVENTOS CONTIGUOS
     % -----------------------------------------
     Lxx = Lx_c;
     Lyy = Ly_c;
     %------------------------------------------
     % SI USAMOS DATA SIMULADA
     % ------------------------
     %           Lxx = Lx_Sim;
     %           Lyy = Ly_Sim;
     %         Mag_S = Mag_Simulada;
     %            f1 = datenum('13-Jan-2004');
     %            f2 = datenum('01-Jan-2020 00:01:00');
     % datestr(Lxx(:,1))
     
     % PRESELECCION DE LOS EVENTOS (INTERVALO DE FECHAS REALES VALIDAS)
     % --------------------------------------------------------------------
     f1 = Lx(1,1);    % FECHA INFERIOR DE LOS DATOS TOTALES
     f2 = Lx(end,end);% FECHA SUPERIOR DE LOS DATOS TOTALES
     % INDICES VALIDOS
     % -----------------
     ind_s = find(Lxx(:,1)>=f1 & Lxx(:,1) <= f2); % INTERVALO DE FECHAS (f1 y f2)
     ind_sf = find(Fo(1,:)>=f1 & Fo(1,:) <= f2);
     Fo = Fo(:,ind_sf);
     [nf,mf,kf] = size(Fo);
     if length(ind_s) > 1
      k_ev = k_ev+1;
      nc = length(ind_s); % NUMERO TOTAL DE SISMOS EN ESTUDIO
      Lxx_u = Lxx(ind_s,:);
      Lyy_u = Lyy(ind_s,:);
      Ep = [];      % EFICIENCIA PORCENTUAL
      S = [];       % ESTATUS
      S1 = [];      % ESTATUS
      MAG_S = {};   % MAGNITUDES SELEC EN EL CICLO
      DEPTH_S = {}; % PROFUNDIDADES SELEC EN EL CICLO
      
      % -------------------------------------------------
      dga = 1;
      Ku = 0;
      Parametro = Nd_a;
      % =================================================
      INDICES = ones(nc,length(Nd_a));
      indices = 1:nc;
      for ku = 1:length(Parametro)
       ind_im = 0;
       for ga=1:dga:length(G)-dg
        Ku = Ku + 1;
        nd_a = Nd_a;        % Parametro(Ku);
        gamma = G(ga:ga+dg);% VALORES DE GAMMA
        % =========================================
        [gamma,E] = Energia_Vs_Gamma(gamma,0);  % VALORES DE LA ENERGIA EN GeV
        % =========================================
        % SALIDAS
        % ========
        Ef = [];     % EFICIENCIA ACIERTOS Y FALLOS
        IND_KC = []; % INDICE CIERTOS
        IND_KF = []; % INDICE FALLOS
        AMP_U = [];  % AMPLITUD ASOCIADA AL PRECURSOR EN CADA ACIERTO
        IM = [];
        count_surf = 0;
        % -----------------------------------------
        for i=1:length(Dh) - Dho
         count_surf = count_surf + 1;
         dh = Dh(i):Dh(i+Dho);
         m_dh = mean(dh/Num_h);
         % -------------------------------------
         Programa_General_Precursores_Relativistas % CALCULA EL PRECURSOR SELECCIONADO (FLUJO POR EJEMPLO)
         Encuentra_Sismos_Sicronizados % SISMOS LOCALIZADOS DIAS ANTES Y DIAS DESPUES DEL MAXIMO DEL PRECURSOR
         % -------------------------------------
         IND_KC = [IND_KC;IND_kc]; % SISMOS ENCONTRADOS CIERTOS
         IND_KF = [IND_KF;IND_kf]; % SISMOS NO ENCONTRADOS
         Ef = [Ef;m_dh pc pf pc+pf p_ind_c p_ind_f p_ind_c + p_ind_f CASO Actual_Station jg length(Gamma_2)];   % EFICIENCIA
         AMP_U = [AMP_U;Amp_U];
         
         % surf(P),shading interp,view([0 90]),pause
         
         Pt(:,:,count_surf) = P'; % EVOLUCION DEL PRECURSOR UMBRALIZADO (SURF)
         Pa(:,count_surf) = As'; % EVOLUCION DEL PRECURSOR SIN UMBRALIZAR (PROBABILIDAD CONDICIONADA)(SURF)
         
         % -------------------------------------
         % SOLO LOS INDICES CIERTOS Y FALLOS
         % -------------------------------------
         As_Ind_c = zeros(1,length(As));
         As_Ind_f = As_Ind_c;
         As_Ind_c(Ind_c) = As(Ind_c);
         As_Ind_f(IND_FP) = As(IND_FP);
         
         AS_IND_C(:,count_surf) = As_Ind_c;
         AS_IND_F(:,count_surf) = As_Ind_f;
         AS_IND_CF(:,count_surf) = As_Ind_c + As_Ind_f;
         % ========
         %  if size(AS_IND_C,2) > 1
         %    clf
         %    hg = 1;
         %    set(figure(1),'Position',[680 32 1197 906],'Color','W')
         %    subplot(4,1,1)
         %    surf(AS_IND_C'),shading interp,view([0 90]),title('AS_IND_C')
         %    xlim([0 Nv])
         %    subplot(4,1,2)
         %    plot((Lxx'-f1)/Md,Lyy','-r','LineWidth',[1]),grid,hold on
         %    plot(smooth(sum(AS_IND_C'),hg),'-b','LineWidth',[2]),
         %    hold off
         %    xlim([0 Nv])
         %    subplot(4,1,3)
         %    surf(AS_IND_F'),shading interp,view([0 90]),title('AS_IND_F')
         %    xlim([0 Nv])
         %    subplot(4,1,4)
         %    plot((Lxx'-f1)/Md,Lyy','-r','LineWidth',[1]),hold on
         %    plot(smooth(sum(AS_IND_F'),hg),'-r','LineWidth',[2])
         %    plot(smooth(sum(AS_IND_C'),hg),'--b','LineWidth',[2])
         %    grid,hold off  % SISMOS NO-ENCONTRADOS
         %    xlim([0 Nv])
         %    Ejes_Visibles(1)
         % %   set(figure(2),'Position',[680 32 1197 906],'Color','W')
         % %   plot((Lxx'-f1)/Md,Lyy','-r','LineWidth',[1]),hold on
         % %   plot(mean(AS_IND_CF').*std(AS_IND_CF')*20),grid,hold off
         %    pause(.01)
         % end
         % ========
         if band_display == 1
          % VISULIZACION DE LAS SALIDAS
          % ---------------------------------------------------------------
          fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------\n')
          disp('     m_dh       pc      pf        pc+pf     p_ind_c    p_ind_f  p_ind_c + p_ind_f      CASO    Actual_Station    jg    length(Gamma_2)')
          disp([Ef])
          fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------\n')
          % ---------------------------------------------------------------
         end
        end
        % =================================================================
        [mf,nf] = max(Ef(:,2));
        Stadistic_Ef = [Stadistic_Ef; mean(G) mean(E) mean(Ef(:,2)) std(Ef(:,2)) Ef(nf,1) Ef(nf,2)]; % ESTADISTICA DE LA EFICIENCIA
        % ETIQUETAS
        % --------------
        % (1) ==> VALOR MEDIO DE GAMMA
        % (2) ==> VALOR MEDIO DE LA ENERGIA
        % (3) ==> EFICIENCIA MEDIA EN TODO EL INTERVALO HORARIO DE INTERES (0 24 H)
        % (4) ==> DESVIACION STANDARD
        % (5) ==> INTERVALO HORARIO MEDIO
        % (6) ==> EFICIENCIA PORCENTUAL
        % -----------------------------------------
        [j_kc,i_kc] = hist(IND_KC,length(IND_KC));
        i_kc = round(i_kc);
        p_kc = (j_kc/length(IND_KC))*100; % PROBABILIDAD local
        IND_KCO = IND_KF;
        IND_KFO = IND_KC;
        % =========================================
        % VISUALIZACION DE EFICIENCIA
        % =========================================
        % ACIERTOS Y FALLOS VERSUS INTERVALO HORARIO
        % -----------------------------------------
        Int_H = Ef(:,1);
        Aciertos = Ef(:,2);
        Fallos = Ef(:,3);
        Total_a = Ef(:,4);
        [jMc,iMc] = max(Aciertos);
        % -----------------------------------------
        Ind_Aciertos = Ef(:,5);
        Ind_Fallos = Ef(:,6);
        Total_f = Ef(:,7);
        [jMf,iMf] = max(Ind_Fallos);
        % -----------------------------------------
        % IND FALSOS POSITIVOS VS INTERVALO HORARIO
        % -----------------------------------------
        % SISMOS NO ENCONTRADOS
        % ELIMINA INDICES REPETIDOS
        % =========================================
        IND_FO = [];
        IND_FC = [];
        while length(IND_KFO)>0
         Ind_R = IND_KFO(1);
         ind_e = find(Ind_R==IND_KFO);
         if length(ind_e)>1
          IND_FO = [IND_FO;IND_KFO(1)];
          IND_KFO(ind_e) = [];
         else
          IND_KFO(1) = [];
         end
        end
        Ind_S = ind_s;
        Ind_S(IND_FO) = [];
        IND_FO = sort(IND_FO);
        % -----------------------------------------
        INDICES(IND_FO,ku) = IND_FO;
        Ind_New_In = [];
        Ind_New_Out = [];
        if ku>1
         Ind_New_In = find(INDICES(:,ku-1)-INDICES(:,ku)<0);
         %Ind_New_Out = find(INDICES(:,ku-1)-INDICES(:,ku)>0);
        else
         Ind_New_In = IND_FO;
         %Ind_New_Out = find(INDICES(:,ku)>0);
        end
        if length(Ind_New_In)>0
         Ind_New_Events = Ind_New_In;
         % -------------------------------------
         St = [1 2 3];     % ESTACIONES
         %
         % ESTIMACION DE LA DISTANCIA A ESTACIONES USADAS Y NUEVOS EVENTOS SEGUN EL PARAMETRO DE CONTROL USADO
         % -------------------------------------
         dx_station = y_station(Station_Used_Now)/grados_km;
         dy_station = x_station(Station_Used_Now)/grados_km;
         Dx_Station_m = mean(dx_station);
         Dy_Station_m = mean(dy_station);
         % -------------------------------------
         % DISTANCIA ENTRE EL CENTRO GEOMETRICO ENTRE ESTACIONES USADAS Y EVENTOS
         % -------------------------------------
         if length(IND_FO)>0
          R_Dx_Station = repmat(Dx_Station_m,length(Ind_New_Events),1);
          R_Dy_Station = repmat(Dy_Station_m,length(Ind_New_Events),1);
          Rx_St_Event = sqrt((R_Dx_Station - Lon_S(Ind_New_Events)).^2 + (R_Dy_Station - Lat_S(Ind_New_Events)).^2);
          Rm_St_Event = (Rx_St_Event)*grados_km; % DISTANCIA MEDIA NUEVOS EVENTOS RESPECTO AL CENTRO GEOMETRICO DE LAS ESTACIONES USADAS
         else
          Rm_St_Event = 0;
         end
         
         % -------------------------------------
         % SALIDA: EFICIENCIA PORCENTUAL EVENTOS CIERTOS ENCONTRADOS
         % ----------------------------------------------------------------
         Ep(Ku) = 100 - (length(Ind_S)/length(ind_s))*100;
         MAG_S(Ku).MAG_S = Mag_S(Ind_New_Events);
         DEPTH_S(Ku).DEPTH_S = Depthkm_S(Ind_New_Events);
         AMPLITUD(Ku).AMPLITUD = AMP_U;
         IMAGE_HISTORY(Ku).IMAGE_HISTORY = IM;
         % -------------------------------------
         %            (1)    (2)    (3)       (4)            (5)            (6)               (7)             (8)       (9)            (10)       (11)           (12)             (13)               (14)               (15)                   (16)              (17)                (18)
         S = [S;  nd_a mean(E) Ep(Ku) mean(Aciertos) mean(Fallos) mean(Ind_Aciertos) mean(Ind_Fallos) Int_H(iMc) Aciertos(iMc) Fallos(iMc) Int_H(iMf) Ind_Aciertos(iMf) Ind_Fallos(iMf) mean(Mag_S(Ind_New_Events)) mean(Lat_S(Ind_New_Events)) mean(Lon_S(Ind_New_Events)) mean(Depthkm_S(Ind_New_Events)) mean(Rm_St_Event)];
         S1 = [S1;Depthkm_S(Ind_New_Events) Rm_St_Event];
         % LEYENDA:
         % --------
         %  (1) ==> VALOR MEDIO DE GAMMA (o del umbral uo)
         %  (2) ==> VALOR MEDIO DE LA ENERGIA EN GeV
         %  (3) ==> EFICIENCIA PORCENTUAL PARA CADA VALOR DE GAMMA
         %  (4) ==> VALOR MEDIO DE LOS SISMOS ACERTADOS
         %  (5) ==> VALOR MEDIO DE LOS SISMOS FALLIDOS
         %  (6) ==> VALOR MEDIO DE LOS INDICES ACERTADOS (CASOS POSITIVOS)
         %  (7) ==> VALOR MEDIO DE LOS INDICES FALLIDOS (FALSOS POSITIVOS)
         %  (8) ==> INTERVALO HORARIO DONDE EL MAXIMO DE ACIERTOS
         %  (9) ==> VALOR MAXIMO DE LOS ACIERTOS (EN INTERVALO 0 A 24 HORAS)
         % (10) ==> VALOR MAXIMO DE LOS FALLOS (EN INTERVALO 0 A 24 HORAS)
         % (11) ==> INTERVALO HORARIO DONDE OCURRE EL MAXIMO DE LOS INDICES ACERTADOS
         % (12) ==> VALOR MAXIMO DE LOS INDICES ACERTADOS
         % (13) ==> VALOR MAXIMO DE LOS INDICES FALLIDOS
         % (14) ==> MEDIA DE LA MAGNITUD PARA LOS EVENTOS ENCONTRADOS
         % (15) ==> MEDIA DE LA LATITUD PARA LOS EVENTOS ENCONTRADOS
         % (16) ==> MEDIA DE LA LONGITUD PARA LOS EVENTOS ENCONTRADOS
         % (17) ==> MEDIA DE LA PROFUNDIDAD SISMICA PARA LOS EVENTOS ENCONTRADOS
         % (18) ==> DISTANCIA PROMEDIO ENTRE ESTACIONES USADAS Y NUEVOS EVENTOS
         % ----------------------------------------------------------------
         
         if band == -3 | band == -1
          set(figure(5),'Position',[5 514 1912 464],'Color','W')
          plot(Int_H,Aciertos,'.-r',Int_H,Aciertos,'or',Int_H,Fallos,'.-b',Int_H,Fallos,'ob',Int_H,Total_a,'.-m',Int_H,Total_a,'om','LineWidth',[2]),hold on
          plot(Int_H(iMc),Aciertos(iMc),'+k','MarkerSize',[12],'LineWidth',[2]),hold off
          xlabel('Time interval from 0 to 24 h')
          ylabel('Eficiency in (%)')
          title(['Efmax: [' num2str(Aciertos(iMc)) '(%)] --- Dhm: [' num2str(Int_H(iMc)) '] hours' ])
          xlim([Dh(1) Dh(end)/Num_h])
          ylim([0 100])
          grid
          Ejes_Visibles(5)
          g = get(5);
          g1 = g.Children;
          set(g1,'Position',[0.0415    0.1250    0.9514    0.8103])
          legend('Thrue','','False','','Total','Location','Best')
          ax1 = gca;        % current axes
          ax1.XTick = [Dh(1):1:Dh(end)/Num_h];
          
          
          set(figure(6),'Position',[5 514 1912 464],'Color','W')
          plot(Int_H,Ind_Aciertos,'.-r',Int_H,Ind_Aciertos,'or',Int_H,Ind_Fallos,'.-b',Int_H,Ind_Fallos,'ob',Int_H,Total_f,'.-m',Int_H,Total_f,'om','LineWidth',[2]),hold on
          plot(Int_H(iMf),Ind_Fallos(iMf),'+k','MarkerSize',[12],'LineWidth',[2]),hold off
          xlabel('Time interval from 0 to 24 h')
          ylabel('Eficiency (Indices) in (%)')
          title(['Indices falses positives: [' num2str(Ind_Fallos(iMf)) '(%)] --- Dhm: [' num2str(Int_H(iMf)) '] hours' ])
          xlim([Dh(1) Dh(end)/Num_h])
          ylim([0 100])
          grid
          g = get(6);
          g1 = g.Children;
          set(g1,'Position',[0.0415    0.1250    0.9514    0.8103])
          Ejes_Visibles(6)
          legend('Thrue','','False','','Total','Location','Best')
          ax1 = gca;        % current axes
          ax1.XTick = [Dh(1):1:Dh(end)/Num_h];
          % STATIONS VISUALIZATION
          % ----------------------
          set(figure(7),'Position',[697 32 1220 946],'Color','w')
          plot(Lon_S(ind_s),Lat_S(ind_s),'og',Lon_S(ind_s),Lat_S(ind_s),'+r'),hold on
          plot(Lon_S(IND_FO),Lat_S(IND_FO),'oc',Lon_S(IND_FO),Lat_S(IND_FO),'+k')
          plot(Lon_S(Ind_New_Events),Lat_S(Ind_New_Events),'*m',Lon_S(Ind_New_Events),Lat_S(Ind_New_Events),'sk','LineWidth',[2],'MarkerSize',[8])
          grid
          for i=1:length(Station_Used_Now)
           text(x_station(St(Station_Used_Now(i)))/grados_km,y_station(St(Station_Used_Now(i)))/grados_km,Label_Station(St(Station_Used_Now(i))),'Color','b')
           text(x_station(St(Station_Used_Now(i)))/grados_km-0.2,y_station(St(Station_Used_Now(i)))/grados_km,'o','Color','r')
          end
          for i=1:length(St)
           text(y_station(i)/grados_km,x_station(i)/grados_km,Label_Station(i),'Color','b')
          end
          for i=1:length(Ind_New_Events)
           text(Lon_S(Ind_New_Events(i))+0.2,Lat_S(Ind_New_Events(i)),num2str(Mag_S(Ind_New_Events(i))),'FontSize',[12])
          end
          plot(Dx_Station_m+.025,Dy_Station_m,'ob',Dx_Station_m+.025,Dy_Station_m,'+r','MarkerSize',[10],'LineWidth',[1])
          title(['Total events: [ ' num2str(nc) ' ] Old events: [ ' num2str(length(IND_FO)) ' ] New events: [ ' num2str(length(Ind_New_Events)) ' ]  Mag: [ ' num2str(mean(Mag_S(Ind_New_Events))) ' ]'])
          ylabel('Latitude')
          xlabel('Longitude')
          hold off
          ylim([20 48])
          xlim([120 150])
          set(gca,'dataaspectratiomode','manual'),set(gca,'dataaspectratio',[1 1 1]);
          Ejes_Visibles(7)
          % --------------------------
          set(figure(8),'Position',[697 32 1220 946],'Color','w')
          plot(Lxx_u'-f1,Lyy_u','-r',Lxx_u'-f1,Lyy_u','--g'),hold on,plot(Lxx_u(IND_FO,:)'-f1,Lyy_u(IND_FO,:)','--r',Lxx_u(IND_FO,:)'-f1,Lyy_u(IND_FO,:)','-c','LineWidth',[2])
          plot(Lxx_u(Ind_New_Events,:)'-f1,Lyy_u(Ind_New_Events,:)','-k',Lxx_u(Ind_New_Events,:)'-f1,Lyy_u(Ind_New_Events,:)','--y','LineWidth',[2])
          grid
          title(['Total events: [ ' num2str(nc) ' ] Old events: [ ' num2str(length(IND_FO)) ' ] New events: [ ' num2str(length(Ind_New_Events)) ' ]  Mag: [ ' num2str(mean(Mag_S(Ind_New_Events))) ' ]'])
          Ejes_Visibles(8)
          
          set(figure(9),'Position',[826 32 1091 946],'Color','W')
          subplot(3,1,1)
          plot(S(:,1),S(:,3),'.-r',S(:,1),S(:,3),'ob')
          grid
          title('Total eficiency')
          ylabel('Eficiency (%)')
          ylim([0 100])
          subplot(3,1,2)
          plot(S(:,1),S(:,4),'.-r',S(:,1),S(:,4),'og'),hold on
          plot(S(:,1),S(:,5),'.-b',S(:,1),S(:,5),'ob'),hold off
          grid
          title('Seismic events as:  hits (red) and falses (blue)')
          ylabel('Eficiency (%)')
          ax1 = gca;  % current axes
          %ax1.XTick = [0:0.5:S(end,1)]
          %ax1.YTick = [0:10:100];
          subplot(3,1,3)
          plot(S(:,1),S(:,6),'.-r',S(:,1),S(:,6),'og'),hold on
          plot(S(:,1),S(:,7),'.-b',S(:,1),S(:,7),'ob'),hold off
          grid
          title('Indices as: hits (red) and falses (blue)')
          ylabel('Eficiency (%)')
          xlabel('Parameter')
          ax1 = gca;  % current axes
          %     ax1.XTick = [0:0.5:S(end,1)]
          %     ax1.YTick = [0:10:100];
          Ejes_Visibles(9)
          
          set(figure(10),'Position',[826 32 1091 946],'Color','W')
          subplot(2,1,1)
          plot(S(:,1),S(:,9),'.-b',S(:,1),S(:,9),'ob'),hold on,plot(S(:,1),S(:,10),'.-r',S(:,1),S(:,10),'or'),hold off,grid
          title('Maximo de Aciertos (b) y fallos (red) (eventos)')
          subplot(2,1,2)
          plot(S(:,1),S(:,12),'.-b',S(:,1),S(:,12),'ob'),hold on,plot(S(:,1),S(:,13),'.-r',S(:,1),S(:,13),'or'),hold off,grid
          title('Maximo de Aciertos (b) y fallos (red) (indices)')
          Ejes_Visibles(10)
          
          set(figure(11),'Position',[826 32 1091 946],'Color','W')
          subplot(3,1,1)
          plot(S(:,1),S(:,14),'.-b',S(:,1),S(:,14),'ob'),grid
          title('Media de las magnitudes encontradas')
          subplot(3,1,2)
          plot(S(:,1),S(:,17),'.-b',S(:,1),S(:,17),'ob'),grid
          title('Media de las profundidades encontradas')
          ylabel('Km')
          subplot(3,1,3)
          plot(S(:,1),S(:,18),'.-b',S(:,1),S(:,18),'ob'),grid
          title('Media de la distancia estacion y eventos nuevos')
          ylabel('Km')
          xlabel('Parameter')
          Ejes_Visibles(11)
          
          set(figure(12),'Position',[826 32 1091 946],'Color','W')
          if size(INDICES,2) >1
           surf(INDICES),colorbar
           title('History events')
           xlabel('Parameter')
           ylabel('Events number')
           zlabel('Events number')
           %    ylim([1 nc])
           %    xlim([1 length(Uoo)])
           % zlim([1 nc])
           grid
           view([0 90])
           cb = colorbar;
           set(cb,'Ticks',[0:5:100])
           Ejes_Visibles(12)
          else
           plot(INDICES,'.--r','LineWidth',[2])
           title('History events')
           xlabel('Parameter')
           ylabel('Events number')
           %    ylim([1 nc])
           %    xlim([1 length(Uoo)])
           % zlim([1 nc])
           grid
           Ejes_Visibles(12)
          end
          pause(.1)
         end
        end
       end
      end
      ST(k_ev).ST = S;
     end
    end
    
    if Caso_Ent_Prod == 0
     Pa = AS_IND_C;
     caso_salida=1;
     Stadistic_Results{Actual_Station,jg,caso_salida} =  struct('Stadistic_Ef',Stadistic_Ef,'Gamma_2',Gamma_2,'Md',Md,'M_dh',M_dh,'Nd_a',Nd_a,'num_h',num_h,'dg',dg,'Ng',Ng,'uo',uo,...
      'Momento_Usado',Momento_Usado,'Total_a',Total_a,'Componente_Tensor_Maxwell',Componente_Tensor_Maxwell,...
      'Station_Used_Now',Station_Used_Now,'iMc',iMc,'Int_H',Int_H,'Aciertos',Aciertos,'Fallos',Fallos,'t',t,...
      'Label_Station',Label_Station,'dho',dho,'Dh',Dh,'nh',nh,'Num_h',Num_h,'f1',f1,'f2',f2,'Nv',Nv,...
      'M_Gamma',M_Gamma,'Max_a',Max_a,'Caso_Ent_Prod',Caso_Ent_Prod,'Pt',Pt,'Pa',Pa);
     %save STADISTIC_ENTRENA_R_Gamma_100_10_1 Stadistic_Entrena
    else
     caso_salida=2;
     Stadistic_Results{Actual_Station,jg,caso_salida} =  struct('Stadistic_Ef',Stadistic_Ef,'Gamma_2',Gamma_2,'Md',Md,'M_dh',M_dh,'Nd_a',Nd_a,'num_h',num_h,'dg',dg,'Ng',Ng,'uo',uo,...
      'Momento_Usado',Momento_Usado,'Total_a',Total_a,'Componente_Tensor_Maxwell',Componente_Tensor_Maxwell,...
      'Station_Used_Now',Station_Used_Now,'iMc',iMc,'Int_H',Int_H,'Aciertos',Aciertos,'Fallos',Fallos,'t',t,...
      'Label_Station',Label_Station,'dho',dho,'Dh',Dh,'nh',nh,'Num_h',Num_h,'f1',f1,'f2',f2,'Nv',Nv,...
      'M_Gamma',M_Gamma,'Max_a',Max_a,'Caso_Ent_Prod',Caso_Ent_Prod,'Pt',Pt,'Pa',Pa);
     % save STADISTIC_PRODUCE_R_Gamma_100_10_1 Stadistic_Produce
    end
    
   end % FIN DE CASO
  end
 end
 
 % SALIDA EN CADA LAZO
 % -------------------------------------------------------------------
 switch Data_Stations
  % std(:); == TYPE 1 NORMALIZACION EN LOS DATOS.
  % std 1;== TYPE 2 NORMALIZACION EN LOS DATOS.
  % STADISTIC_ENTRENA/PRODUCE_ZONA_VERSION % METADATOS
  
  case 1 %JAPON
   % ========= VIEJO FORMATO =====================================
   % save STADISTIC_ENTRENA_R_Japon_100_10_1 Stadistic_Entrena  %std(:);
   % save STADISTIC_PRODUCE_R_Japon_100_10_1 Stadistic_Produce  %std(:);
   % save STADISTIC_ENTRENA_R_Japon_100_10_1b Stadistic_Entrena %ho=1;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);
   % save STADISTIC_PRODUCE_R_Japon_100_10_1b Stadistic_Produce %ho=1;Media 0 std 1 uo=1;Gamma_2 = linspace(1,100,20);
   % save STADISTIC_ENTRENA_R_Japon_100_10_1c Stadistic_Entrena %Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,40);Nd_a = 0
   % save STADISTIC_PRODUCE_R_Japon_100_10_1c Stadistic_Produce %Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,40);Nd_a = 0
   % save STADISTIC_ENTRENA_R_Japon_100_10_1d Stadistic_Entrena %Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,50);Nd_a = 0;nd_b = 30/Md
   % save STADISTIC_PRODUCE_R_Japon_100_10_1d Stadistic_Produce %Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,50);Nd_a = 0;nd_b = 30/Md
   % save STADISTIC_ENTRENA_R_Japon_100_10_1f Stadistic_Entrena %Media 0 std 1;uo=1;Gamma_2 = linspace(25,29,10);Nd_a = 30/Md;nd_b = 30/Md
   % save STADISTIC_PRODUCE_R_Japon_100_10_1f Stadistic_Produce %Media 0 std 1;uo=1;Gamma_2 = linspace(25,29,10);Nd_a = 30/Md;nd_b = 30/Md
   % ========= NUEVO FORMATO =====================================
   % save STADISTIC_Results_Japon  Stadistic_Results %std(:)
   % save STADISTIC_Results_Japona Stadistic_Results %NO SE GUARDO%ho=1;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 30/Md;nd_b = 30/Md
   % save STADISTIC_Results_Japonb Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);Nd_a = 60/Md;nd_b = 60/Md
   % save STADISTIC_Results_Japonc Stadistic_Results %ho=1;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);Nd_a = 60/Md;nd_b = 60/Md
   % save STADISTIC_Results_Japond Stadistic_Results %ho=1;Media 0 std (:);uo=1;Gamma_2 = (1,80,10);Nd_a = 90/Md;nd_b = 90/Md num_h = 1/2
   % save STADISTIC_Results_Japonf Stadistic_Results %ho=1;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);
   % save STADISTIC_Results_Japong Stadistic_Results %Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,40);Nd_a = 0
   % save STADISTIC_Results_Japonh Stadistic_Results %Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,50);Nd_a = 0;nd_b = 30/Md
   % save STADISTIC_Results_Japoni Stadistic_Results %Media 0 std 1;uo=1;Gamma_2 = linspace(25,29,10);Nd_a = 30/Md;nd_b = 30/Md
   % save STADISTIC_Results_Japonz	Stadistic_Results %std(:)	60	60	1 h	(1,100,10)	1	10
   % save STADISTIC_Results_Japonk	Stadistic_Results %std(:)	90	90	1 h	(1,100,10)	1	10
   % save STADISTIC_Results_Japonj	Stadistic_Results %std(:)	90	90	30 m	(1,100,10)	1	10
   % save('STADISTIC_Results_Japonm', Stadistic_Results) %std(:)	90	90	1 h	(1,100,100)	1	10
   % save STADISTIC_Results_Japonp	%std(:)	60	60	1 h	(1,100,40)	1	10
   
  case 2 %EUROPA SILENCIO SISMICO CON SISMOS SIMULADOS DUO
   %save STADISTIC_Results_DOU Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,80,10);Nd_a = 60/Md;nd_b = 60/Md
   %=======================================================================
   %save STADISTIC_Results_DOUa Stadistic_Results %ho=10;Media 0
   %std 1;uo=1;Gamma_2 = linspace(1,80,5);Nd_a = 60/Md;nd_b = 60/Md Entrenamiento2
   % ======================================================================
   %save STADISTIC_Results_DOUb Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 =linspace(1,80,5);Nd_a = 60/Md;nd_b = 60/Md Entrenamiento3
   %save STADISTIC_Results_DOUc Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 = linspace(1,80,5);Nd_a = 60/Md;nd_b = 60/Md Entrenamiento3
   %save STADISTIC_Results_DOUd Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 =linspace(1,10,5);Nd_a = 60/Md;nd_b = 60/Md Entrenamiento3
   %save STADISTIC_Results_DOUd Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 =linspace(1,20,10);Nd_a = 60/Md;nd_b = 60/Md Entrenamiento4 num_h = 1/2;
   %save STADISTIC_Results_DOUf Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 =linspace(1,80,10);Nd_a = 90/Md;nd_b = 90/Md Entrenamiento4 num_h = 1/2;
   %save STADISTIC_Results_DOUg	Stadistic_Results %std (:)	90	90	1 h	(1,100,35)		10	Entrenamiento2 y producion 4
   % save STADISTIC_Results_DOUh	Stadistic_Results %std (:)	60 60	30m	(1,50,10)		10	Entrenamiento4 y producion 4
  case 3 %USA
   % ========= VIEJO FORMATO =====================================
   % save STADISTIC_ENTRENA_R_USA_100_10_1 Stadistic_Entrena % std(:);
   % save STADISTIC_PRODUCE_R_USA_100_10_1 Stadistic_Produce % std(:);
   % save STADISTIC_ENTRENA_R_USA_100_10_1b Stadistic_Entrena % std(:); UO=1 gamma del 1 al 100, 50 valores
   % save STADISTIC_PRODUCE_R_USA_100_10_1b Stadistic_Produce %  std(:);UO=1 gamma del 1 al 100, 50 valores
   % ========= NUEVO FORMATO =====================================
   % save STADISTIC_Results_USA Stadistic_Results % std(:);
   %         save STADISTIC_Results_USAa      Stadistic_Results %ho=1;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 30/Md;nd_b = 30/Md
   %          save STADISTIC_Results_USAb      Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);Nd_a = 90/Md;nd_b = 90/Md
   % ====================================================
   % save STADISTIC_Results_USAc     Stadistic_Results
   %          %ho=10;Media 0 std(:);uo=1;Gamma_2 =
   %          linspace(1,100,10);Nd_a = 90/Md;nd_b = 90/Md num_h = 1/2
   % ====================================================
   % save STADISTIC_Results_USAd     Stadistic_Results
   %          %ho=10;Media 0 std(:);uo=1;Gamma_2 =
   %          linspace(10,15,5);Nd_a = 90/Md;nd_b = 90/Md num_h = 1/2
   % ====================================================
   % save STADISTIC_Results_USAf     Stadistic_Results
   %          %ho=10;Media 0 std(:);uo=1;Gamma_2 =
   %          linspace(25,100,20);Nd_a = 60/Md;nd_b = 60/Md num_h = 1/2
   % ====================================================
   % save STADISTIC_Results_USAg    Stadistic_Results
   %          %ho=10;Media 0 std(:);uo=1;Gamma_2 =
   %          linspace(25,100,20);Nd_a = 60/Md;nd_b = 60/Md num_h = 1
   % ====================================================
   % save STADISTIC_Results_USAh Stadistic_Results %std(:);UO=1 gamma del 1 al 100, 50 valores
  case 4
   % ========= VIEJO FORMATO =====================================
   % save STADISTIC_ENTRENA_R_CHILE_100_10_1a Stadistic_Entrena % std(:);%UO=0.5
   % save STADISTIC_PRODUCE_R_CHILE_100_10_1a Stadistic_Produce % std(:);%UO=0.5
   % save STADISTIC_ENTRENA_R_CHILE_100_10_1b Stadistic_Entrena % std(:);%UO=1 gamma del 1 al 10, 100 valores
   % save STADISTIC_PRODUCE_R_CHILE_100_10_1b Stadistic_Produce % std(:);% UO=1 gamma del 1 al 10, 100 valores
   % save STADISTIC_ENTRENA_R_CHILE_100_10_1c Stadistic_Entrena % std(:);% UO=1 gamma del 1 al 100, 50 valores
   % save STADISTIC_PRODUCE_R_CHILE_100_10_1c Stadistic_Produce % std(:);%UO=1 gamma del 1 al 100, 50 valores
   % save STADISTIC_ENTRENA_R_CHILE_100_10_1d Stadistic_Entrena % std(:);% UO=0.5 gamma del 1 al 50, 10 valores
   % save STADISTIC_PRODUCE_R_CHILE_100_10_1d Stadistic_Produce % std(:);%UO=0.5 gamma del 1 al 50, 10 valores
   % ========= NUEVO FORMATO =====================================
   % save STADISTIC_Results_Chile  Stadistic_Results % std(:);%UO=0.5
   % save STADISTIC_Results_Chilea Stadistic_Results %ho=1;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 30/Md;nd_b = 30/Md
   % save STADISTIC_Results_Chileb Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);Nd_a = 60/Md;nd_b = 60/Md
   % save STADISTIC_Results_Chilec Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 = linspace(1,80,10);Nd_a = 90/Md;nd_b = 90/Md num_h = 1/2
   % save STADISTIC_Results_Chiled Stadistic_Results %std(:) UO=1 gamma del 1 al 10, 100 valores
   % save STADISTIC_Results_Chilef Stadistic_Results %std(:) UO=1 gamma del 1 al 100, 50 valores
   % save STADISTIC_Results_Chileg Stadistic_Results % std(:) UO=0.5 gamma del 1 al 50, 10 valores
   % save STADISTIC_Results_Chileh %	std(:)	60	60	1 h	(1,100,35)	0,5	10
   
  case 5
   % ========= VIEJO FORMATO =====================================
   % save STADISTIC_ENTRENA_R_DUR2_100_10_1 Stadistic_Entrena % std(:);%UO = 1; SIN SE?AL REGISTRADA:
   % save STADISTIC_PRODUCE_R_DUR2_100_10_1 Stadistic_Produce % std(:);%UO = 1; SIN SE?AL REGISTRADA:
   % save STADISTIC_ENTRENA_R_DUR2_100_10_1b Stadistic_Entrena % std(:);%UO = 0.6; gamma 1 al 100, 10 valores SIN SE?AL REGISTRADA:
   % save STADISTIC_PRODUCE_R_DUR2_100_10_1b Stadistic_Produce % std(:);%UO = 0.6; gamma 1 al 100, 10 valores SIN SE?AL REGISTRADA:
   % save STADISTIC_ENTRENA_R_DUR2_100_10_1c Stadistic_Entrena % std(:);%UO = 0.5; gamma 1 al 100, 10 valores; SIN SE?AL REGISTRADA:
   % save STADISTIC_PRODUCE_R_DUR2_100_10_1c Stadistic_Produce % std(:);%UO = 0.5; gamma 1 al 100, 10 valores; SIN SE?AL REGISTRADA:
   % save STADISTIC_ENTRENA_R_DUR2_100_10_1d Stadistic_Entrena % std(:);%UO = 0.3; gamma 1 al 100, 10 valores;
   % save STADISTIC_PRODUCE_R_DUR2_100_10_1d Stadistic_Produce % std(:);%UO = 0.3; gamma 1 al 100, 10 valores;
   % save STADISTIC_ENTRENA_R_DUR2_100_10_1f Stadistic_Entrena % std(:);%UO = 0.3; gamma 1 al 100, 100 valores;
   % save STADISTIC_PRODUCE_R_DUR2_100_10_1f Stadistic_Produce % std(:);%UO = 0.3; gamma 1 al 100, 100 valores;
   % ========= NUEVO FORMATO =====================================
   % save STADISTIC_Results_DUR  Stadistic_Results % std(:);%UO = 1; SIN SE?AL REGISTRADA:
   % save STADISTIC_Results_DURa Stadistic_Results %ho=1;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 30/Md;nd_b = 30/Md
   % save STADISTIC_Results_DURb Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);Nd_a = 60/Md;nd_b = 60/Md
   % save STADISTIC_Results_DURc Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 = linspace(1,80,10);Nd_a = 90/Md;nd_b = 90/Md num_h = 1/2
   % save STADISTIC_Results_DURd Stadistic_Results % std(:);%UO = 0.6; gamma 1 al 100, 10 valores SIN SE?AL REGISTRADA:
   % save STADISTIC_Results_DURf Stadistic_Results % std(:);%UO = 0.5; gamma 1 al 100, 10 valores; SIN SE?AL REGISTRADA:
   % save STADISTIC_Results_DURg Stadistic_Results % std(:);%UO = 0.3; gamma 1 al 100, 10 valores;
   % save STADISTIC_Results_DURh Stadistic_Results % std(:);%UO = 0.3; gamma 1 al 100, 100 valores;
  case 6
   % ========= VIEJO FORMATO =====================================
   % save STADISTIC_ENTRENA_R_PEG2_100_10_1 Stadistic_Entrena  % std(:);
   % save STADISTIC_PRODUCE_R_PEG2_100_10_1 Stadistic_Produce  % std(:);
   % save STADISTIC_ENTRENA_R_PEG2_100_10_1b Stadistic_Entrena % std(:);
   % save STADISTIC_PRODUCE_R_PEG2_100_10_1b Stadistic_Produce % std(:);
   % ========= NUEVO FORMATO =====================================
   % save STADISTIC_Results_PEG Stadistic_Results  % std(:);
   % save STADISTIC_Results_PEGa Stadistic_Results %ho=1;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 30/Md;nd_b = 30/Md
   % save STADISTIC_Results_PEGb Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);Nd_a = 60/Md;nd_b = 60/Md
   % save STADISTIC_Results_PEGc Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 = linspace(1,80,10);Nd_a = 90/Md;nd_b = 90/Md num_h = 1/2
   % save STADISTIC_Results_PEGd Stadistic_Results % std(:);
  case 7
   % ========= NUEVO FORMATO =====================================
   % save STADISTIC_Results_PET  Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);Nd_a = 60/Md;nd_b = 60/Md
   % save STADISTIC_Results_PETa Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,20);Nd_a = 30/Md;nd_b = 30/Md
   % save STADISTIC_Results_PETb Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 60/Md;nd_b = 60/Md
   % save STADISTIC_Results_PETc Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 = linspace(1,80,10);Nd_a = 90/Md;nd_b = 90/Md num_h = 1/2
   % save STADISTIC_Results_PETd Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 = linspace(1,80,10);Nd_a = 60/Md;nd_b = 60/Md %produccion 2
   % save STADISTIC_Results_PETf Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 60/Md;nd_b = 60/Md %produccion 2
   % save STADISTIC_Results_PETg Stadistic_Results %ho=10;Media 0 std (:);uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 60/Md;nd_b = 60/Md
   % save STADISTIC_Results_PETh Stadistic_Results %ho=10;Media 0 std 1;uo=1;Gamma_2 = linspace(1,100,10);Nd_a = 60/Md;nd_b = 60/Md  num_h = 1/2
   % save STADISTIC_Results_PETi	Stadistic_Results %std (:)	90	90	1	(1,100,50)	1	10
   % save STADISTIC_Results_PETj	Stadistic_Results %std (:)	90	0	1	(1,100,20)	1	10	Entrenamiento 2 y produccion 2
   % save STADISTIC_Results_PETk	Stadistic_Results %std (:)	60	60	1	(1,100,40)	1	10
   
 end
end
% FIN
% =========================================================================
