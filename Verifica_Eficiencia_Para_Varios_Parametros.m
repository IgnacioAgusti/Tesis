 % /home/aramis/Escritorio/Programas_News_Ignacio/ProgramasJapon
% Verifica_Eficiencia_Para_Varios_Parametros.m
% 
% VERIFICA EFICIENCIA DE LOS PRECURSORES RELATIVISTAS
% =====================================================

   clear all
   close all
% ====================================================================================================================================================

   Stadistic_Entrena  = {};
    Stadistic_Produce = {}; 
    
for CASO = 0:1                                  % (0) PARA ENTRENAMIENTO (1) PARA PRODUCCION 

   Caso_Ent_Prod = CASO;


% INTERVALOS HORARIOS DE MAYOR EFICIENCIA EN LA LOCALIZACION DE EVENTOS
% =================================================================================================
 if Caso_Ent_Prod == 0
       load Entrena_Etiqueta_Eventos_Japan_II_6M_2000_2019  % DATOS DE ENTRENAMIENTO
       disp('Fase de Entrenamiento')
 else 
       load Produce_Etiqueta_Eventos_Japan_II_6M_2000_2019  % DATOS DE PRODUCCION
       disp('Fase de Produccion')
 end
 
% VISUALIZACION
% -------------------------------------------------------------------------------------------------
  band = -1;    %  (0)  LA PRIMERA VEZ CALCULA Y SALVA LA DATA NORMALIZADA
                %  (1)  CARGA LA DATA NORMALIZADA, DA UN REPORTE BREVE, CALCULA Y VISUALIZA FIGURAS
                % (-1)  CARGA LA DATA NORMALIZADA CALCULA PERO NO VISUALIZA FIGURAS
                % (-2)  CARGA LA DATA NORMALIZADA CALCULA Y VISUALIZA FIGURAS SOBRE LA EFICIENCIA
% =================================================================================================  
% CARGA CAMPO GEOMAGNETICO YA NORMALIZADO
% ------------------------------------------------------ 
  load Campo_Geomagnetico_Normalizado
% ------------------------------------------------------
% load  Data_Geomagnetica_Simulada
% ------------------------------------------------------
       Data_Stations = 1;
[Data,Label_Station] = Data_Geomagnetica(Data_Stations);
% ==============================================================================================
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
% Uoo = linspace(7,0,12);     
  uo = 1;%0.125/4; 
% -----------------------------
% ESTACIONES A UTILIZAR: (1) KKA (2) KNY (3) MMB y (4) USA LAS 3 ESTACIONES SIMULTANEAMENTE.
%                        (OJO) SE PUDEN USAR COMBINACIONES DE LAS
%                         ESTACIONES ==> [1 2], [1 3], [2 3] o [1 2 3]
% --------------------------------------------------------------------------------------------------------------------------------
  Station_Used_Now = [1 2 3 4];
  
% VALOR DE GAMMA A USAR
% -------------------------
               Gamma_2 = 10;%linspace(1,1e1,10);
           %Gamma_2(1) = 1;
% -------------------------------------------     

% CALCULO DE LA EFICIENCIA PARA CADA ESTACION
% -------------------------------------------
    
for Actual_Station = 1:length(Station_Used_Now)
    
       Stadistic_Ef = [];            % ESTADISTICA DE LA EFICIENCIA (para usar en ciclos de un prametro)
             Max_a = [];
                FP = []; 
                
    if Actual_Station <=3
        Station_Used_Now = Actual_Station;  % SOLO UNA ESTACION POR LAZO
    else
        Station_Used_Now = [1 2 3];         % LAS TRES ESTACIONES SIUMULTANEAS
    end
    
disp(['Station_Used_Now: ',num2str(Station_Used_Now)])   
 
 for jg = 1:length(Gamma_2)  
  disp([CASO Actual_Station Gamma_2(jg)])
         %disp(Gamma_2(jg))
         
% for Hj=1:4
%     disp(Hj)
%  if Hj < 4 
%      Station_Used_Now = [Hj];  
%   end
%  if Hj==4


%  end
% -----------------
% VALORES DE GAMMA  
% -----------------   
      Ng = 4;                                  % NUMEROS DE GAMMAS A GENERAR 
 gamma_1 = Gamma_2(jg);                        % VALOR INFERIOR DE GAMMA                         
 gamma_2 = Gamma_2(jg) + Gamma_2(jg)/1e6;      % VALOR SUPERIOR DE GAMMA
       G = linspace(gamma_1,gamma_2,Ng);       % VALORES DE GAMMA
      dg = 3;                                  % TAMA??O DE LA MUESTRA DE GAMMA (GAMMAS USADOS EN LA ESTADISTICA) (debe ser >= 3)
% --------------------------------------------------------------------------------------------------------------------------------------
% NUMERO DE DIAS POR VENTANA 
% ---------------------------
       Md = 14;                               % NUMERO DE VENTANAS (SE CALCULAN AUTOMATICAMENTE A PARTIR DE Md Y TIEMPO TOTAL DE LOS DATOS)
 % -----------------------------------------------------------------------------------------------------------------------------------
 % ACTIVACION DEL PRECURSOR
 % -------------------------- 
   nd_b = 60/Md; %90/Md;                               % EL SISMO OCURRE (nd_b) ANTES DE LA ACTIVACION DEL PRECURSOR 
   Nd_a = 60/Md; %6=120 dias%650/Md;%(50:20*5:95*Md)/Md;         % EL SISMO OCURRE (nd_a) DESPUES DE LA ACTIVACION DEL PRECURSOR  ( (Nd_a) ETIQUETA DE ANCHO DE VENTANA)  
 % --------------------------------------------d---------------------------------------------------------------------------------------------------
 % INTERVALO HORARIO TOTAL
 % ------------------------------
    
  num_h = 1;                                 % NUMERO DE HORAS DEL INTERVALO EN ESTUDIO
  Num_h = 1/num_h;  
    dho = 1;
     Dh = [1:dho:24*Num_h];                      
    Dho = 2;                                  % INTERVALO INCREMENTAL HORARIO (debe ser >=2)
   M_dh = []; 
 for i=1:length(Dh) - Dho,dh = Dh(i):Dh(i+Dho); M_dh(i) = mean(dh); end          % (M_dh) ETIQUETA DE INTERVALO HORARIO
% ------------------------------------------------------------------------------------------------------------------------
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
% -------------------------------------------------------------------------------------
f1 = Lx(1,1);            % FECHA INFERIOR DE LOS DATOS TOTALES
f2 = Lx(end,end);          % FECHA SUPERIOR DE LOS DATOS TOTALES

% INDICES VALIDOS
% -----------------
      ind_s = find(Lxx(:,1)>=f1 & Lxx(:,1) <= f2);           % INTERVALO DE FECHAS (f1 y f2)
     ind_sf = find(Fo(1,:)>=f1 & Fo(1,:) <= f2);
         Fo = Fo(:,ind_sf);
 [nf,mf,kf] = size(Fo);  
if length(ind_s) > 1
  
    k_ev = k_ev+1;       
      nc = length(ind_s);            % NUMERO TOTAL DE SISMOS EN ESTUDIO
   Lxx_u = Lxx(ind_s,:);
   Lyy_u = Lyy(ind_s,:);

       Ep = [];                      % EFICIENCIA PORCENTUAL
        S = [];                      % ESTATUS
       S1 = [];                      % ESTATUS
    MAG_S = {};                      % MAGNITUDES SELECCIONADAS EN EL CICLO
  DEPTH_S = {};                      % PROFUNDIDADES SELECCIONADAS EN EL CICLO
  
% -------------------------------------------------------------------------
            dga = 1;
             Ku = 0;
      Parametro = Nd_a;
% ====================================================================================================================================================
INDICES = ones(nc,length(Nd_a));
indices = 1:nc;
for ku = 1:length(Parametro)
    ind_im = 0;
 for ga=1:dga:length(G)-dg      
             Ku = Ku + 1;
           nd_a = Nd_a;                       % Parametro(Ku);
          gamma = G(ga:ga+dg);                % VALORES DE GAMMA 
      [gamma,E] = Energia_Vs_Gamma(gamma,0);  % VALORES DE LA ENERGIA EN GeV
    
% SALIDAS
% ========
              Ef = [];            % EFICIENCIA ACIERTOS Y FALLOS
          IND_KC = [];            % INDICE CIERTOS
          IND_KF = [];            % INDICE FALLOS
           AMP_U = [];            % AMPLITUD ASOCIADA AL PRECURSOR EN CADA ACIERTO
              IM = [];
      count_surf = 0;
% ------------------------------------------------------------------------------------
for i=1:length(Dh) - Dho
  count_surf = count_surf + 1;
          dh = Dh(i):Dh(i+Dho);
        m_dh = mean(dh/Num_h);
        % -----------------------------------------------------------
          Programa_General_Precursores_Relativistas                  % CALCULA EL PRECURSOR SELECCIONADO (FLUJO POR EJEMPLO)
          Encuentra_Sismos_Sicronizados                              % SISMOS LOCALIZADOS DIAS ANTES Y DIAS DESPUES DEL MAXIMO DEL PRECURSOR
        % -----------------------------------------------------------
        IND_KC = [IND_KC;IND_kc];                                    % SISMOS ENCONTRADOS CIERTOS
        IND_KF = [IND_KF;IND_kf];                                    % SISMOS NO ENCONTRADOS
            Ef = [Ef;m_dh pc pf pc+pf p_ind_c p_ind_f p_ind_c + p_ind_f CASO Actual_Station jg length(Gamma_2)];   % EFICIENCIA
         AMP_U = [AMP_U;Amp_U];
         
        % surf(P),shading interp,view([0 90]),pause

            Pt(:,:,count_surf) = P';             % EVOLUCION DEL PRECURSOR UMBRALIZADO (SURF)
              Pa(:,count_surf) = As';            % EVOLUCION DEL PRECURSOR SIN UMBRALIZAR (PROBABILIDAD CONDICIONADA)(SURF) 
            
            % ----------------------------------------------------------------------------------------------------------------
            % SOLO LOS INDICES CIERTOS Y FALLOS
            % -------------------------------------------
                            As_Ind_c = zeros(1,length(As));
                            As_Ind_f = As_Ind_c;
                     As_Ind_c(Ind_c) = As(Ind_c);
                    As_Ind_f(IND_FP) = As(IND_FP);
                    
              AS_IND_C(:,count_surf) = As_Ind_c;
              AS_IND_F(:,count_surf) = As_Ind_f;
              AS_IND_CF(:,count_surf) = As_Ind_c + As_Ind_f;
%              if size(AS_IND_C,2) > 1
%               clf
%               hg = 1;
%               set(figure(1),'Position',[680 32 1197 906],'Color','W')
%               subplot(4,1,1)
%               surf(AS_IND_C'),shading interp,view([0 90]),title('AS_IND_C')
%               xlim([0 Nv])
%               subplot(4,1,2)
%               plot((Lxx'-f1)/Md,Lyy','-r','LineWidth',[1]),grid,hold on
%               plot(smooth(sum(AS_IND_C'),hg),'-b','LineWidth',[2]),hold off  % SISMOS NO-ENCONTRADOS
%               xlim([0 Nv])
%               subplot(4,1,3)
%               surf(AS_IND_F'),shading interp,view([0 90]),title('AS_IND_F')
%               xlim([0 Nv])
%               subplot(4,1,4)
%               plot((Lxx'-f1)/Md,Lyy','-r','LineWidth',[1]),hold on
%               plot(smooth(sum(AS_IND_F'),hg),'-r','LineWidth',[2])
%               plot(smooth(sum(AS_IND_C'),hg),'--b','LineWidth',[2])
%               grid,hold off  % SISMOS NO-ENCONTRADOS
%               xlim([0 Nv])
%               Ejes_Visibles(1)
% %               set(figure(2),'Position',[680 32 1197 906],'Color','W')
% %               plot((Lxx'-f1)/Md,Lyy','-r','LineWidth',[1]),hold on
% %               plot(mean(AS_IND_CF').*std(AS_IND_CF')*20),grid,hold off
%               pause(.01)
%              end

% VISULIZACION DE LAS SLIDAS
% --------------------------------------------------------------------------------------------------------------------------------------------------------
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------\n')
disp('     m_dh       pc      pf        pc+pf     p_ind_c    p_ind_f  p_ind_c + p_ind_f      CASO    Actual_Station    jg    length(Gamma_2)')
disp([Ef])
fprintf('--------------------------------------------------------------------------------------------------------------------------------------------------\n')
% --------------------------------------------------------------------------------------------------------------------------------------------------------
end
% =================================================================================================
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
% --------------------------------------------------------
 [j_kc,i_kc] = hist(IND_KC,length(IND_KC));
        i_kc = round(i_kc); 
        p_kc = (j_kc/length(IND_KC))*100;             % PROBABILIDAD local     
     IND_KCO = IND_KF;
     IND_KFO = IND_KC;
% =================================================================================================
% VISUALIZACION DE EFICIENCIA
% ==============================
% ACIERTOS Y FALLOS VERSUS INTERVALO HORARIO
% -------------------------------------------
            Int_H = Ef(:,1);
         Aciertos = Ef(:,2);
           Fallos = Ef(:,3);
          Total_a = Ef(:,4); 
        [jMc,iMc] = max(Aciertos); 
% ---------------------------------------------------------------------------
    Ind_Aciertos = Ef(:,5);
      Ind_Fallos = Ef(:,6);
         Total_f = Ef(:,7);  
       [jMf,iMf] = max(Ind_Fallos);
% ---------------------------------------------------------------------------         
% INDICES FALSOS POSITIVOS VERSUS INTERVALO HORARIO
% --------------------------------------------------
% SISMOS NO ENCONTRADOS
% ELIMINA INDICES REPETIDOS
% ==========================
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
 % ---------------------------------------------------    
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
% --------------------------------------------
               St = [1 2 3];     % ESTACIONES
%
% ESTIMACION DE LA DISTANCIA A ESTACIONES USADAS Y NUEVOS EVENTOS SEGUN EL PARAMETRO DE CONTROL USADO
% -----------------------------------------------------------------------------------------------------------
      dx_station = y_station(Station_Used_Now)/grados_km;
      dy_station = x_station(Station_Used_Now)/grados_km; 
    Dx_Station_m = mean(dx_station);
    Dy_Station_m = mean(dy_station);
% ---------------------------------------------------------------------------------------------
% DISTANCIA ENTRE EL CENTRO GEOMETRICO ENTRE ESTACIONES USADAS Y EVENTOS
% ------------------------------------------------------------------------------------------------------------------------------------- 
if length(IND_FO)>0
    R_Dx_Station = repmat(Dx_Station_m,length(Ind_New_Events),1);
    R_Dy_Station = repmat(Dy_Station_m,length(Ind_New_Events),1); 
     Rx_St_Event = sqrt((R_Dx_Station - Lon_S(Ind_New_Events)).^2 + (R_Dy_Station - Lat_S(Ind_New_Events)).^2);
     Rm_St_Event = (Rx_St_Event)*grados_km; % DISTANCIA MEDIA NUEVOS EVENTOS RESPECTO AL CENTRO GEOMETRICO DE LAS ESTACIONES USADAS
else
    Rm_St_Event = 0;
end

% -------------------------------------------------------------------------------------------------------------------------------------
% SALIDA: EFICIENCIA PORCENTUAL EVENTOS CIERTOS ENCONTRADOS
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
                    Ep(Ku) = 100 - (length(Ind_S)/length(ind_s))*100;
           MAG_S(Ku).MAG_S = Mag_S(Ind_New_Events);
       DEPTH_S(Ku).DEPTH_S = Depthkm_S(Ind_New_Events);
     AMPLITUD(Ku).AMPLITUD = AMP_U;
     IMAGE_HISTORY(Ku).IMAGE_HISTORY = IM;
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------           
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
% ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
     
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

% SALIDA
% -------------------------------
           %save Status_Actual_Variando_Umbral S
           %save Status_Actual_Variando_Nd_a S
           % save Status_Actual_Variando_Uo_kak 
           %disp(S)
          % pause(.01)
end
 end
end
ST(k_ev).ST = S;
%disp(S)

end

end
% if Hj==1
%    save Ejemplo_NRA_new_kak
% elseif Hj==2
%   save Ejemplo_NRA_new_kny
% end
% if Hj == 3
%   save Ejemplo_NRA_new_mmb
% elseif Hj == 4
%     save Ejemplo_NRA_new_3St % CONTIENE LAS TRES ESTACIONES
 %end

%end 
%save Stadistic_Gamma_New_III Stadistic_Ef Gamma_2 Md Nd_a num_h dg Ng uo Momento_Usado Componente_Tensor_Maxwell ...
%     Station_Used_Now nd_b iMc Aciertos Pt Int_H Fallos Total_a Fallos iMc Pt Mag Num_h Nv

 
 
% save Stadistic_Gamma_New_III        % SALVA TODAS LAS VARIABLES Y DATA con gamma = 10
% save Stadistic_Gamma_New_IV           % SALVA TODAS LAS VARIABLES Y DATA con gamma = 100

% save Stadistic_Gamma_New_I % GUARDA GAMMA DESDE 1 A 1E3

   if Caso_Ent_Prod == 0
       Pa = AS_IND_C;
       
       Stadistic_Entrena{Actual_Station,jg} =  struct('Stadistic_Ef',Stadistic_Ef,'Gamma_2',Gamma_2,'Md',Md,'M_dh',M_dh,'Nd_a',Nd_a,'num_h',num_h,'dg',dg,'Ng',Ng,'uo',uo,...
                                                      'Momento_Usado',Momento_Usado,'Total_a',Total_a,'Componente_Tensor_Maxwell',Componente_Tensor_Maxwell,...
                                                      'Station_Used_Now',Station_Used_Now,'iMc',iMc,'Int_H',Int_H,'Aciertos',Aciertos,'Fallos',Fallos,'t',t,...
                                                      'Label_Station',Label_Station,'dho',dho,'Dh',Dh,'nh',nh,'Num_h',Num_h,'f1',f1,'f2',f2,'Nv',Nv,...
                                                      'M_Gamma',M_Gamma,'Max_a',Max_a,'Caso_Ent_Prod',Caso_Ent_Prod,'Pt',Pt,'Pa',Pa);
           save STADISTIC_ENTRENA_R_Gamma_10_1 Stadistic_Entrena
   else
       Stadistic_Produce{Actual_Station,jg} =   struct('Stadistic_Ef',Stadistic_Ef,'Gamma_2',Gamma_2,'Md',Md,'M_dh',M_dh,'Nd_a',Nd_a,'num_h',num_h,'dg',dg,'Ng',Ng,'uo',uo,...
                                                       'Momento_Usado',Momento_Usado,'Total_a',Total_a,'Componente_Tensor_Maxwell',Componente_Tensor_Maxwell,...
                                                       'Station_Used_Now',Station_Used_Now,'iMc',iMc,'Int_H',Int_H,'Aciertos',Aciertos,'Fallos',Fallos,...
                                                       't',t,'Label_Station',Label_Station,'dho',dho,'Dh',Dh,'nh',nh,'Num_h',Num_h,'f1',f1,'f2',f2,'Nv',Nv,...
                                                       'M_Gamma',M_Gamma,'Max_a',Max_a,'Caso_Ent_Prod',Caso_Ent_Prod,'Pt',Pt,'Pa',Pa);
            save STADISTIC_PRODUCE_R_Gamma_10_1 Stadistic_Produce
   end
   
 end % FIN DE CASO
 

end


end
 
% SALIDA EN CADA LAZO
% -------------------------------------------------------------------
   save STADISTIC_ENTRENA_R_Gamma_10_1 Stadistic_Entrena
   save STADISTIC_PRODUCE_R_Gamma_10_1 Stadistic_Produce
% -------------------------------------------------------------------
% save IMAGE IMAGE_HISTORY
% ====================================================================================================================================================