% prueba_efectividad_sensibilidad_roc_respaldo.m
% ================================
% DESCRIPCIÓN: Programa las salidas de Precursor_Activation_Gamma_Dependence_new
% y calcula curva ROC para una determinada ventana temporal antes
% (vent_an) o despues  (vent_des).
% 
% SALIDAS: 
% dist_i_gamma(umbral,gamma),umbral_optimo, gamma_optimo
% Positive_predictive_value(cont_umbrales,:) verdaderos pos /verd pos + fal_pos
% False_predictive_value(cont_umbrales,:) verdaderos negativos / verdaderos negativos + falsos negativos
% sensibilidad(cont_umbrales,:)   SENSIBILIDAD X GAMMA
% especificidad(cont_umbrales,:)  ESPECIFICIDAD X GAMMA
%==========================================================================
close all
clearvars  -except T_entrena T Lx Lxx Lx_Cont Ly_Cont Ly IA_MEAN_YS k2 IRISIDS Data_Stations f1 f2 final ini Gamma_2 Nombre_caso Data_Stations
band = 2;
vent_an = 120;%ventana antes
vent_des =0;%ventana despues
% 
%function [umbral_optimo,gamma_optimo,Positive_predictive_value,sal_s,...
% especificidad,sensibilidad]=prueba_efectividad_sensibilidad_roc_respaldo(band,vent_an,vent_des,T_entrena,...
%  T,Lx,Lxx,Lx_Cont,Ly_Cont,Ly,IA_MEAN_YS,k2,Data_Stations,f1,f2,final,Nombre_caso)
% close all
% =========================================================================
Umbr_limit_sup=floor(max(max(IA_MEAN_YS))) +1;
cont_umbrales = 0;
cont_umbrales = 0;
Loc_sismos= Lx'-f1;%todos los sismos
Loc_sismos_ent=Lx_Cont'-f1; %sismos entrenamiento
Mag_prueba = Ly';
count_vn = 0;
ind_vecVN_T{1,1} =0;
ind_vecFN_T = [];
ind_vector_T = [];
display_count = 1;

for Umbral_positivo = 0:1:Umbr_limit_sup
 cont_umbrales = cont_umbrales+ 1;
 
 if band == 1;
  figure
  plot(Loc_sismos,Ly','-b','LineWidth',[2]), hold on %plot todos los sismos
  plot(Loc_sismos_ent,Ly_Cont','--g','LineWidth',[3]) ; hold on;
  xlabel('Tiempo') ;ylabel('Precursor y lineas de sismos')
  plot([0 T(end)],[Umbral_positivo/k2  Umbral_positivo/k2],'-','Color',[0.145,1- 0.28 , 0.124],'LineWidth',[2]), hold on %plot umbral
 end
 % CALCULA EFICIENCIA Y OTRAS COSAS DE LAS SALIDAS DE Precursor_Activation_Gamma_Dependence_new
 eficiencia_antes =[];
 an_logic_T = [];
 eficiencia= [];
 ef2 = [];
 zok = [];
 
 if band == 1;
  display(' --------------------------------------------------------------')
  display(' Abbreviations:')
  display(' activaciones "precursoras": AP')
  display(' numero T de datos: NT')
  display(' activaciones Totales: AT')
  display(' --------------------------------------------------------------')
  display(' ==============================================================')
  display('Calcula eficiencia con data real:')
 end
 
 FP_logic_T =[];
 %VARIABLES VERDADEROS NEGATIVOS:===================================
 non = []; %variable que contendra logicos para verdaderos negativos
 non_T =[]; %variable que contendra logicos para verdaderos negativos totales
 % =================================================================
 last_sis = size(Lx,1); %INDICE del ultimo sismo
 %ESTO SOLO SE MODIFICA SI T(end)>= sis
 
 
 %% VA X CADA GAMMA EN LA VARIABLE IA_MEAN_YS
 for zoo=1:size(IA_MEAN_YS,1) %for va por los gammas
  if band == 1;
   %figure(1)
   % get(figure(1),'Position')
   %plot(Loc_sismos),Ly','-b','LineWidth',[2]), hold on %plot todos los sismos
   %plot((Lx_Cont'-Lxx(1,1)),Ly_Cont','--g','LineWidth',[3]) ; hold on;
   %plot([0 Loc_sismos(end)],[Umbral_positivo/k2  Umbral_positivo/k2],'-g','LineWidth',[2]), hold on %plot todos los sismos
  end
  % ------ jo, indices mayores a Umbral_positivo en los datos de IA_MEAN_YS -----------
  [~,jo]=find((IA_MEAN_YS(zoo,:)) > Umbral_positivo);%salidas: jo, indices mayores a Umbral_positivo
  % ------ jo_ig, indices iguales a Umbral_positivo en los datos de IA_MEAN_YS -----------
  %[~,jo_ig]=find((IA_MEAN_YS(zoo,:)) == Umbral_positivo);
  % ------- jo_me, indices menores a Umbral positivo en los datos
  [~,jo_me]=find((IA_MEAN_YS(zoo,:)) <= Umbral_positivo);
  % -------
  
  %verdaderos negativos:
  primero = 0; % para moverme entre los sismos
  %fin = 0;
  non_T = zeros(size(Loc_sismos,2),length(T(jo_me)));
  % ====================
  if band == 1;
   plot(T,IA_MEAN_YS(zoo,:)/k2,'.','LineWidth',[3])
  end
  an_logic_T = zeros(size(Loc_sismos,2),length(T(jo)));
  %% SE MUEVE POR CADA SISMO EN EL REGISTRO
  for i =1:size(Loc_sismos,2) %for por los sismos
   sis=Loc_sismos(1,i); %toma un sismo particular
   %% VERDADEROS POSITIVOS POR CADA GAMMA:
   % ================ Busca en el vector T(jo) ===============
   an_logic = and( T(jo)>= sis-vent_an,T(jo) <= sis+vent_des); %LOGICA PARA ANTES Y DESPUES
   % Guarda los indices logicos en la matrix an_logic_T
   an_logic_T = [an_logic_T;an_logic]; %Guarda los indices
   if sum(an_logic) == 0
    ind_vector_T(i,zoo).datos = 0;
   else
    ind_vector_T(i,zoo).datos = (jo(an_logic));
   end
   if band == 1;
    plot(T(jo(an_logic)),(IA_MEAN_YS(zoo,(jo(an_logic)))/k2),'*r','LineWidth',[3])
   end
   % ================= FIN DE VERDADEROS POSITIVOS POR GAMMA =
   %% VERDADEROS NEGATIVOS POR CADA GAMMA:
   if  i >= size(Loc_sismos,2) %Loc_sismos(1,i+1)>= T(end)
    fin = T(end);
   else
    fin = Loc_sismos(1,i+1)-vent_an;
   end
   % ================ Busca en el vector T(jo_me)=============
   % [~,jo_me]=find((IA_MEAN_YS(zoo,:)) <= Umbral_positivo);
   log1 = T(jo_me)< (sis-vent_an);
   log2 = T(jo_me) > primero;
   
   non = or (and(log1, log2),and(T(jo_me)>(sis+vent_des),T(jo_me) < fin));%Logico verdaderos negativos:
   non_T = [non_T;non]; %Guarda los indices verdaderos negativos para ese gamma y lo almacena
   
   
   
   if primero >= T(end)
    primero = T(end);
   else
    primero = sis+vent_des;
   end
   
   if (zoo== 1)
    ind_v_n=(jo_me(non));
    for vn =1:length(ind_v_n)
     vector_log_vn = IA_MEAN_YS(:,ind_v_n(vn));
     max_vec_log_vn= max(IA_MEAN_YS(:,ind_v_n(vn)));
     if (max_vec_log_vn<Umbral_positivo)
      %,not(isempty (ind_v_n(vn)))))
      if not(isempty(ind_v_n(vn)))
       %display([i,vn,cont_umbrales])
       ind_vecVN_T{i,vn,cont_umbrales}= ind_v_n(vn) ;
       
      end
      
      
      if band == 1;
       plot(T(ind_v_n(vn)),vector_log_vn/k2,'.m','LineWidth',[10])%MARCA PUNTOS verdaderos negativos para todos los gammas
       % pause(0.5)
      end
     else
      ind_vecVN_T{i,vn,cont_umbrales}= 0 ;
     end
     %pause(0.01)
    end
   end
   
   %% CALCULA LOS FALSOS NEGATIVOS. ==========================
   
   an_logic_FN = and( T(jo)>= sis-vent_an,T(jo) <= sis+vent_des); %LOGICA PARA ANTES Y DESPUES
   an_logic_FN2= and( T(jo_me)>= sis-vent_an,T(jo_me) <= sis+vent_des);
   log_control=sum(an_logic_FN);
   log_contol2= (sum(IA_MEAN_YS(:,jo_me(an_logic_FN2))>Umbral_positivo)>=1);
   if  T(end)>= sis %cuando llega al limite del calculo del precursor, se detiene
    if log_control == 0
     zok2(zoo,i) = 0;  % VENTANA ESTA VACIA TOTAL algun  GAMMA
     ind_vecFN_Txgamma{zoo,i} = jo_me(an_logic_FN2);
     if  and(log_control == 0,sum(log_contol2)==0) %%isempty(an_logic)%and( isempty(an_logic), isempty(IA_MEAN_YS(zoo,(jo(an_logic))))
      zok(zoo,i) = 0; % VENTANA ESTA VACIA TOTAL DE TODOS LOS GAMMAS
      ind_vecFN_T{zoo,i} = jo_me(an_logic_FN2);
     else
      zok(zoo,i) = 1; % VENTANA NO ESTA VACIA
      ind_vecFN_T{i,zoo}= 0;
     end
    else
     ind_vecFN_Txgamma{zoo,i} = 0;
     zok(zoo,i) = 1; % VENTANA NO ESTA VACIA
     zok2(zoo,i) = 1;
     ind_vecFN_T{i,zoo}= 0;
    end
   else
    last_sis = i -1;
    break
    %i = size(Loc_sismos,2) +1; %fuerza la salida del loop
   end
   % FIN DE CALCULA LOS FALSOS NEGATIVOS.=====================
   %pause(0.1)
   %display(an_logic)
 end
  
  %% EFICIENCIA PRECURSOR:(verdaderos positivos)
  % -------- EXTRA: ----------------
  indT_VP_cell{zoo}=unique([ind_vector_T(:, zoo).datos]);
  % CONTIENE UNA ESTRUCTURA CON LOS INDICES POSITIVOS
  [io,~]=find(sum(an_logic_T)>0); %suma para quitar redundancias
  suma_aciertos_T=sum(io);%cantidad total de aciertos x gamma
  eficiencia(zoo) = suma_aciertos_T / length(T(jo))*100;
  % =============================================================
  
  %% EFICIENCIA  FALSOS POSITIVOS. ==============================
  % CALCULA LOS FALSOS POSITIVOS. =================
  indT_FP_cell{zoo}=setdiff(jo,indT_VP_cell{zoo}); %indices de los FP totales en el vector T
  % size(indT_FP_cell)
  %EFICIENCIA FALSOS POSITIVOS
  [io2,~]=find(sum(an_logic_T)==0);
  suma_fallos_T(zoo)=sum(io2);
  %total = suma_fallos_T(zoo) + suma_aciertos_T;
  eficiencia_falsos_positivos(zoo) = suma_fallos_T(zoo) / length(T(jo))*100;
  
  %% EFICIENCIA VERDADEROS NEGATIVOS.============================
  [io,zop]=find(sum(non_T)>0); %suma para quitar redundancias
  suma_non_T(zoo)=sum(io);%cantidad total de non para un gamma
  eficiencia_non(zoo) = suma_non_T(zoo)/ length(T(jo_me))*100;
  % ====================================================================
  %Cantidad de disparos divido en la cantidad de datos total:
  %Para dar un sentido de cuanto se suele disparar en el total de los
  %datos
  %activaciones totales(falsos y positivos) / numero total de datos
  ef2(zoo) = length(T(jo))/ length(T) *100;%
  %activaciones "precursoras" / numero total de datos
  ef3(zoo) = suma_aciertos_T/length(T) *100;
  if band == 1;
   display(['Gamma:' num2str(zoo) ' eficiencia: ' num2str(round(eficiencia(zoo),2)) '% //' ...
    ' AP / NT: ' num2str(round(ef3(zoo),2)) '% //' ...
    ' AT/ NT: ' num2str(round(ef2(zoo),2)) '%' ])
  end
  
  verd_pos(cont_umbrales,zoo)=suma_aciertos_T; %CONTEO TOTAL DE VERDADEROS POSITIVOS X GAMMA X UMBRAL
  Suma_verdpos_falpos(zoo)=length(T(jo)); %VERDADEROS POSITIVOS + FALSOS POSITIVOS
  Fal_pos(cont_umbrales,zoo) = suma_fallos_T(zoo); %FALSOS POSITIVOS
  
  %VERDADEROS NEGATIVOS:
  % activaciones totales(falsos y positivos) / numero total de datos
  ef2_non(zoo) = length(T(jo_me))/ length(T) *100;
  %activaciones NO "precursoras" / numero total de datos
  ef3_non(zoo) = suma_non_T(zoo)/length(T) *100;
  %
  if band == 1;
   display(['Gamma:' num2str(zoo) ' verd negat: ' num2str(round(eficiencia_non(zoo),2)) '% //' ...
    ' ANP / NT: ' num2str(round(ef3_non(zoo),2)) '% //' ...
    ' ANT/ NT: ' num2str(round(ef2_non(zoo),2)) '%' ])
  end
 end
 % GUARDADO DE SALIDAS POSIBLEMENTE IMPORTANTES:
 % verd_postitivos %VERDADEROS POSITIVOS X GAMMA
 eficiencia_p1= eficiencia;%eficiencia calculo directo de la ef con parametros de entrada
 ef2_p1=ef2;%ef2 calculo directo de la ef con parametros de entrada
 ef3_p1=ef3;%ef3 calculo directo de la ef con parametros de entrada
 
 %% FALSOS NEGATIVOS X GAMMA:
 %ind_vecFN_T{zoo,i}
 %cuenta la cantidad de SISMOS CON FALSOS NEGATIVOS X GAMMA
 %        count_sismico = zeros(1,size(IA_MEAN_YS,1));
 Suma_FN = [];
 FN_T = [];
 for zoo = 1:size(zok,1)
  %Nota importante: si aquí sale un erro va a ser pq zok no tiene
  %las dimensiones correctas o empieza a asignar puntos que no
  %existan dependiendo del umbral
  indxgammaVN = unique(nonzeros([ind_vecFN_Txgamma{zoo,:}]));
  for i =1:size(zok,2) %for por los sismos
   %control = 0;
   if (zok(zoo,i))== 0 % si es 0 significa no detecto señal antes/despues del sismo
    
   end
  end
  %Suma_FN = prueba123;
  if isempty(indxgammaVN)% (zok(zoo,:))~= 0% or(isempty(Suma_FN),(zok(zoo,i))~=0)
   FN_T(zoo) = 0;
   % pause
  else
   %[io,~] = find( (Suma_FN(zoo,:)) > 0);
   
   FN_T(zoo)=length(indxgammaVN); %falsos negativos en el precursor x gamma
  end
  %     Ef_FN(zoo) = FN_T(zoo)/length(T(jo_me)).*100;%eficiencia % falsos negativos en el precursor x gamma
 end
 %         % SALIDA:falsos_neg contiene la cantidad de FALSOS NEGATIVOS X GAMMA
 %% CUENTA SISMOS NO DETECTADOS X GAMMA:
 
 %ind_vecFN_T{zoo,i}
 %cuenta la cantidad de SISMOS CON FALSOS NEGATIVOS X GAMMA
 count_sismico = zeros(1,size(IA_MEAN_YS,1));
 
 for zoo = 1:size(zok2,1)
  for i =1:size(zok2,2) %for por los sismos
   if (zok2(zoo,i))== 0 % si es 0 significa no detecto señal antes/despues del sismo
    count_sismico(zoo) =  count_sismico(zoo) + 1;
   end
  end
 end
 sis_no_detect=count_sismico;
 ef_sismica_Ef = sis_no_detect/last_sis .*100;% de sismos no detectados
 %% cuenta la cantidad de SISMOS CON FALSOS NEGATIVOS X TODOS LOS GAMMAS
 indT_FN_cell = unique(nonzeros([ind_vecFN_T{:,:}]));
 count_FALSO_NEGATIVO_TOTAL = 0;

 for i =1:size(zok,2)%for por los sismos
  if sum(zok(:,i))==0
   sis=Loc_sismos(1,i); %toma un sismo
   count_FALSO_NEGATIVO_TOTAL =   count_FALSO_NEGATIVO_TOTAL + 1;
   %                 %como ya se que todos los gammas tendran 0 activacion:
   for zoo = 1:size(zok,1)
    % ------ jo, indices mayores a Umbral_positivo en los datos de IA_MEAN_YS -----------
    [~,jo_me]=find((IA_MEAN_YS(zoo,:)) <= Umbral_positivo);
    % ================= Busca en el vector T(jo)
    an_logic = and( T(jo_me)>= sis-vent_an, T(jo_me) <= sis+vent_des); %LOGICA PARA ANTES Y DESPUES
    if band == 1;
     plot(T(jo_me(an_logic)),(IA_MEAN_YS(zoo,(jo_me(an_logic)))/k2),'*k','LineWidth',[1])%MARCA PUNTOS NO PRECURSORES
    end
   end
   if band == 1;
    plot([Loc_sismos(1,i) Loc_sismos(1,i)],[0 Mag_prueba(2,i)],'-k','LineWidth',[3]), hold on %plot sismo no detectado
    %                     %pause
   end
  end
 end
 if indT_FN_cell == 0
  FN_T_2 = 0;
  Ef_FN_2 = 0;
  Fals_neg_pr=0;
 else
  FN_T_2 = length(indT_FN_cell);%conteo del precursor falso negativo
  Ef_FN_2 = FN_T_2/length(T(jo_me)).*100;    %porcentaje del precursor falso negativo
  Fals_neg_pr =  count_FALSO_NEGATIVO_TOTAL/last_sis  .*100;%porcentaje de sismos falsos negativos
 end
 if band == 1;
  display(['PORCENTAJE de sismos no detectados x gamma: ' num2str(round(ef_sismica_Ef,2))])
  display(['PORCENTAJE de sismos no detectados X TODOS LOS GAMMAS: '  num2str(round(Fals_neg_pr,2))  ])
  display(['TOTAL DE SISMOS NO DETECTADOS: '  num2str(round( count_FALSO_NEGATIVO_TOTAL))  ])
  display(['TOTAL DE SISMOS: '  num2str(round(last_sis))  ])
 end
 sal_s.ef_xgamma{cont_umbrales} = 100 - ef_sismica_Ef; % %de sismos si detectados x gamma
 sal_s.ef_tgamma(cont_umbrales) = 100 - Fals_neg_pr;%de sismos si detectados T gamma
 sal_s.total = last_sis;
 sal_s.tsno_detec(cont_umbrales)= count_FALSO_NEGATIVO_TOTAL;
 
 
 %% VERDADEROS POSITIVOS TOTALES
 for zoo = 1:size(ind_vector_T,2)
  indT_VP_cell{zoo}=unique(nonzeros([ind_vector_T(:, zoo).datos]));
 end
 if zoo > 1
  ind_verd_pos_T=unique([indT_VP_cell{1}]);%indices verdaderos positivos totales
 else
  ind_verd_pos_T=unique([indT_VP_cell{:}]);
 end
 if band == 1
  
  plot(T(ind_verd_pos_T),IA_MEAN_YS(:,ind_verd_pos_T)/k2,'.g','MarkerSize',20)
  
 end
 VP_T_COUNT = find(ind_verd_pos_T>0);
 % FIN DE LOS VERDADEROS POSITIVOS TOTALES, salida:ind_verd_pos_T
 % ==================================================================
 %% FALSOS POSITIVOS TOTALES:
 
 ind_fal_pos_T=unique(nonzeros([indT_FP_cell{:}]));%indices verdaderos positivos totales

 
 if band ==1
  plot(T(ind_fal_pos_T),IA_MEAN_YS(:,ind_fal_pos_T)/k2,'om','MarkerSize',5)
  plot(T(ind_fal_pos_T),IA_MEAN_YS(:,ind_fal_pos_T)/k2,'xm','MarkerSize',5)
 end
 [FP_T_COUNT,~] = find(ind_fal_pos_T>0);
 %% VERDADEROS NEGATIVOS TOTALES:
 if size(ind_vecVN_T,3) == cont_umbrales
  ind_VN_T = unique(nonzeros([ind_vecVN_T{:,:,cont_umbrales}]));
 else
  indT_FN_cell = 0;
 end

 [~,VN_T_COUNT] = find(ind_VN_T>0);
 if band ==1
  plot(T(ind_VN_T),IA_MEAN_YS(:,ind_VN_T)/k2,'ob','MarkerSize',5)
 end
 
 %SALIDA: count_sismico contiene la cantidad de FALSOS NEGATIVOS X TODOS LOS GAMMAS
 %% TOTALES de todos los gammas: VP_T FP_T FN_T_3 VN_T
 VP_T(cont_umbrales) = length(VP_T_COUNT) ; %verdaderos positivos totales
 
 VN_T(cont_umbrales) = length(VN_T_COUNT); %verdaderos negativos TOTALES
 %       VN_T(cont_umbrales) = sum(suma_non_T(cont_umbrales,:)); %verdaderos negativos TOTALES suma_non_T
 
 FN_T_3(cont_umbrales) =  FN_T_2; %Falsos negativos totales;
 
 %      FP_T(cont_umbrales) = length(T) - FN_T_3(cont_umbrales) - VN_T(cont_umbrales) - VP_T(cont_umbrales) ; %falsos positivos totales
 FP_T(cont_umbrales) = length(FP_T_COUNT); %falsos positivos totales
 %FN_T_3(cont_umbrales) = length(T) FN_T_2; %Falsos negativos totales;
 
 
 verd_t(cont_umbrales) =VP_T(cont_umbrales) + VN_T(cont_umbrales);
 falsos_t(cont_umbrales) =FP_T(cont_umbrales) +FN_T_3(cont_umbrales);
 total2 = verd_t + falsos_t;
 %
 
 %Esta seccion es pq si detecta exitosamente un sismo y otros puntos
 %en la ventana no son VP entonces no sabe como catalogar los que
 %esten por abajo del umbral, así que dire que si 1 punto es asociado
 %a un sismo, entonces toda la ventana tambien será asociada al sismo
 
 if not(total2(cont_umbrales) == length(T))
  
  if display_count == 1
   display('aviso: total es diferente a T')
   display('PENDIENTE TIENES ESTO ACTIVADO: Si se desea que no cuente de esta manera, comentar:')
   display_count = 2;
  end
  if 0 <(length(T)-total2(cont_umbrales))
   %
   VP_T(cont_umbrales) =  VP_T(cont_umbrales) + length(T) - total2(cont_umbrales);
   verd_t(cont_umbrales) =VP_T(cont_umbrales) + VN_T(cont_umbrales);
   total2 = verd_t + falsos_t;
  else
   %                display('tienes mas totales que T')
   %                return
  end
 end
 % fin TOTALES de todos los gammas
 %% SALIDAS X GAMMA:
 % calculo de PPV, FPV, Sensibilidad y especificidad X gammas
 Positive_predictive_value(cont_umbrales,:)=verd_pos(cont_umbrales,:)./(verd_pos(cont_umbrales,:) +Fal_pos(cont_umbrales,:));%verdaderos pos /verd pos + fal_pos
 False_predictive_value(cont_umbrales,:)=suma_non_T./(suma_non_T +FN_T); % verdaderos negativos / verdaderos negativos + falsos negativos
 sensibilidad(cont_umbrales,:) = verd_pos(cont_umbrales,:)./(FN_T+verd_pos(cont_umbrales,:)); %SENSIBILIDAD X GAMMA
 especificidad(cont_umbrales,:) = suma_non_T./ ( suma_non_T + Fal_pos(cont_umbrales,:));% verdaderos negativos /verdaderos negativos + falsos positivos X GAMMA
 % FIN DE SALIDAS X GAMMA
 if band == 1;
  
  display(['SENSIBILIDAD x gamma : '  num2str(round(sensibilidad(cont_umbrales,:).*100,1))  ])
  display(['ESPECIFICIDAD x gamma : '  num2str(round(especificidad(cont_umbrales,:).*100,1))  ])
  display(['TASA DE FALSOS NEGATIVOS : '  num2str(round(100-(sensibilidad(cont_umbrales,:).*100),1))  ])
  display(['TASA DE FALSOS POSITIVOS : '  num2str(round(100-(especificidad(cont_umbrales,:).*100),1))  ])
 end
 
 %cuenta la cantidad de SISMOS CON VERDADEROS POSITIVOS X TODOS LOS GAMMAS
 count_verd_pos = 0;
 
 for i =1:size(zok,2)%for por los sismos
  if sum(zok(:,i))> 0
   count_verd_pos =   count_verd_pos + 1;
   if band == 1;
    plot([Loc_sismos(1,i) Loc_sismos(1,i)],[0 Mag_prueba(2,i)],'-g','LineWidth',[3]), hold on %plot sismo no detectado
   end
   %             display(count_verd_pos)
   %             pause(0.5)
  end
 end
 if band == 1;
  display(['TOTAL DE SISMOS  DETECTADOS: '  num2str(round( count_verd_pos))  ])
 end
 %
 %     plot(1-especificidad(cont_umbrales),sensibilidad(cont_umbrales),'or');hold on; plot([0 1],[0 1],'-r')
 %     plot(1-especificidad(cont_umbrales),sensibilidad(cont_umbrales),'xb')
 %% calcula punto en la curva ROC mas cercano al [0,1]
 dist_i_gamma = zeros(size(sensibilidad,1),size(sensibilidad,2));
 for gamma = 1:size(sensibilidad,2)
  for umbral = 1:size(sensibilidad,1)
   X = [0,1;1-especificidad(umbral,gamma),sensibilidad(umbral,gamma)];
   dist_i_gamma(umbral,gamma)=pdist(X,'euclidean');
  end
 end
 [valor,pos]=min(dist_i_gamma);
 [umbral_optimo, gamma_optimo] = find(dist_i_gamma == min(valor));
 % DESCOMENTAR LUEGO LA LINEA DE ARRIBA.
 %salida: dist_i_gamma(umbral,gamma),umbral_optimo, gamma_optimo
 
end %FIN DEL LAZO PARA LOS UMBRALES
%%  SALIDAS GRAFICADAS X GAMMA
% dist_i_gamma(umbral,gamma),umbral_optimo, gamma_optimo
% sal_s.
% Positive_predictive_value(cont_umbrales,:) %verdaderos pos /verd pos + fal_pos
% False_predictive_value(cont_umbrales,:)  % verdaderos negativos / verdaderos negativos + falsos negativos
% sensibilidad(cont_umbrales,:)   %SENSIBILIDAD X GAMMA
% especificidad(cont_umbrales,:)  %ESPECIFICIDAD X GAMMA
if band == 2
 ho_plot = 5; %factor de suavización para las figuras ROC
 figure
 hold on; plot([0 1],[0 1],'-r')
 legend_asig =[];
 % Place axes at (0.1,0.1) with width and height of 0.8
 ax = gca;
 grid on
 ax.GridAlpha = 0.8;
 ax.GridLineStyle = '--';
 % Main plot
 for gamma = 1:size(sensibilidad,2)
  %guarda las figuras para poder hacer una legenda personalizada
  ggg_p(gamma)=plot(smooth(1-especificidad(:,gamma),ho_plot),smooth(sensibilidad(:,gamma),ho_plot),'LineWidth',2);
  %calculo el area bajo cada curva:
  Q(gamma)=trapz(especificidad(:,gamma),sensibilidad(:,gamma));
 end
 Ejes_Visibles(1)
 
 legend_asig  = strsplit(num2str(round(Gamma_2),2));
 repgamma = repmat(cellstr('\gamma'),[1 size(sensibilidad,2)]);
 repAUC = repmat(cellstr(' AUC: '),[1 size(sensibilidad,2)]);
 AUCgammas=cellstr(strsplit((num2str(round(Q,2)))));
 legend_asig=strcat(legend_asig,repgamma,repAUC,AUCgammas);
 legend([ggg_p],[legend_asig], 'Location', 'northeastoutside')
 % plot([0.5 0.5],[0 1],'--r') %Linea referencia
 % plot([0 1],[0.5 0.5],'--r') %Linea referencia
 
 title (['  Curva ROC para cada \gamma, ' char(Nombre_caso(Data_Stations)) ' ' num2str(vent_an) ' días antes'])
 %title ([' ROC curve for each \gamma,' num2str(vent_an) ' days before']) %Si se desea que no diga el sitio
 xlim([0 1]);ylim([0 1]);xlabel('Razón de falsos positivos') ;ylabel(['Razón de verdaderos positivos']);
 plot([0 1],[0 1],'--r','LineWidth',2)
 grid on
 
 %% Plot el coeficiente AUC
 
 % Place second set of axes on same plot
 
 handaxes2 = axes('Position', [0.5 0.2 0.26 0.29]);
 set(handaxes2, 'Box', 'off')
 
 % Pone la figura
 plot(Gamma_2,Q,'LineWidth',2), hold on
 ylim([0 1])
 ax = gca;
 grid on
 ax.GridAlpha = 0.8;
 ax.GridLineStyle = '--';
 title('AUC de curvas ROC')
 xlabel('Valores \gamma en la curva ROC') ;ylabel(['AUC']);
 hold on
 plot([0 100],[0.5 0.5],'--r','LineWidth',2)
 grid minor
 %%
 figure
 plot(False_predictive_value,Positive_predictive_value,'.r');hold on
 plot(False_predictive_value,Positive_predictive_value,'ob')
 plot([0.5 0.5],[0 1],'--r')
 plot([0 1],[0.5 0.5],'--r')
 xlim([0 1]);ylim([0 1]);xlabel('FPV ') ;ylabel(['PPV'])
 %% CALCULO DE LA PRECISION X GAMMA:[ min=0,max =4]
 precision_xgamma=sensibilidad+especificidad+Positive_predictive_value+False_predictive_value;
 % =========================================================================
 %% calculo de PPV, FPV, Sensibilidad y especificidad TODOS los gammas
 % ENTRADAS NECESARIAS:
 % =========================================================================
 % VP_T
 % FP_T
 % FN_T_3
 % VN_T
 % =========================================================================
 Positive_predictive_value_T = VP_T./(VP_T +FP_T);
 False_predictive_value_T = VN_T./(VN_T +FN_T_3);
 sensibilidad_T = VP_T./(FN_T_3+VP_T);
 especificidad_T = VN_T./ (VN_T + FP_T );
 %% ESPECIFICIDAD Y SENSIBILIDAD TOTAL DE TODOS LOS GAMMAS
 % ho_plot = 10;
 figure
 plot(smooth(1-especificidad_T,10),smooth(sensibilidad_T,10),'r','LineWidth',2);hold on;
 plot([0.5 0.5],[0 1],'--r')
 plot([0 1],[0.5 0.5],'--r')
 plot([0 1],[0 1],'-r')
 for umbral = 1:size(sensibilidad_T,2)
 end
 grid on
 %title ('Curva ROC for all \gamma')
 title (['Curva ROC para todos los \gamma, ' char(Nombre_caso(Data_Stations)) ' ' num2str(vent_an) ' días antes'])
 xlim([0 1]);ylim([0 1]);xlabel('Razón de falsos positivos') ;ylabel(['Razón de verdaderos positivos']);
 Ejes_Visibles(3)
 % FIN DE ESPECIFICIDAD Y SENSIBILIDAD TOTAL DE TODOS LOS GAMMAS
 % =========================================================================
 %% PRECISION TOTAL (SEN + ESP + TPV +FPV)
 % figure
 % plot(1:1:size(sensibilidad_T,2),(sensibilidad_T+especificidad_T+Positive_predictive_value_T+False_predictive_value_T),'b','LineWidth',4);
 % hold on
 % plot(1:1:size(sensibilidad_T,2),(sensibilidad_T+especificidad_T+Positive_predictive_value_T+False_predictive_value_T),'go','MarkerSize',4);
 % for umbral = 1:size(sensibilidad_T,2)
 % plot((umbral),(sensibilidad_T(umbral)+especificidad_T(umbral)),'.b');
 % plot((umbral),(Positive_predictive_value_T(umbral)+False_predictive_value_T(umbral)),'og');
 % end
 % legend('(sensibilidad + especificidad + PPV + FPV)','(sensibilidad + especificidad + PPV + FPV)','(sensibilidad + especificidad)','(Positive predictive value + False predictive value)');
 % xlabel('Umbrales'); grid on; title ('Funcion de la presicion [0 = min ;medium = 2; 4 = max] para todos los gammas')
 % ylabel('(Sen + Esp + PPV + FPV )')
 % FIN DE PRECISION TOTAL
 % =========================================================================
 %%
 figure
 for umbral = 1:size(sensibilidad,1)
  %plot(1-especificidad(umbral,gamma),Positive_predictive_value(umbral,gamma),'x')
  plot(umbral,sal_s.ef_xgamma{umbral}/100,'.');hold on
  plot(umbral,sal_s.ef_xgamma{umbral}/100,'o') ;hold on
  % display(['Umbral ' num2str(sal_s.ef_xgamma{umbral})])
 end
 %end
 
 title ('Analisis x CADA gamma;')
 xlim([0 umbral]);ylim([0 1]);xlabel('umbrales ') ;ylabel(['% de sismos detectados']);
 
 %%
 if or(band == 1,band == 0)
  set(figure(1),'Position',[ 500   367   868   318])
  set(figure(2),'Position',[ 502    53   440   230])
  set(figure(3),'Position',[ -2    42   503   299])
  set(figure(5),'Position',[942    43   423   240])
  set(figure(4),'Position',[1   425   500   259])
 end
end % FIN


