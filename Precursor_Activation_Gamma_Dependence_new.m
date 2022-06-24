% Precursor_Activation_Gamma_Dependence.m
% ====================
% DESCRIPCIÓN: Programa toma las salidas del programa Optimal_Time_Interval_Produce .m
% y visualiza los resultados.
% SALIDA: Figuras y variable IA_MEAN_YS que contiene el precursor
% relativista
%=============================
close all
% =========================================================================
Mag_max = 10;
nd = 20*10.0;  % NUMERO DE DIAS BEFORE EVENT
d_gamma = 10;  % INCREMENTO DE GAMMA (PARA VISUALIZAR
so = 0; % CUANTAS DESVIACIONES SE TOMAN EN CUENTA PARA CONSIDERAR LOS MAXIMOS DEL PRECURSOR 
% -------------------------------------------------------------------------
load Precursor_Activation_Parameters % CONTIENE LAS ENTRADAS NECESARIAS
% ========================================================
%  band = 0; % band = 1 VISUALIZA TODAS LAS FIGURAS
% ========================================================
linas_anuales = 0; % 1,0 Visualiza lineas cada 365 dias en los datos
% ========================================================
anota = 0; %1,0 Visualiza cuadros de texto de calibracion y produccion
% ========================================================
% Ajusta la altura total del precursor para mejorar la visualizacion
k2 = 5; % igual a 5 si se quiere el default
Lxx = Lx;
Lyy = Ly;
% ----------------------------------------------------------
band = 1;
% CASOS
% -----------------------------------
switch Data_Stations
 case 1
  %CASO JAPON:
  ind = find(Lyy(:,2)>=9.1);
  Event_Lx =  Lxx(ind,:)-f1;
  Event_Ly =  Lyy(ind,:);
  % -----------------------------------
  % INTERVALO
  % ------------------
  n1 = 1;
  n2 = 5;
  % ------------------
 case 2
  %f1=datestr(round(min(Lx)));
  %f1=(round((Lx(1))));
  %     CASO EUROPA SILENCIO:
  ind2 = find(Lyy(:,2)>=6.5);
  ind = ind2(5);
  Event_Lx =  Lxx(ind,:)-f1;
  Event_Ly =  Lyy(ind,:);
  % ------------------
  n1 = 1;
  n2 = 4;
  % ------------------
  nd = 10.0;
  k2 = 5; % igual a 5 si se quiere el default
 case 3
  % TERREMOTO DE USA
  % -----------------------------------
  %CASO USA:
  ind2 = find(Lyy(:,2)>=7.2);
  ind = ind2(1);
  Event_Lx =  Lxx(ind,:)-f1;
  Event_Ly =  Lyy(ind,:);
  % nd = 30*5;                 % NUMERO DE DIAS BEFORE EVENT
  %CASO USA:
  %        ind2 = find(Lyy(:,2)==6.8);
  %        ind = ind2(1);
  %   Event_Lx =  Lxx(ind,:)-f1;
  %   Event_Ly =  Lyy(ind,:);
  % ------------------
  n1 = 1;
  n2 = 5;
  % ------------------
  nd = 20*10.0;
 case 4
  %   CASO CHILE:
  ind2 = find(Lyy(:,2)==8.0);
  ind = ind2(1);
  Event_Lx =  Lxx(ind,:)-f1;
  Event_Ly =  Lyy(ind,:);
  % ------------------
  n1 = 1;
  n2 = 10;
  % ------------------
 case 5
  %   CASO EUROPA DUR:
  ind2 = find(Lyy(:,2)>=5.1);
  ind = ind2(10);
  Event_Lx =  Lxx(ind,:)-f1;
  Event_Ly =  Lyy(ind,:);
  % -----------------------------------
  % ------------------
  n1 = 1;
  n2 = 5;
  % ------------------
  %  nd = 30*5.0;
 case 6
  %     CASO EUROPA PEG:
  ind2 = find(Lyy(:,2)==6.8);
  ind = ind2(1);
  Event_Lx =  Lxx(ind,:)-f1;
  Event_Ly =  Lyy(ind,:);
  % ------------------
  n1 = 1;
  n2 = 10;
  % ------------------
  nd = 20*10.0;
  %
  % %     CASO EUROPA PEG:
  %        ind2 = find(Lyy(:,2)==6.8);
  %        ind = ind2(1);
  %   Event_Lx =  Lxx(ind,:)-f1;
  %   Event_Ly =  Lyy(ind,:);
  %    % ------------------
  %      n1 = 1;
  %      n2 = 10;
  %   % ------------------
  %   nd = 20*10.0;
 case 7
  %     CASO JAPON PET:
  ind2 = find(Lyy(:,2)>=7.1);
  ind = ind2(2);
  Event_Lx =  Lxx(ind,:)-f1;
  Event_Ly =  Lyy(ind,:);
  % ------------------
  n1 = 1;
  n2 = 10;
  % ------------------
  nd = 21*10.0;
  %k2 = 14; % igual a 5 si se quiere el default
  %k2 = 14; % igual a 5 si se quiere el default
 case 8
  %CASO JAPON:
  ind2 = find(Lyy(:,2)>=6);
  ind = ind2(7);
  Event_Lx =  Lxx(ind,:)-f1;
  Event_Ly =  Lyy(ind,:);
  % -----------------------------------
  % INTERVALO
  % ------------------
  n1 = 1;
  n2 = 5;
  % ------------------
  k2 = 10;
end
% factor de suavizacion sobre la variable IA_MEAN_YS
ho = 5;
IA_MEAN_YS = [];
% =========================================================================
for i=1:size(IA_MEAN_Y,1),
 IA_MEAN_YS = [IA_MEAN_YS;smooth(IA_MEAN_Y(i,:),ho)'];
end
IM = mean(IA_MEAN_YS(n1:n2,:)/k2);
IA_M = mean(IA_MEAN_YS);
IA_M = -Mag_max*IA_M/max(IA_M);
IA_M = -(Mag_max+IA_M);
Days_Before = Event_Lx(1,1) - nd;
Days_After = Event_Lx(1,1) + nd ;
%Days_After = Event_Lx(1,1) + nd - 100;
Total_Days = round((Days_After-Days_Before)/2);
vx = round(Days_Before/Md):round(Days_After/Md); % INTERVALO TEMPORAL
max_tx = mean(std(IA_MEAN_YS(:,vx)'))*so;
[ty,tx] = max(IA_MEAN_YS(:,vx)'); % LOCALIZACION DE LOS MAXIMOS
ind_ty = find(ty>=max_tx); % ELIMINACION DE PRECURSORES INFERIORES A SU DESVIACION
tx = tx(ind_ty);
ty = ty(ind_ty);
Gamma_2 = Gamma_2(ind_ty);
IA_MEAN_YS = IA_MEAN_YS(ind_ty,:);
Tx = T(vx);                                      % POSICION DE LOS MAXIMOS
Dt = (Tx(tx) - Event_Lx(1))'; % DAYS BEFORE EVENT
[gy,gx] = max(ty);
Gamma_Opt = Gamma_2(gx);  % VALOR OPTIMO DE GAMMA
Gx = Gamma_Opt*ones(1,10);
Gy = linspace(0,1.1*gy,10);
[g_opt,E_opt] = Energia_Vs_Gamma(Gamma_Opt,0);
% =========================================================================
if band == 1
 % ------------------------------------------------------------------------
 set(figure(1),'Position',[112         212        1148         470],'Color','W')
 % ------------------------------------------------------------------------
 for i=1:size(IA_MEAN_YS,1)
  plot((Lx'-Lxx(1,1)),Ly','-b','LineWidth',[2]), hold on
  plot((Lx_Cont'-Lxx(1,1)),Ly_Cont','--g','LineWidth',[3])
  plot(Event_Lx,Event_Ly,'-r','LineWidth',[3])
  %Plot el precursor
  plot(T,IA_MEAN_YS(1:i,:)/k2,'-','LineWidth',[2])
  plot((Lx_Cont'-Lxx(1,1)),Ly_Cont','--g','LineWidth',[3])
  plot(Event_Lx,Event_Ly,'-r','LineWidth',[3])
  hold off
  xlabel(['Fecha inicial y final: ' datestr(round(f1)) ' --- hasta --- ' ...
   datestr(round(f2)) ' --- Tiempo total: [\approx' ...
   num2str(round((T(end)/365.25))) '] años --- [eje x en días]'])
  ylabel('Eventos (azul y verde) y Precursor (Curvas)')
  title(['Señal Precursora, ' char(Nombre_caso(Data_Stations))...
   ', Intervalos Gamma: [' num2str(min(Gamma_2),2) ':' ...
   num2str(max(Gamma_2)/length(Gamma_2),3) ':'...
   num2str(max(Gamma_2),3) ']'])
  xlim([0 T(end)+15])
  grid
  Ejes_Visibles(1)
  g = get(figure(1));
  g1 = g.Children;
  set(g1(1),'Position',[0.0345 0.1241 0.9612 0.8237],'Box','on');
  ax1 = gca;                   % current axes
  ax1.XTick = [0:200:T(end)];
  ax1.GridLineStyle = '--';
  ax1.GridColor = [0 0 0];
  ax1.GridAlpha = 1;
  ax1.LineWidth = 0.125;
 end
%% Cuadritos de texto de calibracion y produccion
 if anota == 1
  figure(1)
  % Pos de los 2 cuadros de texto:
  dim1 = [0.0349    0.8949    0.3217    0.0537]; %Esto se puede editar
  dim2 = [0.4394    0.8932    0.5558    0.0531]; %Esto se puede editar
  %Cuadro de texto para la calibracion
  anot1=annotation('textbox', dim1, 'String', 'Datos sísmicos de calibración',...
   'BackgroundColor','green','FontSize',12,'FaceAlpha',0.5,'Horizontalalignment','center',...
   'FontWeight','bold');
  % Cuadro de texto para la produccion
  anot2=annotation('textbox', dim2, 'String', 'Datos sísmicos de producción',...
   'BackgroundColor','blue','FontSize',12,'FaceAlpha',0.5,'Horizontalalignment','center',...
   'FontWeight','bold');
  % En caso se quiera modificar la pos a mano, luego se puede cambiar anot1 y
  % 2 con estos comandos:
  % get(anot1,'Position')
  % get(anot2,'Position')
 end
 
 %%
 % ------------------------------------------------------------------------
 set(figure(2),'Position',[5 382 1885 556],'Color','W')
 % ------------------------------------------------------------------------
 plot(T,IA_MEAN_YS/k2,'LineWidth',[2])
 hold on
 plot((Lx'-Lxx(1,1)),Ly','-b','LineWidth',[2])
 plot((Lx_Cont'-Lxx(1,1)),Ly_Cont','--g','LineWidth',[3])
 plot(Event_Lx,Event_Ly,'-r','LineWidth',[3])
 plot(T,IM,'-k','LineWidth',[8])
 plot(T,IM,'-r','LineWidth',[2])
 hold off
 xlabel(['Fecha inicial y final: ' datestr(round(f1)) ' --- hasta --- ' ... 
  datestr(round(f2)) ' --- Tiempo total: [\approx' ...
  num2str(round((T(end)/365.25))) '] años --- [eje x en días]'])
 ylabel('Eventos (azul y verde) y Precursor (curvas)')
 title(['Gamma: [' num2str(Gamma_2(i)) ']'])
 xlim([0 T(end)+10])
 grid
 g = get(figure(2));
 g1 = g.Children;
 set(g1,'Position',[0.0345 0.1241 0.9612 0.8237],'Box','on');
 ax1 = gca;                   % current axes
 ax1.XTick = [0:200:T(end)];
 ax1.GridLineStyle = '--';
 ax1.GridColor = [0 0 0];
 ax1.GridAlpha = 1;
 ax1.LineWidth = 0.125;
 Ejes_Visibles(2)
 if size(IA_MEAN_YS,1) == 10
  % CASO 10 GAMMAS
  % -----------------------------------------------------------------------
  set(figure(3),'Position',[1 28 1893 910],'Color','W')
  % -----------------------------------------------------------------------
  for i=1:10
   subplot(5,5,i)
   plot((Lx'-Lxx(1,1)),Ly','-b','LineWidth',[2]) ,hold on
   plot((Lx_Cont'-Lxx(1,1)),Ly_Cont','--g','LineWidth',[2])
   plot(Event_Lx,Event_Ly,'-r','LineWidth',[3])
   plot(T,IA_MEAN_YS(i,:)/k2,'-k','LineWidth',[5])
   plot(T,IA_MEAN_YS(i,:)/k2,'-r','LineWidth',[2])
   plot(T,IA_MEAN_YS(i,:)/k2,'-m','LineWidth',[1])
   hold off
   grid
   ylim([0 22])
   xlim([0 T(end)])
   title(['\gamma =' num2str(Gamma_2(i))],'FontSize',[12])
   %             axis off
   %             if i==9 | i==10,axis on,end
  end
  load Position_Fig_3
  g = get(figure(3));
  g1 = g.Children;
  for i=1:10,set(g1(i),'Position',[G3(i,:)],'Box','on');
   if i==9 | i==10,%set(g1(i),'Position',[G3(i,:)],'Box','on');
    xlabel(['Fecha inicial y final: ' datestr(round(f1)) ' --- hasta --- ' ...
     datestr(round(f2)) ' --- Tiempo total: [\approx' ...
     num2str(round((T(end)/365.25))) '] años --- [eje x en días]'])
   end
  end
  set(g1(1).XLabel,'Position',[-190.3790 -9.7333 0])
  Ejes_Visibles(3)
 end
 if size(IA_MEAN_YS,1) <= 20
  % CASO 20 GAMMAS
  % -----------------------------------------------------------------------
  set(figure(4),'Position',[1 28 1893 910],'Color','W')
  % -----------------------------------------------------------------------
  for i=1:size(IA_MEAN_YS,1)
   ym = max(IA_MEAN_YS(:))/k2;
   subplot(10,2,i)
   plot((Lx'-Lxx(1,1)),Ly','-b','LineWidth',[2]) ,hold on
   plot((Lx_Cont'-Lxx(1,1)),Ly_Cont','--g','LineWidth',[2])
   plot(Event_Lx,Event_Ly,'-r','LineWidth',[3])
   plot(T,IA_MEAN_YS(i,:)/k2,'-k','LineWidth',[5])
   plot(T,IA_MEAN_YS(i,:)/k2,'-r','LineWidth',[2])
   plot(T,IA_MEAN_YS(i,:)/k2,'-m','LineWidth',[1])
   hold off
   grid
   ylim([0 ym])
   xlim([0 T(end)])
   title(['\gamma =' num2str(Gamma_2(i))],'FontSize',[12])
   %axis off
   %if i==9 | i==10,axis on,end
  end
  load  Position_Fig_3_20
  g = get(figure(4));
  g1 = g.Children;
  for i=1:size(IA_MEAN_YS,1),set(g1(i),'Position',[G3(i,:)],'Box','off');
   set(g1(i),'Position',[G3(i,:)],'Box','on');
   if i==18,
    xlabel(['Fecha inicial y final: ' datestr(round(f1)) ' --- hasta --- '...
     datestr(round(f2)) ' --- Tiempo total: [\approx' ...
     num2str(round((T(end)/365.25))) '] años --- [eje x en días]'])
   end
  end
  set(g1(1).XLabel,'Position',[-6.9672e+03 -7.0440 0])
  Ejes_Visibles(4)
 end
 
 if size(size(IA_MEAN_X),1)>1
  % -----------------------------------------------------------------------
  set(figure(5),'Position',[5 382 1885 556],'Color','W')
  % -----------------------------------------------------------------------
  plot(smooth(mean(IA_MEAN_X/k2),4),'-k','LineWidth',[6]),hold on
  plot(smooth(mean(IA_MEAN_X/k2),4),'-g','LineWidth',[3])
  plot(smooth(mean(IA_MEAN_X/k2),4),'-r','LineWidth',[1]),hold off
  grid
  Ejes_Visibles(5)
  g = get(5);
  g1 = g.Children;
  set(g1(1),'Position',[0.0473 0.1100 0.9278 0.8150])
  
  % DETERMINACION DE LOS PICOS
  % ------------------------------------
  Xa = zeros(size(IA_MEAN_YS));
  Ya = zeros(size(IA_MEAN_YS));
  for i=1:size(IA_MEAN_YS,1),
   [ya,xa] = findpeaks(IA_MEAN_YS(i,:));
   plot(T,IA_MEAN_YS(1:i,:)),hold on, plot(T(xa),ya,'*r'),pause
   Xa(i,1:length(xa)) = T(xa);
   Ya(i,1:length(xa)) = ya;
  end
 end
 % ------------------------------------------------------------------------
 set(figure(6),'Position',[6 32 1883 906],'Color','W')
 % ------------------------------------------------------------------------
 lx = T;
 
 hold on,
 for i=1:size(IA_MEAN_YS,1)
  ly = (20*i+20)*ones(1,size(IA_MEAN_YS,2));
  plot(lx,ly,'--k')
  plot(T,IA_MEAN_YS(i,:)'+20*i+20,'-','LineWidth',[1])
  text(lx(end),ly(end),['\gamma = ' num2str(round(Gamma_2(i)))],'Color','k','FontWeight','Bold','FontSize',[12])
 end
 plot((Lx'-Lxx(1,1)),70*Ly','-b','LineWidth',[2])
 plot((Lx_Cont'-Lxx(1,1)),70*Ly_Cont','--g','LineWidth',[2])
 plot(Event_Lx,70*Event_Ly,'-r','LineWidth',[3])
 hold off
 xlabel(['Fecha inicial y final: ' datestr(round(f1)) ' --- hasta --- ' datestr(round(f2)) ' --- Tiempo total: [\approx' num2str(round((T(end)/365.25))) '] años --- [eje x en días]'])
 g = get(figure(6));
 g1 = g.Children;
 set(g1,'Position',[0.0398 0.1100 0.9092 0.8337])
 set(g1(1).XLabel,'Position',[2.0040e+03 -21.9426 -1])
 xlim([T(1) T(end)+10])
 grid
 Ejes_Visibles(6)
 % ========================================================================
 set(figure(7),'Position',[5 385 1885 553],'Color','W')
 % --------------------------------------------------------------
 
 if size(IA_MEAN_YS,1) < 2
  set(figure(7),'Position',[5 362 1885 576],'Color','W')
  plot((Lx'-Lxx(1,1)),Ly'-Mag_max,'-b','LineWidth',[2]),hold on
  plot((Lx_Cont'-Lxx(1,1)),Ly_Cont'-Mag_max,'--g','LineWidth',[2])
  plot(Event_Lx,Event_Ly-Mag_max,'-r','LineWidth',[3])
  plot(Event_Lx, 1.2*max(IA_MEAN_YS(:))*(Event_Ly),'-k','LineWidth',[4])
  plot(Event_Lx, 1.2*max(IA_MEAN_YS(:))*(Event_Ly),'--r','LineWidth',[2])
  plot(T,IA_MEAN_YS,'-k','LineWidth',[10]),hold on
  plot(T,IA_MEAN_YS,'-r','LineWidth',[5])
  plot(T,IA_MEAN_YS,'-g','LineWidth',[2])
  hold off
  title(['Evolución temporal del precursor \gamma = [' num2str(Gamma_2) ']' ] )
  ylabel('Amplitud del precursor [au]','FontSize',[14])
  ylim([-Mag_max 1.2*max(IA_MEAN_YS(:))])
  g = get(7);
  g1 = g.Children;
  set(g1(1),'Position',[0.0377 0.1100 0.9480 0.8150],'Box','on')
 else
  [xa,ya] = meshgrid(linspace(T(1),T(end),size(IA_MEAN_YS,2)),linspace(Gamma_2(1),Gamma_2(end),size(IA_MEAN_YS,1)));
  surf(xa,ya,IA_MEAN_YS),shading interp,view([0 90])
  colorbar,colormap jet,alpha 0.5
  axis xy,axis on
  hold on,plot((Lx'-Lxx(1,1)),Ly'-Mag_max,'-b','LineWidth',[2])
  plot((Lx_Cont'-Lxx(1,1)),Ly_Cont'-Mag_max,'--g','LineWidth',[2])
  plot(Event_Lx,Event_Ly-Mag_max,'-r','LineWidth',[3])
  plot(Event_Lx, 1.2*max(IA_MEAN_YS(:))*(Event_Ly),'-k','LineWidth',[4])
  plot(Event_Lx, 1.2*max(IA_MEAN_YS(:))*(Event_Ly),'--r','LineWidth',[2])
  plot(T,IA_M,'-k','LineWidth',[8])
  plot(T,IA_M,'-r','LineWidth',[3])
  ylim([-Mag_max Gamma_2(end)])
  hold off
  ylabel('valores de \gamma  ','FontSize',[14])
  title('Evolución en el tiempo del precursor y sus dependencias con el valor gamma')
  %title('Precursor time evolution and its gamma values dependences')
  ylim([-10 Gamma_2(end)])
  g = get(7);
  g1 = g.Children;
  set(g1(1),'Position',[0.9620    0.1103    0.0141    0.8156],'Box','on')
  set(g1(2),'Position',[0.0377    0.1465    0.9093    0.7785],'Box','on')
  set(g1(2).XLabel,'Position',[1.9943e+03 -13.1462 0])
  set(g1(2).YLabel,'Position',[-97.6563 10.3596 0])
 end
 xlabel(['Fecha inicial y final: ' datestr(round(f1)) ' --- hasta --- ' datestr(round(f2)) ' --- Tiempo total: [\approx' num2str(round((T(end)/365.25))) '] años --- [eje x en días]'])
 xlim([T(1) T(end)])
 ax1 = gca;                   % current axes
 ax1.XTick = [0:200:T(end)];
 ax1.GridLineStyle = '--';
 ax1.GridColor = [0 0 0];
 ax1.GridAlpha = 1;
 ax1.LineWidth = 0.125;
 Ejes_Visibles(7)
 
 if size(IA_MEAN_YS,1)>1
  % -----------------------------------------------------------------------
  set(figure(8),'Position',[612 32 1278 906],'Color','W')
  % -----------------------------------------------------------------------
  subplot(1,3,1)
  plot(Tx-Tx(1),IA_MEAN_YS(:,vx)','LineWidth',[2])
  hold on
  %plot(Tx-Tx(1),mean(IA_MEAN_YS(25:35,vx)),'-k','LineWidth',[2])
  plot((Lx'-Lxx(1,1))-Tx(1),10*Ly','-b','LineWidth',[0.125])
  plot(Tx(tx)-Tx(1),ty,'*k','MarkerSize',[8])
  plot(Tx(tx)-Tx(1),ty,'.g','MarkerSize',[8])
  plot(Event_Lx-Tx(1),10*Event_Ly,'-r','LineWidth',[3])
  grid
  hold off
  title('Evento sismico aislado para varios valores de \gamma ')
  xlabel('Tiempo (días)')
  ylabel('Amplitud del precursor (ua)')
  xlim([0 Tx(end)-Tx(1)])
  ax1 = gca;                   % current axes
  ax1.XTick = [0:500:Tx(end)-Tx(1)];
  ax1.GridLineStyle = '--';
  ax1.GridColor = [0 0 0];
  ax1.GridAlpha = 1;
  ax1.LineWidth = 0.125;
  subplot(1,3,2)
  [yd,xd] = hist(Dt,200);
  bd = bar(xd,yd);
  bx = zeros(1,10);
  by = linspace(0,max(ya(:)),10);
  hold on,plot(bx,by,'-b','LineWidth',[2]),hold off
  set(bd,'EdgeColor',[1 0 0],'BarWidth',[2],'LineWidth',[2],'LineStyle','-')
  grid
  xlabel('Días antes (<0) días después (>0)','FontSize',[10])
  ylabel('Contador')
  ax1 = gca;                   % current axes
  ax1.XTick = [-nd:20:nd];
  ax1.GridLineStyle = '--';
  ax1.GridColor = [0 0 0];
  ax1.GridAlpha = 1;
  ax1.LineWidth = 0.125;
  subplot(1,3,3)
  plot(Gamma_2,Dt,'-k','LineWidth',[8]),hold on
  plot(Gamma_2,Dt,'-r','LineWidth',[4])
  plot(Gamma_2,Dt,'oc','LineWidth',[4])
  plot(Gamma_2,Dt,'.m','MarkerSize',[14])
  hold off
  g = get(8);
  g1 = g.Children;
  set(g1(1),'Position',[0.4844 0.0640 0.5007 0.4636],'Box','on')
  set(g1(2),'Position',[0.4844 0.5982 0.5003 0.3620],'Box','on')
  set(g1(3),'Position',[0.0382 0.0618 0.3966 0.8985],'Box','on')
  ax1 = gca;                   % current axes
  ax1.XTick = [0:5:Gamma_2(end)];
  %ax1.YTick = [-nd:100:nd];
  ax1.GridLineStyle = '--';
  ax1.GridColor = [0 0 0];
  ax1.GridAlpha = 1;
  ax1.LineWidth = 0.125;
  ylim([-1.2*max(abs(Dt)) 1.2*max(abs(Dt))])
  xlim([0 Gamma_2(end)])
  xlabel('Valores \gamma ')
  ylabel('Días antes (<0) --- días después (>0)','FontSize',[12])
  grid
  Ejes_Visibles(8)
  % -------------------------------------------------------------------------------
  set(figure(9),'Position',[5 385 1885 553],'Color','W')
  plot(Gamma_2,ty,'-b','LineWidth',[2]),grid,hold on
  plot(Gamma_2,smooth(smooth(ty,6),3),'-r','LineWidth',[6])
  plot(Gamma_2,smooth(smooth(ty,20),10),'-g','LineWidth',[6])
  plot(Gamma_2,smooth(smooth(ty,20),10),'-k','LineWidth',[2])
  plot(Gx,Gy,'--m','LineWidth',[2])
  hold off
  xlim([Gamma_2(1) Gamma_2(end)])
  xlabel('Valores \gamma ')
  ylabel('Amplitud del precursor (ua)')
  title(['Amplitud del precursor [' num2str(Total_Days) '] dias antes/despues al evento --- Gamma-Opt: [' num2str(round(Gamma_Opt)) '] Opt-Energía: [' num2str(E_opt) '] GeV'])
  g = get(9);
  g1 = g.Children;
  set(g1,'Position',[0.0419 0.1100 0.9464 0.8150],'Box','on')
  ax1 = gca;                   % current axes
  ax1.XTick = [0:3:Gamma_2(end)];
  ax1.GridLineStyle = '--';
  ax1.GridColor = [0 0 0];
  ax1.GridAlpha = 1;
  ax1.LineWidth = 0.125;
  Ejes_Visibles(9)
 end
end
% OPTIMAL GAMMA VALUES
% ------------------------
disp('-------------------------------------------------------------------')
disp('Valores Optimos de Gamma: ')
disp(Gamma_2')
disp('-------------------------------------------------------------------')
disp(['Mean Optimal Gamma Value: ' num2str(mean(Gamma_2))])
disp(['Std Optimal Gamma Value: ' num2str(std(Gamma_2))])
% -------------------------------------------------------------------------
close(figure(3))
close(figure(4))
%% Muestra la cantidad de años en el grafico con lineas verticales naranjas
if linas_anuales == 1
 figure(1)
 hold on;
 TOO = 365.25;
 for i = 1:round((f2-f1)/365.25)
  if TOO+365.25 <= round(f2-f1+365)
   plot([TOO, TOO], [0, 10],'Color',[1 0.5 0],'LineWidth',1.3,'LineStyle','-'),hold on;
   TOO=TOO+365.25;
  end
 end
 switch Data_Stations
  case {7,4,2} %para el caso de 4, los datos no llegan hasta el sismo
   ts = datenum([2011 3 11]) - f1; %terremoto 9.1
   for i =1:2
    figure(i)
    hold on
    plot([ts, ts], [0, 9.1],'Color',[1 0 0],'LineWidth',4,'LineStyle','-');
    plot([ts, ts], [0, 9.1],'Color',[1 0 1],'LineWidth',4,'LineStyle','--'),hold on;
   end
 end
end
%%
% FIN DEL PROGRAMA



