% ESTACIONES INTERMAGNET
% ====================
% Entradas:
% Estaciones_intermagnet:
% [Contiene la localización de las estaciones INTERMAGNET]
% IEBQuakesExport31767_m5:
% [Contiene la informacion relacionada a los sismos]
% DESCRIPCIÓN: Muestra distintas representaciones graficas de los archivos:
% Con eso da la opción de filtrar por radio, por magnitud, por año inicial
% y final, por mes y por caso, tambien da la opción de:
% caso = 1 se abre un mapa mundial y se centra en un terremoto cercano 
% al click, busca 3 estaciones cercanas y te dice la distancia 
% terremoto- estaciones.
% caso = 2 Al hacer click en el mapa, se centra en la estación más cercana 
% al click y marca todos los terremotos dentro del radio seleccionado. 
% caso = 3 Muestra un mapa general mundial y sus sismos.
% caso = 4 Dando el nombre directo de 3 estaciones INTERMAGNET, 
% permite manejar los radios individuales de cada estación y muestra un 
% mapa de los sismos cercanos además de generar las tablas de salida que 
% pueden ser usados posteriormente. 
% SALIDA: T % Tabla con datos sismicos de la zona.
%=============================
clear all
close all
%% Filtrado de entradas
%Entradas
%==========================================================================
radio=300;                       % En kilometros
Mo=6;                            % Filtro de magnitud minima
Depth_max =1000;                 % PROFUNDIDAD MAXIMA EN KM
T1 = 2004;                       % TIEMPO INICIAL DE LA DATA
T2 = 2015;                       % TIEMPO FINAL DE LA DATA
M1 = 1;                          % MES INFERIOR
M2 = 13;                         % MES SUPERIOR
caso=4;
% ======================================================================

Time1 = T1; 
Time2 = T2;

if caso == 4
r1=500;%radio estacion 1
r2=500;%radio estacion 2
r3=500;%radio estacion 3
display('Caso 4. Estaciones:')
Label_Station={'KAK','KNY','MMB'}
end
% generando un area en el cual se tomarán los sismos dentro de ese area
% solamente
% =========================================================================
imprimir=0; %0 = no guarda imagenes ; 1= si guarda
if imprimir==1
     fileName='MapaTerremotos_Estaciones_Generado'; %MapaMundiGuardado
     fileNameC1='ZonaCaso1';
     fileNameC2='ZonaCaso2';
     save Etiquetas_de_losSismos_vs_Tiempo_Japon2 YearS MonthS DayS 
 end
% ==========================================================================
% Carga las entradas necesarias para procesar:
load ('Estaciones_intermagnet');load('IEBQuakesExport31767_m5'); 
latitud=90-Colatitute  ;        %90-colatitud=latitud
Region_Total = Region;
Long_Intermagnet=EastLongitude;
ind = find(Long_Intermagnet>=180);
Long_Intermagnet(ind) = Long_Intermagnet(ind) - 360;
Long_Intermagnet_New = Long_Intermagnet;
% ==========================================================================
ind_m = [];
% -------------------------------------------------
% LAZO PARA LOCALIDADES Y MAGNITUDES SELECCIONADAS
% -------------------------------------------------
ind_o = find(Mag >= Mo & Depthkm <= Depth_max & Year >= T1 & Year <= T2 ...
    & Month >= M1 & Month <= M2);
ind_m = [ind_m;ind_o];
Ind_M = ind_o;
Number_total_of_Events = length(Ind_M);

% SOLO LOS DATOS FILTRADOS
% ------------------------------------
Day  = Day(Ind_M);
Depthkm = Depthkm(Ind_M);
IRISID  = IRISID(Ind_M);
Lat  = Lat(Ind_M);
Lon  = Lon(Ind_M);
Mag = Mag(Ind_M);
Month = Month(Ind_M);
Region = Region(Ind_M);
TimeUTC = TimeUTC(Ind_M);
Timestamp = Timestamp(Ind_M);
Year = Year(Ind_M);
% ---------------------------------
%% MUESTRA LOS MAPAS PARA EMPEZAR
% Muestra las estaciones 
figure1 = figure(1);
axes1 = axes('Parent',figure1);
set(axes1,'ClippingStyle','rectangle','Color',...
    [0.8313725490 0.815686274509 0.78431372549],'DataAspectRatio',[1 1 1]);
% - - - - - - - -
% pinta estaciones:
geoshow('landareas.shp', 'FaceColor', [0 0 0]);hold on;  ylim([-90 85]);
xlim([-180 180]);grid
geoshow(latitud,Long_Intermagnet,'DisplayType','point','Marker', '.',...
    'MarkerEdgeColor','g','MarkerSize',17)
title('Longitude(x axes),Latitude(y axes), g=Stations');
title({['Posicion de estaciones Intermagnet, ' ...
    'Numero de terremotos mostrados: ' num2str(length(Mag)) ]});
xlabel('Longitude ');
ylabel('Latitude');
% pinta terremotos: de acuerdo a su magnitud:
color2=linspace(0, 1, length(Mag));% generates N points between X1 and X2.
%% GRAFICA TERREMOTOS:
for i=1:length(Mag)
     i=length(Mag)-i +1;
     geoshow(Lat(i),Lon(i),'DisplayType','point','Marker', '.',...
         'MarkerEdgeColor',[(0.18+(Mag(i)/10 - 0.09999)) color2(i) 0],...
         'MarkerSize',exp(Mag(i)/2.45))
end
% Imprime la figura:
if imprimir==1;
    ImprimeFigura(fileName)
end
% contra más rojo,mayor magnitud,se puede seguir ajustando
%% CASO 2:
% CASO 2:
for resumir2=1:1
    if caso==2
        close (figure(2))
        [xo,yo] = ginput(1); % toma un punto del mapa 
        [Numero_Estacion] = Busca_Estacion_cercana...
            (latitud,Long_Intermagnet,xo,yo);
        % ================================================
        % Una vez obtenido el sitio,procedemos a buscar terremotos 
        % en el radio especificado:
        ELON=Long_Intermagnet(Numero_Estacion);
        ELAT=latitud(Numero_Estacion);
        [T_Dentro_Radio_indices,Distancias_respectivas] = ...
            Terremotos_Radio(Lon,Lat,ELON,ELAT,radio);      
        %% VISUALIZA UN CIRCULO Y BORRA LO DEMÁS
        % Visualiza la region "filtrada" con un circulo
        %un circulito:
        theta = linspace(0,2*pi);
        r = radio/111.1;
        xc = Long_Intermagnet(Numero_Estacion);
        yc = latitud(Numero_Estacion);
    % Plot a circle centered at the point (xc,yc) with a radius equal to r.
        x = r*cos(theta) + xc;
        y = r*sin(theta) + yc;
        % ------- visualiza solo lo que este dentro del circulo:
        figure(2)
        figure1 = figure(2);
        axes1 = axes('Parent',figure1);
        set(axes1,'ClippingStyle','rectangle','Color',...
            [0.831372549019608 0.815686274509804 0.784313725490196],...
            'DataAspectRatio',[1 1 1]);
        geoshow('landareas.shp', 'FaceColor', [0 0 0]);hold on;
        xlim([(x(50)-1) (x(100)+1)]);ylim([y(75) y(25)]);
        grid
        geoshow(latitud,Long_Intermagnet,'DisplayType','point',...
            'Marker', 'x','LineWidth',10,'MarkerEdgeColor',...
            [1.0000    0.8431         0],'MarkerSize',24)
        geoshow(latitud(Numero_Estacion), ...
            Long_Intermagnet(Numero_Estacion),'DisplayType','point',...
            'Marker', '.','LineWidth',3,'MarkerEdgeColor','g',...
            'MarkerSize',12)
        text(Long_Intermagnet,latitud+1.15,IAGA,'Color',[ 0 0.4980 0],...
            'FontSize',13,'HorizontalAlignment','Center',...
            'Backgroundcolor',[0.9412 0.9412 0.9412])
        title(['Geographical position and events, central Station:' ...
            char(IAGA(Numero_Estacion)) ', Date:'  num2str(Time1) ' to '...
            num2str(Time2)])
        %- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        DayS=Day(T_Dentro_Radio_indices);
        MonthS=Month(T_Dentro_Radio_indices);
        RegionS=Region(T_Dentro_Radio_indices);
        TimeUTCS=TimeUTC(T_Dentro_Radio_indices);
        TimestampS=Timestamp(T_Dentro_Radio_indices);
        YearS = Year(T_Dentro_Radio_indices);
        DepthkmS=Depthkm(T_Dentro_Radio_indices);
        IRISIDS=IRISID(T_Dentro_Radio_indices);
        LatS=Lat(T_Dentro_Radio_indices);
        LonS=Lon(T_Dentro_Radio_indices);
        MagS=Mag(T_Dentro_Radio_indices);
        % - - - - - - -  - Plot terremotos:
        %% Plot terremotos
        for i=1:length(MagS)
           i=length(MagS)-i +1;
           geoshow(LatS(i),LonS(i),'DisplayType','point','Marker', '.',...
             'MarkerEdgeColor',[(0.29+(MagS(i)/10 - 0.09)) color2(i) 0],...
             'MarkerSize',(MagS(i)^2.2))
           geoshow(LatS(i),LonS(i),'DisplayType','point','Marker', '.',...
           'MarkerEdgeColor',[(0.29+(MagS(i)/10 - 0.09)) 0.5 color2(i)],...
            'MarkerSize',(MagS(i)^1.8))
        end
        plot(x,y,'--','LineWidth',3)
        % TABLA DISTANCIAS FINAL:
        if length(Distancias_respectivas)>=1 %Evita un error
            T=table(MagS,round(Distancias_respectivas),LatS,LonS...
                ,DayS,MonthS,YearS,TimeUTCS,TimestampS,DepthkmS,IRISIDS,...
                'VariableNames',{'Magnitud' 'DistanciasAEstacionKm' ...
                'Latitud' 'Longitud' 'Dia' 'Mes' 'Year' 'TimeUTC' ...
                'Timestamp' 'DepthkmS' 'IRISIDS'});
        else
            T=table(num2str(MagS),round(Distancias_respectivas),...
                num2str(LatS),num2str(LonS),DayS,MonthS,YearS,TimeUTCS,...
                num2str(TimestampS),IRISIDS,DepthkmS,'VariableNames'...
            ,{'Magnitud' 'DistanciasAEstacionKm' 'Latitud' 'Longitud' ...
           'Dia' 'Mes' 'Year' 'TimeUTC' 'Timestamp' 'IRISIDS' 'DepthkmS'});
        end
        display(T)
        if imprimir==1; %Imprime la figura
            ImprimeFigura(fileNameC2)
        end
       Ejes_Visibles(2); xlabel(['Longitude (' char(176) ')']); 
       ylabel(['Latitude(' char(176) ')']);
    end
    Parametros_entrada.radio = radio;
    Parametros_entrada.Magnitud_Minima = Mo;
    Parametros_entrada.Profundidad_Max = Depth_max;
    Parametros_entrada.Year_Ini = T1;
    Parametros_entrada.Year_Final = T2;
end
% PARA GUARDAR LAS SALIDAS DESCOMENTAR ESTO:
if caso == 2;
  Dist_respectivas=round(Distancias_respectivas);
%  save ('Etiquetas_KNYEntrenamiento2','T','Parametros_entrada'...
%       ,'MagS','Dist_respectivas','LatS','LonS','DayS','MonthS',...
%       'YearS','TimeUTCS','TimestampS','IRISIDS','DepthkmS')
%  save ('Etiquetas_KNYProduccion2','T','Parametros_entrada'...
%      ,'MagS','Dist_respectivas','LatS','LonS','DayS','MonthS',...
%      'YearS','TimeUTCS','TimestampS','IRISIDS','DepthkmS')
end

%% CASO 1:
for resumir1=1:1
    if caso==1
    figure(1)
    [xo,yo]=ginput(1);
    % Busca el terremoto m?s cercano al click:
    [Numero_Terremoto]=Busca_Terremoto_cercano(Lat,Lon,xo,yo);
    LatS1=Lat(Numero_Terremoto);
    LonS1=Lon(Numero_Terremoto);
    ELON1=Long_Intermagnet;
    ELAT1=latitud;
    [E_Dentro_Radio_indices,Distancias_respectivasE] = ...
      Estaciones_Radio(LatS1,LonS1,ELON1,ELAT1,radio);
    % SALIDA:
    % TABLA DISTANCIAS FINAL:
    display({'Magnitud:' num2str(Mag(Numero_Terremoto))});
    if length(Distancias_respectivasE)>=1
  T2=table((IAGA(E_Dentro_Radio_indices),round(Distancias_respectivasE),...
    Long_Intermagnet(E_Dentro_Radio_indices),...
    ELAT1(E_Dentro_Radio_indices),Country(E_Dentro_Radio_indices),...
  'VariableNames',{'IAGA' 'Distancias_km' 'Longitud' 'Latitud' 'Country'});
    else
   T2=table(IAGA(E_Dentro_Radio_indices),round(Distancias_respectivasE),...
    num2str(Long_Intermagnet(E_Dentro_Radio_indices)),...
 num2str(ELAT1(E_Dentro_Radio_indices)),Country(E_Dentro_Radio_indices),...
  'VariableNames',{'IAGA' 'Distancias_km' 'Longitud' 'Latitud' 'Country'});
    end
    T2 = sortrows(T2,'Distancias_km','ascend');
    display(T2)
    % FIN DE LA SECCION DE TABLAS
    if imprimir==1; %Imprime la figura:
       ImprimeFigura(fileNameC1)
    end
    end
end
% Las salidas de ambos casos son las tablas T Y T2 para el caso 2 y 1
% respectivamente. 
%% CASO 3
% Este caso simplemente es otra manera de ver los datos
% NO LO RECOMIENDO:
for resumir3=1:1
if caso==3 %Visualiza todo el mundo y pinta las estaciones con su nombre 
close (figure(3))
figure(3)
ax = worldmap('World');
land = shaperead('landareas', 'UseGeoCoords', true);
geoshow(ax, land, 'FaceColor', [0.5 0.7 0.5])
title('Longitude(x axes),Latitude(y axes)')
scatterm(latitud,EastLongitude,10,'filled')
scaleruler('units','km')
textm(latitud+1,EastLongitude+1,IAGA,'FontSize',10);
%------------------------------------------------------------------
scatterm(Lat,Lon,9,'filled','r');
plotm(latitud,EastLongitude,'xb','LineWidth',2)
% para cosas mas dinamicas,revisa: Text Properties Clipping
C = num2cell(Year);
D = num2cell([1:p]);
textm(Lat,Lon,C,'FontSize',7);%GET(H)
end
end
%% CASO 4
if caso == 4;
% BUSCA LAS POSICIONES DE LAS ESTACIONES:
idx = [];
Ind_station = [];
for j=1:length(Label_Station);
    for i=1:length(IAGA)
        idx = strfind(IAGA(i),Label_Station(j));
        if (cell2mat(idx)== 1)
            % Contiene el indice en la variable IAGA
            Ind_station = [Ind_station,i];
            i = length(IAGA);
        end
    end
end
clearvars j idx i
% Genera area sobre las estaciones:
r_todo = [r1 r2 r3]; % Radio 1, 2 y 3 ejem;('KNY','KAK','MMB')
ELON=Long_Intermagnet((Ind_station));
ELAT=latitud((Ind_station));
clearvars Eventos_f
EventF.Stat  = [(Label_Station(1)); (Label_Station(2)); (Label_Station(3))];
EventF.Stat = struct('Ind',[],'Dist',[]);
for std_f =1:length(Label_Station)
    display(['Buscando sismos en estacion:' Label_Station(std_f)])
    [T_Dentro_Radio_indices,Distancias_respectivas] = ...
        Terremotos_Radio(Lon,Lat,ELON(std_f),ELAT(std_f),r_todo(std_f));
    % Esta estructura contiene todos los datos filtrados de los sismos.
    EventF.Stat(std_f).Ind = T_Dentro_Radio_indices;
    EventF.Stat(std_f).Dist = Distancias_respectivas;
end
% Crea 3 circulos:
theta = linspace(0,2*pi);
r=[];
for std_f= 1: length(Label_Station)
    r = r_todo(std_f)/111.1;
    xc = Long_Intermagnet(Ind_station(std_f));
    yc = latitud(Ind_station(std_f));
    % Plot a circle centered at the point (xc,yc) with a radius equal to r.
    x(std_f,:) = r*cos(theta) + xc;
    y(std_f,:) = r*sin(theta) + yc;
end
% Encuentra los extremos de las estaciones para ajustar la figura:
Max_x= max(x(:));Max_y= max(y(:));
Min_x= min(x(:));Min_y= min(y(:));
% Prepara figura:
figure1 = figure(2);
axes1 = axes('Parent',figure1);
set(axes1,'ClippingStyle','rectangle','Color',...
    [0.831372549 0.815686274 0.78431372],'DataAspectRatio',[1 1 1]);
geoshow('landareas.shp', 'FaceColor', [0 0 0]);hold on;
xlim([(Min_x-2) (Max_x+2)]);
ylim([Min_y-2 Max_y+2]);
grid,geoshow(latitud(Ind_station),Long_Intermagnet(Ind_station),...
    'DisplayType','point','Marker', 'x','LineWidth',10,'MarkerEdgeColor'...
    ,[1.0000 0.8431 0],'MarkerSize',24)
geoshow(latitud(Ind_station),Long_Intermagnet(Ind_station),...
    'DisplayType','point','Marker', '.','LineWidth',3,'MarkerEdgeColor',...
    'g','MarkerSize',12)
text(Long_Intermagnet(Ind_station),latitud(Ind_station)+1.15,...
    IAGA(Ind_station),'Color',[ 0 0.4980 0],'FontSize',9,...
    'HorizontalAlignment','Center','Backgroundcolor',[0.9412 0.9412 0.941])
title('Region selected,r=Earthquake g=Stations')
for std_f= 1: length(Label_Station)
    DayS=Day(EventF.Stat(std_f).Ind);
    MonthS=Month(EventF.Stat(std_f).Ind);
    RegionS=Region(EventF.Stat(std_f).Ind);
    TimeUTCS=TimeUTC(EventF.Stat(std_f).Ind);
    TimestampS=Timestamp(EventF.Stat(std_f).Ind);
    YearS = Year(EventF.Stat(std_f).Ind);
    DepthkmS=Depthkm(EventF.Stat(std_f).Ind);
    IRISIDS=IRISID(EventF.Stat(std_f).Ind);
    LatS=Lat(EventF.Stat(std_f).Ind);
    LonS=Lon(EventF.Stat(std_f).Ind);
    MagS=Mag(EventF.Stat(std_f).Ind);
    % CREA tablas de salida de cada estacion:
    if  length(EventF.Stat(std_f).Dist)>=1
        switch std_f
            case 1,
                T1 = table(MagS,round(EventF.Stat(std_f).Dist),LatS,LonS...
                    ,DayS,MonthS,YearS,TimeUTCS,TimestampS,DepthkmS,...
                    IRISIDS,'VariableNames',{'Magnitud' ...
                    'DistanciasAEstacionKm' 'Latitud' 'Longitud' 'Dia'...
                    'Mes' 'Year' 'TimeUTC' 'Timestamp' 'DepthkmS'...
                    'IRISIDS'});
            case  2,
                T2 = table(MagS,round(EventF.Stat(std_f).Dist),LatS,LonS...
                    ,DayS,MonthS,YearS,TimeUTCS,TimestampS,DepthkmS,...
                    IRISIDS,'VariableNames',{'Magnitud' ...
                    'DistanciasAEstacionKm' 'Latitud' 'Longitud' 'Dia'...
                    'Mes' 'Year' 'TimeUTC' 'Timestamp' 'DepthkmS'...
                    'IRISIDS'});
            case  3,
                T3 = table(MagS,round(EventF.Stat(std_f).Dist),LatS,LonS...
                    ,DayS,MonthS,YearS,TimeUTCS,TimestampS,DepthkmS,...
                    IRISIDS,'VariableNames',{'Magnitud' ...
                    'DistanciasAEstacionKm' 'Latitud' 'Longitud' 'Dia'...
                    'Mes' 'Year' 'TimeUTC' 'Timestamp' 'DepthkmS'...
                    'IRISIDS'});
                Parametros_entrada.radio = [r1 r2 r3];
                Parametros_entrada.Magnitud_Minima = Mo;
                Parametros_entrada.Profundidad_Max = Depth_max;
                Parametros_entrada.Year_Ini = T1;
                Parametros_entrada.Year_Final = T2;
        end
        % fin de creacion de tablas T1,T2,T3 y Parametros_entrada
    end
    % - - - - - - -  - Plot terremotos:
    for i=1:length(MagS)
        hold on
        i=length(MagS)-i +1;
        geoshow(LatS(i),LonS(i),'DisplayType','point','Marker', '.',...
           'MarkerEdgeColor',[(0.18+(MagS(i)/10 - 0.09999)) color2(i) 0]...
           ,'MarkerSize',(MagS(i)^1.7))
    end
    plot(x(std_f,:),y(std_f,:),'--','LineWidth',1.5)
    xlabel('Longitude ');
    ylabel('Latitude');
end
P1 = [ Long_Intermagnet(Ind_station(1)) latitud(Ind_station(1))];
P2 = [ Long_Intermagnet(Ind_station(2)) latitud(Ind_station(2))] ;
P3 = [ Long_Intermagnet(Ind_station(3)) latitud(Ind_station(3))] ;
plot([P1(1),P2(1),P3(1),P1(1)],[P1(2),P2(2),P3(2),P1(2)],'B-',...
    'LineWidth',1.5)

% plot distancias estaciones: (lineas rectas)
[arclen,~]=distance(P1,P2);
dist_std(1) = deg2km(arclen);
[arclen,~]=distance(P1,P3);
dist_std(2) = deg2km(arclen);
[arclen,~]=distance(P2,P3);
dist_std(3) = deg2km(arclen);
dist_std=round(dist_std);

text((P1(1)+P2(1))/2,(P1(2)+P2(2))/2,{[num2str(dist_std(1)) ' km']},...
    'Color','b','FontSize',7.5,'HorizontalAlignment'...
    ,'Center','Backgroundcolor',[0.9412    0.9412    0.9412])
text((P1(1)+P3(1))/2,(P1(2)+P3(2))/2,{[num2str(dist_std(2)) ' km']},...
    'Color','b','FontSize',7.5,'HorizontalAlignment'...
    ,'Center','Backgroundcolor',[0.9412    0.9412    0.9412])
text((P2(1)+P3(1))/2,(P2(2)+P3(2))/2,{[num2str(dist_std(3)) ' km']},...
    'Color','b','FontSize',7.5,'HorizontalAlignment'...
    ,'Center','Backgroundcolor',[0.9412    0.9412    0.9412])
%   legend:
legend_text={['Station: ' char(Label_Station(1)) ' Radius: ' num2str(r1)... 
   ' km';'Station: ' char(Label_Station(2)) ' Radius: ' num2str(r2) ...
  ' km';'Station: ' char(Label_Station(3)) ' Radius: ' num2str(r3) ' km']};
text(Max_x+6.5,Max_y,legend_text,'Color','b',...
    'FontSize',8,'HorizontalAlignment'...
    ,'Center','Backgroundcolor',[0.9412    0.9412    0.9412])
% titulo:
title(['Geographical position, Stations: ' char(Label_Station(1)) ', ' ...
    char(Label_Station(2)) ', ' char(Label_Station(3)) ', Date: ' ...
    num2str(Time1) ' to '  num2str(Time2) '.'])
% fin del caso 4. SALIDAS: T1 T2 T3;
%% Filtro de T1 T2 T3
% Se filtran las salidas T1, T2, T3 para generar T sin sismos repetidos
% Se busca si existen eventos repetidos en las listas T1,T2,T3
% Una vez detectado que existen eventos repetidos, necesita eliminarlo 
% de la lista que tenga mayor distancia.
% Tabla T1 y T2
% revisa repetidos en tabla t1 y t2
intersecta = intersect(T1.IRISIDS,T2.IRISIDS);
Mapa_indice = [];
for i=1:length(intersecta)
    ind2=find(intersecta(i)==T1.IRISIDS);
    Mapa_indice=[Mapa_indice; ind2];
end
% display (length(intersecta))
T1(Mapa_indice,:) = []; %vacia el repetido de la tabla 1
% Tabla T3 y T1
% revisa repetidos en tabla t1 y t2
Mapa_indice=[];
for i=1:length(intersecta)
    ind2=find(intersecta(i)==T3.IRISIDS);
    Mapa_indice=[Mapa_indice; ind2];
end
T3(Mapa_indice,:) = []; %vacia el repetido de la tabla 1
% Tabla T3 y T2
% revisa repetidos en tabla t1 y t2
intersecta=intersect(T3.IRISIDS,T2.IRISIDS); 
Mapa_indice=[];
for i=1:length(intersecta)
    ind2=find(intersecta(i)==T3.IRISIDS);
    Mapa_indice=[Mapa_indice; ind2];
end
display (length(intersecta))
T3(Mapa_indice,:) = []; %vacia el repetido de la tabla 1

%% Se crea la tabla total T
T = [];
T = [T1; T2; T3];
T.DistanciasAEstacionKm=[T.DistanciasAEstacionKm(:)...
    T.DistanciasAEstacionKm(:) T.DistanciasAEstacionKm(:)];

% Tabla T tiene variable DistanciasAEstacionKm que necesita ser modificada
% para proporcionar distancia a las 3 estaciones

%% DISTANCIAS ESTACIONES-TERREMOTOS:
grados_km =111.1;
for std = 1:length(Label_Station)
 y_station= ELON(std).*grados_km;
 x_station=ELAT(std).*grados_km;
 for indtablat=1:length(T.Latitud)
  dx_Terremotoalternativo = (x_station -(T.Latitud(indtablat)*grados_km));
  dy_Terremotoalternativo = (y_station -(T.Longitud(indtablat)*grados_km));
  d_st_ev_ = round(sqrt(dx_Terremotoalternativo.^2 + ...
      dy_Terremotoalternativo.^2));
  T.DistanciasAEstacionKm(indtablat,std) = d_st_ev_;
  end
end
% SALIDA T
end
% FIN
% =====================================