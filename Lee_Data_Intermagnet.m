% LEE DATA INTERMAGNET
% ====================
% DESCRIPCIÓN: Este programa permite la lectura de una determinada carpeta
% con los documentos descomprimidos sacados de INTERMAGNET, elimina los
% espurios, lleva la señal a media = 0 y lo guarda en con el formato de
% nombre deseado.
% SALIDA: 'Japon_MMB_2004_2021' Bx Bxo By Byo Bz Bzo F Fo Bxo_Mo Byo_Mo ...
%    Bzo_Mo Bxo_Sp Byo_Sp Bzo_Sp Porcentaje_Datos_Malos;
%=============================
band=0;
clear all;close all;clc
band=0;
% DIRECTORIO:
File= 'C:\Tesis\Data\Japon_2001_2015\Japon_MMB\Japon_MMB_2004_2021';
eval (['cd ' File]);
% EMPIEZA A CARGAR LOS DATOS
F = [];
a = dir ('*.min');
h = waitbar(0,'Please wait...preparing files');
for i=1:length(a)
    b = a(i).name;
    c = b(4:11);
    y = str2double(c(1:4));
    m = str2double(c(5:6));
    d = str2double(c(7:8));
    F(i) = datenum([y,m,d]);
    waitbar(i/length(a),h)
end
close(h)
[~,jo] = sort(F);
a = a(jo);
Bxo = zeros(length(F),1434);
Byo =  zeros(length(F),1434);
Bzo =  zeros(length(F),1434);
Fo =  zeros(length(F),1434);

h = waitbar(0,'Please wait...generating format of data Bx,By and Bz');

for i=1:length(a)
    %eval(['cd ' a(i).name])
    c = importdata(a(i).name);
    C = c(27:end,:);
    cm = cell2mat(C);
    fo = datenum(datestr(cm(:,1:24)));
    bx = str2num(cm(:,33:42));
    by = str2num(cm(:,42:50));
    bz = str2num(cm(:,52:60));
    vb = 1:1434;%length(bx); % SI REPORTA ERROR USAR 1437
    Bxo = [Bxo;bx(vb)'];
    Byo = [Byo;by(vb)'];
    Bzo = [Bzo;bz(vb)'];
    Fo = [Fo;fo(vb)'];
    clc
    disp([round((i/length(a))*100)])
    
end
% VISUALIZACION
% ======================
if band == 1
    Fs = Fo(:);
    t = (Fs - Fs(1))/365;               % TIEMPO EN AÑOS
    set(figure(1),'Position',[3 414 1021 249],'Color','W')
    plot(t,Bxo(:),'-r'),hold on
    plot(t,Byo(:),'-g')
    plot(t,Bzo(:),'-b'),grid,hold off
    title('Bx (r),By (g) and Bz (b)')
    xlabel('Time in years')
    ylabel ('Amplitude in nT')
end
cd 'C:\Tesis\Data\Japon_2001_2015'; %salida de los documentos
[Porcentaje_Datos_Malos]=Porcentaje_Indices_malos(Bxo,Byo,Bzo);
[Bxo_Sp,Byo_Sp,Bzo_Sp]=Limpia_Espurios_data_bxo_byo_bzo(Bxo,Byo,Bzo,band);
[Bxo_Mo]=Corrige_saltos_data(Bxo_Sp'); %Corrige el salto en la data
[Byo_Mo]=Corrige_saltos_data(Byo_Sp'); %Corrige el salto en la data
[Bzo_Mo]=Corrige_saltos_data(Bzo_Sp'); %Corrige el salto en la data
Bxo_Sp=Bxo_Sp';  Byo_Sp=Byo_Sp';  Bzo_Sp=Bzo_Sp';
Bx = Bxo_Sp(:);
By = Byo_Sp(:);
Bz = Bzo_Sp(:);
F = Fo(:);
Bxo = Bxo';
Byo = Byo';
Bzo = Bzo';
Fo = Fo';
% SALIDA:
save 'Japon_MMB_2004_2021' Bx Bxo By Byo Bz Bzo F Fo Bxo_Mo Byo_Mo ...
    Bzo_Mo Bxo_Sp Byo_Sp Bzo_Sp Porcentaje_Datos_Malos;
% FIN
% =========================================================================