% PORCENTAJE INDICES DAÑADOS
% ====================
% DESCRIPCIÓN: Calcula el porcentaje de índices espurios en un tren de 
% datos otorgado por intermagnet
% SALIDA: Porcentaje de datos dañados.
% =============================

function [Porcentaje]=Porcentaje_Indices_malos(Bxo,Byo,Bzo)
    Pr=[];
    ind_x = find(Bxo(:) >= 90000);
    ind_y = find(Byo(:) >= 90000);
    ind_z = find(Bzo(:) >= 90000);
    Pr(1)=size(ind_x,1)/(size(Bxo,1)*size(Bxo,2));
    Pr(2)=size(ind_y,1)/(size(Byo,1)*size(Byo,2));
    Pr(3)=size(ind_z,1)/(size(Bzo,1)*size(Bzo,2));
    Porcentaje= Pr.*100;
end
% FIN
% =============================