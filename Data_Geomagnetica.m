% DATA GEOMAGNETICA DE LAS ESTACIONES A ESTUDIAR
% ==============================================
function [Data,Label_Station] = Data_Geomagnetica(n);
%clc
Data = [];
Label_Station = [];
% clc
% disp('(1) Japon (2) Europa (3) Otras estaciones')
switch n
 case 1
  disp('Data seleccionada correspondiente a Japon')
  % JAPON
  % -----------------------------------
  Data_1 = 'Japon_KAK_2004_2015.mat';
  Data_2 = 'Japon_KNY_2004_2015.mat';
  Data_3 = 'Japon_MMB_2004_2015.mat';
  % -----------------------------------------------------------------
  Label_Station = {'KAK','KNY','MMB'};% ETIQUETAS DE LAS ESTACIONES JAPON
  % -----------------------------------------------------------------
 case 12
  disp('Data seleccionada correspondiente a Japon(solo KNY)')
  % JAPON
  % -----------------------------------
  Data_1 = 'Japon_KNY_2004_2021.mat';
  Data_2 = 'Japon_KNY_2004_2021.mat';
  Data_3 = 'Japon_KNY_2004_2021.mat';
  % -----------------------------------------------------------------
  Label_Station = {'KNY','KNY','KNY'};% ETIQUETAS DE LAS ESTACIONES JAPON
  % -----------------------------------------------------------------
 case 13
  disp('Data seleccionada correspondiente a Japon(solo MMB)')
  % JAPON
  % -----------------------------------
  Data_1 = 'Japon_MMB_2004_2021.mat';
  Data_2 = 'Japon_MMB_2004_2021.mat';
  Data_3 = 'Japon_MMB_2004_2021.mat';
  % -----------------------------------------------------------------
  Label_Station = {'MMB','MMB','MMB'};% ETIQUETAS DE LAS ESTACIONES JAPON
  % -----------------------------------------------------------------
 case 2
  disp('Data seleccionada correspondiente a Europa (DOU-MAB-WNG)')
  Data_1 = 'DOU_2005_2013.mat';
  Data_2 = 'MAB_2005_2013.mat';
  Data_3 = 'WNG_2005_2013.mat';
  % -----------------------------------------------------------------------
  Label_Station = {'DOU','MAB','WNG'};% ETIQUETAS DE LAS ESTACIONES EUROPA
  % -----------------------------------------------------------------------
 case 3
  disp('Data seleccionada correspondiente a USA')
  Data_1 = 'USA_NEW_1998_2007_P3.mat';
  Data_2 = 'USA_SIT_1998_2007_P3.mat';
  Data_3 = 'USA_VIC_1998_2007_P3.mat';
  % -----------------------------------------------------------------------
  Label_Station = {'NEW','SIT','VIC'}; % ETIQUETAS DE LAS ESTACIONES USA
  % -----------------------------------------------------------------------
 case 4
  disp('Data seleccionada correspondiente a CHILE')
  Data_1 = 'Chile_HUA_2002_2011.mat';
  Data_2 = 'Chile_HUA_2002_2011.mat';
  Data_3 = 'Chile_HUA_2002_2011.mat';
  % -----------------------------------------------------------------------
  Label_Station = {'HUA','HUA','HUA'};  % ETIQUETAS DE LAS ESTACIONES CHILE
  % -----------------------------------------------------------------------
 case 5
  
  disp('Data seleccionada correspondiente a SUR DE EUROPA (DUR)')
  Data_1 = 'DUR_2016_2020.mat';
  Data_2 = 'LON_2016_2020.mat';
  Data_3 = 'THY_2016_2020.mat';
  % -----------------------------------------------------------------------
  Label_Station = {'DUR','LON','THY'};  % ETIQUETAS DE LAS ESTACIONES (DUR)
  % -----------------------------------------------------------------------
 case 6
  
  disp('Data seleccionada correspondiente a SUR DE EUROPA (PEG-PAG-IZN)')
  Data_1 = 'PEG_2016_2020.mat';
  Data_2 = 'PAG_2016_2020.mat';
  Data_3 = 'IZN_2016_2020.mat';
  % -----------------------------------------------------------------------
  Label_Station = {'PEG','PAG','IZN'};  % ETIQUETAS DE LAS ESTACIONES (PEG)
  % -----------------------------------------------------------------------
 case 7
  
  disp('Data seleccionada correspondiente muy al norte de Japon (PET)')
  Data_1 = 'PET_2007_2020.mat';
  Data_2 = 'PET_2007_2020.mat';
  Data_3 = 'PET_2007_2020.mat';
  % -----------------------------------------------------------------------
  Label_Station = {'PET','PET','PET'}; % ETIQUETAS DE LAS ESTACIONES EUROPA
  % -----------------------------------------------------------------------
 otherwise
  disp('Elija otra estacion')
  return
end
% -------------------------------------------------------------------------
Data = [Data_1;Data_2;Data_3];    % ETIQUETAS DE LOS DIRECTORIOS
% -------------------------------------------------------------------------

