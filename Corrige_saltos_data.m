function [BxN]=Corrige_saltos_data(Bio)
mx=mean(Bio);
mr = repmat(mx,size(Bio,1),1);
BxN = Bio-mr;
%save 'EYR_2003_2006' Bx Bxo By Byo Bz Bzo F Fo