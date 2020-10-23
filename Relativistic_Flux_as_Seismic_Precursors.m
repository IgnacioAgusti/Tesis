% Relativistic_Flux_as_Seismic_Precursors.m
%
% 16/06/2020
% RELATIVISTIC SEISMIC PRECURSORS (MODIFICACION GENERAL)
% --------------------------------------------------------
% GENERA EL FLUJO RELATIVISTA TOTAL
% ===============================================================

function  Ft = Relativistic_Flux_as_Seismic_Precursors(BX,BY,BZ,Gamma,x_st,y_st,z_st,x_ev,y_ev,z_ev,Mag,band);  

% disp(Gamma)
% ENTRADAS
% ================================================================================
     c = 1;          % VELOCIDAD DE LA LUZ NORMALIZADA
    kx = 1;          % CONTRIBUCION DEL FLUJO EN X
    ky = 1;          % CONTRIBUCION DEL FLUJO EN Y
    kz = 1;          % CONTRIBUCION DEL FLUJO EN Z
    mu = 5e-2;       % COEFICIENTE DE DECAIMIENTO DE LA MAGNITUD
    so = 1e3;        % REGION DE INFLUENCIA SISMICA EN km
  alfa = 1;          % CONSTANTE DE UNIDADES
     g = [];
    gx = [];
    gy = [];
    gz = [];
    Lg = [];
    MG = [];
    ST = {};
   L_x = [];
   L_y = [];
% ================================================================================
% ENERGIA (FACTOR GAMMA)
% =========================
  beta = sqrt(1-1/Gamma^2);
% OBTENCION DE PRECURSORES
% --------------------------------
% CORRECCION POR: DILATACION DE MAGNITUD Y CONTRACCION DE DISTANCIAS (dz)
% ========================================================================
            dx = (x_st - x_ev);
            dy = (y_st - y_ev);
            dz = (z_st - z_ev);
       d_st_ev = sqrt(dx.^2 + dy.^2 + Gamma^(-2)*dz.^2);
% ========================================================================

% FLUJO POR ESTACION
%tic
%close all

ho = 10;         % FILTRAJE

for st=1:size(BX,3)
    
  B1 = BX(:,:,st);  
  B2 = BY(:,:,st);
  B3 = BZ(:,:,st);

%   subplot(3,1,1)
%   plot(B1')
%    subplot(3,1,2)
%   plot(B2')
%    subplot(3,1,3)
%   plot(B3')
%   pause
  [mb,nb] = size(B1);

  
if ho == 1
    B1 = (B1(:));
    B1 = reshape(B1,mb,nb,1);
    B2 = (B2(:));
    B2 = reshape(B2,mb,nb,1);
    B3 = B3(:);
    B3 = reshape(B3,mb,nb,1);
    
else
      %fastsmooth es una funcion mas rapida y con el mismo resultado   %fastsmooth(Bxn(:),10,1,0);
%   B1 = fastsmooth(B1(:),ho,1,0);
%   B1 = reshape(B1,mb,nb,1);
%   B2 = fastsmooth(B2(:),ho,1,0);
%   B2 = reshape(B2,mb,nb,1);
%   B3 = fastsmooth(B3(:),ho,1,0);
%   B3 = reshape(B3,mb,nb,1);
    B1 = smooth(B1(:),ho);
    B1 = reshape(B1,mb,nb,1);
    B2 = smooth(B2(:),ho);
    B2 = reshape(B2,mb,nb,1);
    B3 = smooth(B3(:),ho);
    B3 = reshape(B3,mb,nb,1);
end
% MODELO COMPLETO
% -----------------------------------------------------------
   %A = c*alfa*beta*(1 + Gamma*mu*Mag - d_st_ev(st)^2/(2*so^2));
   %  disp(A),pause
  A = 1; % NO HACE USO DEL MODELO (CALIBRACION DEL EXPERIMENTO)
% --------------------------------------------------------------
  B = Gamma*(B1.^2 + B2.^2 + B3.^2*Gamma^(-2));
  
    % flujo fx
    % ---------    
          Bxz = B1.*B3; 
           fx = A.*(Bxz./B)*kx;
   gx(:,:,st) = abs(gradient(fx));
    % flujo fy
    % ---------
          Byz = B2.*B3;
           fy = A.*(Byz./B)*ky;
   gy(:,:,st) = abs(gradient(fy));
    % flujo fz
    % ---------
          Bxy = B1.^2 + B2.^2;
           Bz = B/Gamma;
           fz = -A.*(Bxy./Bz)*kz;
   gz(:,:,st) = abs(gradient(fz));
       %fz = fz/max(abs(fz(:)));
    % flujo total
    % --------------------------------------
              %ft = sqrt(fx.^2 + fy.^2 + fz.^2);
              ft = fx + fy + fz; 
%                m = mean(ft(:));
%                m = repmat(m,size(ft));
%                s = std(ft(:)); 
%                s = repmat(s,size(ft));
               %ft = (ft - m)./s;     
        g(:,:,st) = abs(gradient(ft));
      %   g(:,:,st) = ft; 
% ========================================================================
if band == 1
% VISUALIZACION
% -----------------
        np = 10;
        Nd = 60*24*365.25;
        St = {'St 1','St 2', 'St 3'};
        lx = length(gx(:))*ones(1,np);
        ly = linspace(min(g(:)),max(g(:)),np);
       L_x = [L_x;lx];
       L_y = [L_y;ly];
        mb = length(B1(:));
        Mb = max(abs([B1(:);B2(:);B3(:)]));
        tb = (1:mb)/Nd;
        tz = (1:st*mb)/Nd;
        Mf = max(abs([gx(:);gy(:);gz(:)]));
       ltz = length(g(:));
        mg = min(g(:));
        Mg = max(g(:));
    Lg(st) = length(g(:))/Nd;
    MG(st) = mean([mg Mg]);
    ST(st) = St(st);
% ------------------------------------------------------------
set(figure(1),'Position',[5 32 1885 906],'Color','W')
subplot(5,1,1)
    plot(tb,[B1(:),B2(:),B3(:)],'LineWidth',[2]),title(['Bx, By and Bz for station: ' num2str(st) ' and gamma value: ' num2str(Gamma)]),grid
    axis([0 mb/Nd -Mb Mb])
subplot(5,1,2)
    plot(tz,gx(:),'-b','LineWidth',[2]),title(['fx and station: ' num2str(st)]),grid
    axis([0 tz(end) 0 Mf])
 subplot(5,1,3)
    plot(tz,gy(:),'-b','LineWidth',[2]),title(['fy and station: ' num2str(st)]),grid
    axis([0 tz(end) 0 Mf])
subplot(5,1,4)
    plot(tz,gz(:),'-b','LineWidth',[2]),title(['fz and station: ' num2str(st)]),grid
    axis([0 tz(end) 0 Mf])
subplot(5,1,5)
    plot(tz,g(:),'-b','LineWidth',[2]),title(['Ft and stations: ' num2str(1) ' to ' num2str(st) ])
    axis([0 tz(end) 0 Mg])
    hold on,plot(L_x'/Nd,L_y','-r','LineWidth',[4]),hold off
    xlabel('Time in years (x3)')
    text(Lg(1:st)-Lg(1)/2,MG(1:st),ST(1:st),'FontWeight','Bold','Color','r','FontSize',[14])

        go = get(1);
        g1 = go.Children;
    set(g1(1),'Position',[0.0244 0.1100 0.9560 0.1124])
    set(g1(2),'Position',[0.0265 0.2827 0.9544 0.1124])
    set(g1(3),'Position',[0.0276 0.4553 0.9554 0.1124])
    set(g1(4),'Position', [0.0286 0.6280 0.9560 0.1124])
    set(g1(5),'Position',[0.0297 0.8007 0.9554 0.1124])
    grid
    Ejes_Visibles(1)
    pause
end
end
%toc
% NORMALIZACION DEL FLUJO PARA TODAS LAS ESTACIONES
% --------------------------------------------------
     m = mean(g(:));
     m = repmat(m,size(g));
     s = std(g(:)); 
     s = repmat(s,size(g));
    Ft = abs((g - m)./s);
% ==================================================================================     
if band == 1
    fp1 = Ft(:,:,1);
    fp2 = Ft(:,:,2);
    fp3 = Ft(:,:,3);
     Mf = max([fp1(:);fp2(:);fp3(:)]);
     tz = (1:length(fp1(:)))/Nd;
     
    set(figure(2),'Position',[5 32 1885 906],'Color','W')
    subplot(3,1,1)
    
    plot(tz,fp1(:),'LineWidth',[2]),title(['Fp for station ' num2str(1) ' and gamma value: ' num2str(Gamma)]),grid
    axis([0 tz(end) 0 Mf])
    subplot(3,1,2)
     
    plot(tz,fp2(:),'LineWidth',[2]),title(['Fp for station ' num2str(2) ' and gamma value: ' num2str(Gamma)]),grid
    axis([0 tz(end) 0 Mf])
    subplot(3,1,3)
     
    plot(tz,fp3(:),'LineWidth',[2]),title(['Fp for station ' num2str(3) ' and gamma value: ' num2str(Gamma)]),grid
    axis([0 tz(end) 0 Mf])
    xlabel('Time in years (x3)')
        go = get(2);
        g1 = go.Children;
    set(g1(1),'Position',[0.0276 0.1100 0.9507 0.2157])
    set(g1(2),'Position',[0.0286 0.4111 0.9480 0.2143])
    set(g1(3),'Position',[0.0276 0.7107 0.9491 0.2143])
    Ejes_Visibles(2)
end
% ==================================================================================     


