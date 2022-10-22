% 24/04/2020
% RELATIVISTIC SEISMIC PRECURSORS (MODIFICACION GENERAL)
% 16/06/220
% Relativistic_Flux_as_Seismic_Precursors_Tensor.m
% --------------------------------------------------------
% GENERA EL FLUJO RELATIVISTA TOTAL
% ===============================================================

function  [W,Rp,Sxy,Sxz,Syz,Tz] = Relativistic_Flux_as_Seismic_Precursors_Tensor(BX,BY,BZ,Gamma,x_st,y_st,z_st,x_ev,y_ev,z_ev,Mag,band);                                                                                  

% ENTRADAS
% ================================================================================
     c = 1;          % VELOCIDAD DE LA LUZ NORMALIZADA
    mu = 5e-2;       % COEFICIENTE DE DECAIMIENTO DE LA MAGNITUD
    so = 1e3;        % REGION DE INFLUENCIA SISMICA EN km
  alfa = 1;          % CONSTANTE DE UNIDADES
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

 ho = 1;   % FILTRAJE

 kd = (1/(4*pi));
   
for st=1:size(BX,3)
    
  B1 = BX(:,:,st);  
  B2 = BY(:,:,st);
  B3 = BZ(:,:,st);
  
  [mb,nb] = size(B1);
  if ho == 1
    B1=B1(:);
    B2=B2(:);
    B3=B3(:);
  else
    B1 = smooth(B1(:),ho);
    B2 = smooth(B2(:),ho);
    B3 = smooth(B3(:),ho);
  end
  B1 = reshape(B1,mb,nb,1);
  B2 = reshape(B2,mb,nb,1);
  B3 = reshape(B3,mb,nb,1);
% 
% MODELO COMPLETO DE DECAIMIENTO DE LA AMPLITUD DEL PRECURSOR
% -----------------------------------------------------------
        % A = c*alfa*beta*(1 + Gamma*mu*Mag - d_st_ev(st)^2/(2*so^2)); 
          A = 1; % CALIBRACION DEL EXPERIMENTO
% ----------------------------------------------------------- 
        B = Gamma*(B1.^2 + B2.^2 + B3.^2*Gamma^(-2));
% ========================================================================
% TENSOR ELECTROMAGNETCO
% NORMALIZACION POR ESTACION DE CADA COMPONENTE DEL TENSOR 
% ----------------------------------------------------------
% COMPONENTE (i=0,j=0)
% -----------------------------
                w = (A)*(kd*Gamma^2*(beta^2 + 1)*(B1.^2 + B2.^2) + B3.^2);              % ENERGIA TOTAL (DEPENDE DE LA POSICION ESTACION-EVENTO (MODULADA POR A))
%                 m = mean(w(:));
%                 m = repmat(m,size(w));
%                 s = std(w(:)); 
%                 s = repmat(s,size(w));
%                 w = (w - m)./s;
%        gw(:,:,st) =  abs(log2(w));
         gw(:,:,st) =  gradient(w);
% ========================================================================
% TENSOR DE STRESS DE MAXWELL
% ----------------------------------------------------------------------------
% COMPONENTES DIAGONALES (i=j)
% PRESION DE RADIACION
% ----------------------------
      
Sigma_xx = kd*(1/2*B.^2 - (Gamma*beta)^2*B2.^2 - B1.^2);            % COMPONENTE X
Sigma_yy = kd*(1/2*B.^2 - (Gamma*beta)^2*B1.^2 - B2.^2);            % COMPONENTE Y
Sigma_zz = kd*(1/2)*(1/2*B.^2 - B3.^2);                             % COMPONENTE Z
      rp = (A)*(1/3)*(Sigma_xx + Sigma_yy + Sigma_zz);              % STRESS VOLUMETRICO (RADIATION PRESSURE)(MODULADA POR (A))
%                 m = mean(rp(:));
%                 m = repmat(m,size(rp));
%                 s = std(rp(:)); 
%                 s = repmat(s,size(rp));
%                rp = (rp - m)./s;
%        gp(:,:,st) = abs(log2(rp));
         gp(:,:,st) = gradient(rp);
 
% COMPONENTES CRUZADAS (i>j)
%
% STRESS DE CIZALLADURA (SHEAR STRESS)(FLUJO DE MOMENTUM)
% ----------------------------------------------------------
         Sigma_xy = 2*(A)*kd*(1-beta^2)*Gamma^2*B1.*B2;  % MODULADA POR (A)
%                 m = mean(Sigma_xy(:));
%                 m = repmat(m,size(Sigma_xy));
%                 s = std(Sigma_xy(:)); 
%                 s = repmat(s,size(Sigma_xy));
%       gxy(:,:,st) = abs(log2(Sigma_xy));
        gxy(:,:,st) = gradient(Sigma_xy);
% ------------------------------------------------------------------- 
         Sigma_xz = 2*(A)*kd*Gamma*B1.*B3;
%                 m = mean(Sigma_xz(:));
%                 m = repmat(m,size(Sigma_xz));
%                 s = std(Sigma_xz(:)); 
%                 s = repmat(s,size(Sigma_xz));
%          Sigma_xz = (Sigma_xz - m)./s;
%       gxz(:,:,st) = abs(log2(Sigma_xz));
        gxz(:,:,st) =  gradient(Sigma_xz);
% ------------------------------------------------------------
           Sigma_yz = 2*(A)*kd*Gamma*B2.*B3;
%                 m = mean(Sigma_yz(:));
%                 m = repmat(m,size(Sigma_yz));
%                 s = std(Sigma_yz(:)); 
%                 s = repmat(s,size(Sigma_yz));
%          Sigma_yz = (Sigma_yz - m)./s;
%       gyz(:,:,st) = abs(log2(Sigma_yz));
        gyz(:,:,st) =  gradient(Sigma_yz);

% TRAZA DEL TENSOR
% ----------------------------------------
                 Tz = w - Sigma_xx - Sigma_yy - Sigma_zz;
%                 m = mean(Tz(:));
%                 m = repmat(m,size(Tz));
%                 s = std(Tz(:)); 
%                 s = repmat(s,size(Tz));
%                Tz = (Tz - m)./s;
%       gtz(:,:,st) = abs(log2(Tz));
        gtz(:,:,st) = gradient(Tz);
        
% VISUALIZACION
% -----------------
% ========================================================================
if band == 1
% ------------------------------------------------------------
        np = 10;
        Nd = 60*24*365.25;
        St = {'St 1','St 2', 'St 3'};
        lx = length(gtz(:))*ones(1,np);
        ly = linspace(min(gtz(:)),max(gtz(:)),np);
       L_x = [L_x;lx];
       L_y = [L_y;ly];
        mb = length(B1(:));
        Mb = max(abs([B1(:);B2(:);B3(:)]));
        tb = (1:mb)/Nd;
        tz = (1:st*mb)/Nd;
        Mf = max(abs([gw(:);gp(:);gxy(:);gxz(:);gyz(:);gtz(:)]))/8;
       ltz = length(gtz(:));
        mg = min(gtz(:));
        Mg = max(gtz(:));
    Lg(st) = length(gtz(:))/Nd;
    MG(st) = mean([mg Mg]);
    ST(st) = St(st);
% ------------------------------------------------------------
set(figure(1),'Position',[5 32 1885 906],'Color','W')
subplot(7,1,1)
    plot(tb,[B1(:),B2(:),B3(:)],'LineWidth',[2]),title(['Bx, By and Bz for station: ' num2str(st) ' and gamma value: ' num2str(Gamma)]),grid
    axis([0 mb/Nd -Mb Mb])
subplot(7,1,2)
    plot(tz,gw(:),'-b','LineWidth',[2]),  title(['W and station:   ' num2str(st)]),grid
    axis([0 ltz/Nd 0 Mf])
 subplot(7,1,3)
    plot(tz,gp(:),'-b','LineWidth',[2]), title(['Rp and station:  ' num2str(st)]),grid
    axis([0 ltz/Nd 0 Mf])
subplot(7,1,4)
    plot(tz,gxy(:),'-b','LineWidth',[2]),title(['Sxy and station: ' num2str(st)]),grid
    axis([0 ltz/Nd 0 Mf])
subplot(7,1,5)
    plot(tz,gxz(:),'-b','LineWidth',[2]),title(['Sxz and station: ' num2str(st)]),grid
    axis([0 ltz/Nd 0 Mf])
subplot(7,1,6)
    plot(tz,gyz(:),'-b','LineWidth',[2]),title(['Syz and station: ' num2str(st)]),grid
    axis([0 ltz/Nd 0 Mf])
subplot(7,1,7)
    plot(tz,gtz(:),'-b','LineWidth',[2]), title(['Tz and station:  ' num2str(st)]),grid    
    axis([0 ltz/Nd 0 Mf])
    hold on,plot(L_x'/Nd,L_y','-r','LineWidth',[4]),hold off
    text(Lg(1:st)-Lg(1)/2,MG(1:st),ST(1:st),'FontWeight','Bold','Color','r','FontSize',[14])
    xlabel('Time in years (x3)')
        go = get(1);
        g1 = go.Children;
    set(g1(1),'Position',[0.0297 0.1100 0.9528 0.0755])
    set(g1(2),'Position',[0.0308 0.2588 0.9517 0.0481])
    set(g1(3),'Position',[0.0313 0.3800 0.9528 0.0481])
    set(g1(4),'Position',[0.0324 0.5013 0.9512 0.0481])
    set(g1(5),'Position',[0.0329 0.6226 0.9512 0.0481])
    set(g1(6),'Position',[0.0334 0.7439 0.9517 0.0481])
    set(g1(7),'Position',[0.0340 0.8651 0.9517 0.0481])
    Ejes_Visibles(1)
    pause
end
        
end

% NORMALIZACION PARA TODAS LAS ESTACIONES
% --------------------------------------------------
% ENERGIA
% --------
          m = mean(gw(:));
          s = std(gw(:));
          m = repmat(m,size(gw));
          s = repmat(s,size(gw));
          W = abs((gw - m)./s); 
% -------------------------------------------------------------------------
% PRESION DE RADIACION
% --------------------
          m = mean(gp(:));
          s = std(gp(:));
          m = repmat(m,size(gp));
          s = repmat(s,size(gp));
         Rp = abs((gp - m)./s);
% -------------------------------------------------------------------------
% FLUJO DE MOMENTUM
% ----------------------
          m = mean(gxy(:));
          s = std(gxy(:));
          m = repmat(m,size(gxy));
          s = repmat(s,size(gxy));
        Sxy = abs((gxy - m)./s); 
% -------------------------------------------------------------------------
          m = mean(gxz(:));
          s = std(gxz(:));
          m = repmat(m,size(gxz));
          s = repmat(s,size(gxz));
        Sxz = abs((gxz - m)./s); 
% -------------------------------------------------------------------------
          m = mean(gyz(:));
          s = std(gyz(:));
          m = repmat(m,size(gyz));
          s = repmat(s,size(gyz));
        Syz = abs((gyz - m)./s);
% -------------------------------------------------------------------------
% TRAZA
% ----------------------
          m = mean(gtz(:));
          s = std(gtz(:));
          m = repmat(m,size(gtz));
          s = repmat(s,size(gtz));
         Tz = abs((gtz - m)./s);
% -------------------------------------------------------------------------
% ==================================================================================     
