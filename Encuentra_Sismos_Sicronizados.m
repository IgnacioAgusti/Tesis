% Encuentra_Sismos_Sicronizados.m =========================================
% SISMOS LOCALIZADOS DIAS ANTES Y DIAS DESPUES DEL MAXIMO DEL PRECURSOR
% VERIFICA EFICIENCIA DE LOS PRECURSORES RELATIVISTAS
% ===================================================
% DESCRIPCIÓN: Programa determina aciertos y fallos en el precursor para un
% cierto tiempo. 
% SALIDA:
% p_ind_c INDICES ACERTADOS
% p_ind_f INDICES FALSOS POSITIVOS
%==========================================================================

% =========================================================================
%  uo =  2.5;  % UMBRAL CONSIDERADO PARA EL PRECURSOR
% ================================================== m_dh = mean(dh);

% ind_s = find(Lxx(:,1)>=f1 & Lxx(:,1) <= f2);         % INTERVALO DE
% FECHAS (f1 y f2)
%    nc = length(ind_s);                               % NUMERO TOTAL DE
%    SISMOS EN ESTUDIO
% Lxx_u = Lxx(ind_s,:); Lyy_u = Lyy(ind_s,:);
%     nd_b = 0;    % DIAS ANTES DEL EVENTO   (EL SISMO OCURRE (nb_d) ANTES
%     DEL MAXIMO DEL PRECURSOR) nd_a = 90;   % DIAS DESPUES DEL EVENTO (EL
%     SISMO OCURRE (nb_a) DESPUES DEL MAXIMO DEL PRECURSOR)
% -------------------------------------------------------------------------
Amp_U = [];

if band == 1
 % close all
 set(figure(5),'Position',[5 447 1912 531],'Color','W')
end

% =========================================================================
% ELIMINA INDICES FUERA DEL RANGO
% -------------------------------------------------------------------------
As = smooth(As,1);
ind = find(As>=uo);
% Ind_in = find(ind - nd_b >= 0 & ind + nd_a <=  Nv*Md);
Ind_in = find(ind - nd_b > 0 & ind + nd_a <  Nv);
ind = ind(Ind_in);
% --------------------------------------------------------------------
Lxx_c = [];  Lyy_c = [];% SISMOS ENCONTRADOS EN EL RANGO DIAS ANTES Y DIAS DESPUES
Lxx_f = [];  Lyy_f = [];% SISMOS NO-ENCONTRADOS EN EL RANGO DIAS ANTES Y DIAS DESPUES

Ind_c = [];ind_c = [];
Ind_f = [];ind_f = [];

lx = t - f1;
ly = ones(1,length(t))*uo;
kc = 0;
kf = 0;

IND_kc = [];
IND_kf = [];

ind_kc = 0;
for  ind_kc=1:nc
 Io = find(Lxx_u(ind_kc,1) >=  t(round(ind - nd_b)) & Lxx_u(ind_kc,1) <= t(round(ind + nd_a))); % CRITERIO DE BUSQUEDA
 
 if length(Io)>0
  kc = kc + 1;
  Lxx_c = [Lxx_c;Lxx_u(ind_kc,:)]; % SISMOS CIERTOS EN EL RANGO BUSCADO
  Lyy_c = [Lyy_c;Lyy_u(ind_kc,:)];
  Ind_c = [Ind_c;ind(Io)];         % INDICES CORRESPONDIENTES EN EL PRECURSOR (As)
  ind_c = [ind_c;Io'];             % INDICES CORRESPONDIENTES ??
  IND_kc = [IND_kc;ind_kc];
  % ====================================
  %     plot(t-f1,As,'-b','LineWidth',[1]), hold on
  %     plot(t(Ind_c)-f1,As(Ind_c),'+g','LineWidth',[1]),
  %     plot(t(Ind_c)-f1,As(Ind_c),'.g','LineWidth',[1])
  %     plot(t(IND_FP)-f1,As(IND_FP),'.r',t(IND_FP)-f1,As(IND_FP),'+r','LineWidth',[1])
  %
  %     plot(Lxx_u'-f1,Lyy_u','-r')
  %     plot(Lxx_c'-f1,Lyy_c','--g','LineWidth',[2])  % SISMOS ENCONTRADOS
  %     plot(Lxx_f'-f1,Lyy_f','--c','LineWidth',[2])  % SISMOS
  %     NO-ENCONTRADOS plot(lx,ly,'-m'),grid hold off
  % =======================================================================
  % AMPLITUD MEDIA,MAXIMA,DESVIACION,MAGNITUD,DISTANCIA A
  % EVENTO-ESTACON,MEDIA DE INTERVALO HORARIO, DIFERENCIA
  % FECHA(REAL-OBTENIDA) Y ANCHO DE LA VENTANA EN DIAS
  Dr_C = ones(1,ind_kc);
  
  Am = [mean(As(ind(Io)));max(As(ind(Io)));std(As(ind(Io)));Lyy_u(ind_kc,end);Dr_C(ind_kc);m_dh;mean(t(Ind_c)) - Lxx_u(ind_kc,1);nd_a]';
  Amp_U = [Amp_U;Am];
  % =======================================================================
 else
  kf = kf + 1;
  Lxx_f = [Lxx_f;Lxx_u(ind_kc,:)]; % SISMOS CIERTOS EN EL RANGO BUSCADO
  Lyy_f = [Lyy_f;Lyy_u(ind_kc,:)];
  Ind_f = [Ind_f;ind(Io)];
  ind_f = [ind_f;Io'];
  IND_kf = [IND_kf;ind_kc];
 end
end

% STATUS ====== SISMOS ENCONTRADOS Y FALLOS
% --------------------------------
pc = round((kc/nc)*100);                  % SISMOS ACERTADOS
pf = round((kf/nc)*100);                  % SISMO FALLIDOS
% --------------------------------

% =============================================== ELIMINA INDICES REPETIDOS
% ===============================================
IND_C = [];
while length(ind_c)>0
 ind_e = find(ind_c==ind_c(1));
 if length(ind_e)>1
  IND_C = [IND_C;ind_c(1)];
  ind_c(ind_e) = [];
 else
  IND_C = [IND_C;ind_c(1)];
  ind_c(1) = [];
 end
end
% -------------------------------------------------------------------------
IND_F = [];
while length(ind_f)>0
 ind_e = find(ind_f==ind_f(1));
 if length(ind_e)>1
  IND_F = [IND_F;ind_f(1)];
  ind_f(ind_e) = [];
 else
  IND_F = [IND_F;ind_f(1)];
  ind_f(1) = [];
 end
end
% -------------------------------------------------------------------------
% INDICES FALSOS POSITIVOS ---------------------------
ind_fp = [];
for i=1:length(IND_C)
 ind_fp = [ind_fp;find(IND_C(i)==ind)];
end
ind(IND_C) = [];
IND_FP = ind;
% -------------------------------------------------------------------------
ind_im = ind_im + 1;
IM(ind,ind_im) = 1;%As;
% IM(:,ind_im) = As;
% -------------------------------------------------------------------------
if band == 1
 plot(t-f1,As,'-b','LineWidth',[1]),
 hold on
 plot(t(Ind_c)-f1,As(Ind_c),'+g','LineWidth',[1]),
 plot(t(Ind_c)-f1,As(Ind_c),'.g','LineWidth',[1])
 %   plot(t(Ind_f)-f1,As(Ind_f),'+r','LineWidth',[1])
 plot(t(IND_FP)-f1,As(IND_FP),'.r',t(IND_FP)-f1,As(IND_FP),'+r','LineWidth',[1])
 plot(Lxx_u'-f1,Lyy_u','-r')
 plot(Lxx_c'-f1,Lyy_c','--g','LineWidth',[2])  % SISMOS ENCONTRADOS
 plot(Lxx_f'-f1,Lyy_f','--c','LineWidth',[2])  % SISMOS NO-ENCONTRADOS
 plot(lx,ly,'-m'),grid
 xlabel('Time in days')
 ylabel('Magnitude')
 title(['Time interval:  [' num2str(dh) '] mdh: [ ' num2str(m_dh) ' ] Threshold uo: [' num2str(uo) '] ndb: [' num2str(nd_b) '] nda: [' num2str(nd_a) '] (days)' ])
 %ylim([0 40])
 xlim([0 Nv*Md])
 g = get(5);
 g1 = g.Children;
 set(g1(1),'Position',[0.0329    0.0914    0.9582    0.8451])
 hold off
 Ejes_Visibles(5)
 %     set(figure(5),'Position',[419 32 1498  931],'Color','W')
 %     imshow(IM',[0 25]),colormap jet,colorbar axis on,grid
 pause
end
% ============================================ AMPLITUDES DEL PRECURSOR
% LOCALIZADAS EN LA ACTIVACION
% ------------------------------------------------------ Amp_U = []; if
% length(Ind_c)>1
%   Amp_U = sum(As(Ind_c));
%     t_U = mean(t(Ind_c)); %repmat('==',length(Ind_c),1)
%   %disp([num2str(Amp_U) ' == ' datestr(t_U) '==' datestr(Lxx(Ind_c,1))])
% end INDICES ACERTADOS Y FALLOS ==========================================
p_ind_c = round(length(IND_C)/length(Ind_in)*100);           % INDICES ACERTADOS
p_ind_f = round( length(IND_FP)/length(Ind_in)*100 );        % INDICES FALSOS POSITIVOS
% ----------------------------------------------

% disp([pc pf pc+pf p_ind_c p_ind_f p_ind_c + p_ind_f])
% =========================================================================

% % POSICON DE LOS EVENTOS (1ERO,2DO,...,) %
% ------------------------------------------ [md,nd] = size(t_enc); n1 = 1;
% n2 = 1; if n2-n1 > 1
%     subplot(1,2,1)
%     plot(mean(t_dado(1:nd,n1:n2)'),t_enc(1:nd),'*-'),xlim([0 4000]),grid
%     subplot(1,2,2) hist(t_enc(1:nd)'-
%     mean(t_dado(1:nd,n1:n2)')',30),grid,pause(1) xlim([-nd_a nd_a])
% else
%     subplot(1,2,1) plot((t_dado(1:nd,n1:n2)'),t_enc(1:nd),'*-'),xlim([0
%     4000]),grid subplot(1,2,2) hist(t_enc(1:nd)'-
%     (t_dado(1:nd,n1:n2)')',30),grid,pause(1) xlim([-nd_a nd_a]) [xc,yc,p]
%     = polinomio2((t_dado(1:nd,n1:n2)')',t_enc(1:nd)',1,100,0);
%           Pol = [Pol;p];
% end %
% =========================================================================