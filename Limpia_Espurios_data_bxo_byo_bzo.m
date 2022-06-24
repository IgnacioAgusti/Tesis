% LIMPA_ESPURIOS_DATA_BXO_BYO_BZO
% ====================
% DESCRIPCIÓN: Dado un tren de datos Bxo, Byo, Bzo en el formato:
% size(B(x,y,z)o) = Minutos x día
% Encontrará los valores espurios marcado por intermagnet, una vez
% encontrado, determinará si dicho espurio tiene vecinos espurios, si no
% tiene, el punto aislado espurio es sustituido por la media de los puntos
% adyacentes, si tiene vecinos espurios, es sustuido por datos buenos del
% día pasado.
% SALIDA: B(x,y,z)o En el mismo formato de entrada y sin espurios en el
% tren de datos.
%=============================
function [Bxo,Byo,Bzo]=Limpia_Espurios_data_bxo_byo_bzo(Bxo,Byo,Bzo,band)
Ind = [];ind_x=[];ind_y=[];ind_z=[];close all;
%h = waitbar(0,'Espere...limpiando las espigas de Bx,By and Bz');

%para aquellos que no tengan un punto final o al principio,se
%sustituye la datala del día anterior (la faltante)

for j=1:size(Bxo,1)
    ind_x = find(Bxo(j,:) >= 90000);
    ind_y = find(Byo(j,:) >= 90000);
    ind_z = find(Bzo(j,:) >= 90000);
    
    if band==1
        %       waitbar(j/size(Bxo,1),h)
        subplot(2,3,1);hold on
        plot(Bxo(j,:));title(' datos no corregido Bx');xlabel('Minutos');
        ylabel('campo magnetico');
        subplot(2,3,2);hold on
        plot(Byo(j,:));title(' datos no corregido By');xlabel('Minutos');
        ylabel('campo magnetico');
        subplot(2,3,3);hold on;
        plot(Bzo(j,:));title(' datos no corregido Bz');xlabel('Minutos');
        ylabel('campo magnetico');
        pause(0.0001);
    end
    if ((length(ind_x) >0)  | (length(ind_y) >0) |(length(ind_z) >0))
        Ind = [Ind; j];
        %función sustituye un vector de datos con data corrupta
        % por la media de los puntos más cercanos
        %PARA BXO
        %_______________________________________________________
        n=0;
        while ((length(ind_x) >0))
            n=n+1;
            if( (length(ind_x) >0) )
                ig=Bxo(j,:);
                a=ind_x;z=[];
                %Los siguientes condiconales son para cuando no hay punto
                % final o inicial para hacer la media
                if ((ig(1,size(ig,2))== 99999) | ((ig(1,1)== 99999)))
                    if (ig(1,size(ig,2))== 99999)
                        %cuando eso pasa,se rellena con la data del día 
                        %anterior a la misma hora
                        Bxo(j,a(1):size(Bxo,2))=Bxo(j-n,a(1):size(Bxo,2));
                    else
                        Bxo(j,1:length(a))=(Bxo(j-n,1:length(a)));
                    end
                else
                    %función que recibe un vector de datos con indices
                    %corruptos y devuelve un vector del mismo tamaño 
                    %pero limpio (usando la media de los vecinos)
                    
                    for t=1:length(a)
                        % ya que el vector esta ordenado al valor previo 
                        % estamos seguros que es correcto
                        prev=ig(a(t) -1);
                        z=t;
                        % si el proximo valor tambien esta corrupto
                        % vemos el de más arriba hasta encontrar uno limpio
                        while((z+1)<= length(a) && (a(z)+1 ==a(z+1)))
                            z=z+1;
                        end
                        next=ig(a(z)+1);
                        ig(a(t))=(prev+next)/2;
                    end
                    Bxo(j,:)=ig;
                end
            end
            ind_x = find(Bxo(j,:) >= 99998);
            
        end
        
        %PARA BYO
        %______________________________________________________
        n=0;
        while ((length(ind_y) >0))
            n=n+1;
            if( (length(ind_y) >0) )
                ig=Byo(j,:);a=ind_y;z=[];prev=[];next=[];
                if ((ig(1,size(ig,2))== 99999) | ((ig(1,1)== 99999)))
                    if (ig(1,size(ig,2))== 99999)
                        Byo(j,a(1):size(Byo,2))=Byo(j-n,a(1):size(Byo,2));
                    else
                        Byo(j,1:(length(a)))=(Byo(j-n,1:length(a)));
                    end
                else
                    for t=1:length(a)
                        %ya que el vector esta ordenado al valor previo 
                        %estamos seguros que es correcto
                        prev=ig(a(t) -1);
                        z=t;
                        %si el proximo valor tambien esta corrupto vemos 
                        % el de más arriba hasta encontrar uno limpio
                        while((z+1)<= length(a) && (a(z)+1 ==a(z+1)))
                            z=z+1;
                        end
                        next=ig(a(z)+1);
                        ig(a(t))=(prev+next)/2;
                    end
                    Byo(j,:)=ig;ig=[];
                end
            end
            ind_y = find(Byo(j,:) >= 99998);
            
        end
        %PARA BZO
        %______________________________________________________
        n=0;
        while ((length(ind_z) >0))
            n=n+1;
            
            if( (length(ind_z) >0) )
                
                ig=Bzo(j,:);a=ind_z;z=[];prev=[];next=[];
                if ((ig(1,size(ig,2))== 99999) | ((ig(1,1)== 99999)))
                    if (ig(1,size(ig,2))== 99999)
                        Bzo(j,a(1):size(Bzo,2))=Bzo(j-n,a(1):size(Bzo,2));
                    else
                        Bzo(j,1:(length(a)))=(Bzo(j-n,1:length(a)));
                    end
                else
                    for t=1:length(a)
                        %ya que el vector esta ordenado al valor previo 
                        %estamos seguros que es correcto
                        prev=ig(a(t) -1);
                        z=t;
                        %si el proximo valor tambien esta corrupto vemos 
                        %el de más arriba hasta encontrar uno limpio
                        while((z+1)<= length(a) && (a(z)+1 ==a(z+1)))
                            z=z+1;
                        end
                        next=ig(a(z)+1);
                        ig(a(t))=(prev+next)/2;
                    end
                    Bzo(j,:)=ig;
                end
            end
            
            ind_z = find(Bzo(j,:) >= 99998);
            
        end
    end
    %Muestra los datos sin espurios
    if band==1
        hold on;subplot(2,3,4);
        plot(Bxo(j,:));title('datos sin espurios Bx');xlabel('Minutos');
        ylabel('campo magnetico');
        hold on;subplot(2,3,5);
        plot(Byo(j,:));title('datos sin espurios By');xlabel('Minutos');
        ylabel('campo magnetico');
        hold on;subplot(2,3,6);hold on
        plot(Bzo(j,:));title('datos sin espurios Bz');xlabel('Minutos');
        ylabel('campo magnetico');
    end
    
end
end
% FIN
% =========================================================================