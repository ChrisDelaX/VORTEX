%1: profil rectangulaire simple 
%   -> 1 couche (L=1)                
%2: profil rectangulaire avec couche antireflet dans les creux et sur les bosses
%   -> 3 couches (L=3)
%3: profil rectangulaire avec couche antireflet sur les bosses uniquement
%   -> 2 couches (L=2)
%4: profil rectangulaire avec couche antireflet continue
%   -> 2 couches (L=2)

%5: profil rectangulaire simple avec couche d'arret continue
%   -> 2 couches (L=2)
%6: profil rectangulaire avec couche d'arret continue et antireflet dans les creux et sur les bosses
%   -> 4 couches (L=4)
%7: profil rectangulaire avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> 3 couches (L=3) 
%8: profil rectangulaire avec couche d'arret continue et antireflet continue
%   -> 3 couches (L=3) 


%9: profil rectangulaire simple avec couche d'arret discontinue
%   -> 2 couches (L=2)
%10: profil rectangulaire avec couche d'arret discontinue et antireflet dans les creux et sur les bosses
%   -> 4 couches (L=4)
%11: profil rectangulaire avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> 3 couches (L=3)
%12: profil rectangulaire avec couche d'arret discontinue et antireflet continue
%   -> 3 couches (L=3) 


%13: profil trapézoidal simple
%   -> L couches
%14: profil trapézoidal avec couche antireflet dans les creux, sur les bosses et adhérent aux parois obliques
%   -> L+2 couches
%15: profil trapézoidal avec couche antireflet sur les bosses uniquement
%   -> L+1 couches
%16: profil trapézoidal avec couche antireflet continue
%   -> L+1 couches

%17: profil trapézoidal simple avec couche d'arret continue
%   -> L+1 couches
%18: profil trapézoidal avec couche d'arret continue et antireflet dans les creux, sur les bosses et adhérent aux parois obliques
%   -> L+2+1 couches
%19: profil trapézoidal avec couche d'arret continue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%20: profil trapézoidal avec couche d'arret continue et antireflet continue
%   -> L+1+1 couches


%21: profil trapézoidal simple avec couche d'arret discontinue
%   -> L+1 couches
%22: profil trapézoidal avec couche d'arret discontinue et antireflet dans les creux, sur les bosses et adhérent aux parois obliques
%   -> L+2+1 couches
%23: profil trapézoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
%   -> L+1+1 couches
%24: profil trapézoidal avec couche d'arret discontinue et antireflet continue
%   -> L+1+1 couches
%25: profil LETI 1
%26: profil LETI 2 (couche enterrée)
%
%27: profil Fab (réseau blasé)
%   -> L couches

% L=nombre de couche de discrétisation du profil non rectangulaire en sus
% des 2 milieux externes
% Lb=période du réseau
% F_t(l)=facteur de remplissage de chaque couche avec l=1,2,3...L en partant du substrat
% d_t(l)=épaisseur de chaque couche avec l=1,2,3...L en partant du substrat
% E1_t(l)=indice du milieu 1 (le plus bas) de chaque couche
% E2_t(l)=indice du milieu 2 (le plus haut) de chaque couche

switch p
    
    case 1
        %1: profil rectangulaire simple 
        %   -> 1 couche (L=1)         
        L=1;
        d_t(1:L)=d;
        F_t(1:L)=F;
        E1_t(1:L)=E1;       E2_t(1:L)=E2;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
            prof_inv(l,:)=1./Eprof;

        end
        

    case 2
        %2: profil rectangulaire avec couche antireflet dans les creux et sur les bosses
        %   -> 3 couches (L=3)
        L=3;
        d_t(1:L)=d_AR;
        d_t(2)=d-2*d_AR;  %d=épaisseur totale
        F_t(1:L)=F;
        E1_t(1:L)=E1;       E2_t(1:L)=E2;
        E1_t(1)=E_AR;       E2_t(3)=E_AR;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end
        
    case 3
        %3: profil rectangulaire avec couche antireflet sur les bosses uniquement
        %   -> 2 couches (L=2)  
        L=2;
        d_t(1)=d-d_AR;  %d=épaisseur totale
        d_t(2)=d_AR;
        F_t(1:L)=F;
        E1_t(1:L)=E1;       E2_t(1)=E2;
                            E2_t(2)=E_AR;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end

    case 4
        %4: profil rectangulaire avec couche antireflet continue
        %   -> 2 couches (L=2)
        L=2;
        d_t(1)=d-d_AR;  %d=épaisseur totale
        d_t(2)=d_AR;
        F_t(1)=F;
        F_t(2)=1;
        E1_t(1:L)=E1;       E2_t(1)=E2;
                            E2_t(2)=E_AR;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end

    case 5
        %5: profil rectangulaire simple avec couche d'arret continue
        %   -> 2 couches (L=2)
        L=2;
        d_t(1)=d_arret;
        d_t(2)=d-d_arret;
        F_t(1)=1;
        F_t(2)=F;
        E1_t(1:L)=E1;       E2_t(1)=E_arret;
                            E2_t(2)=E2;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end

    case 6
        %6: profil rectangulaire avec couche d'arret continue et antireflet dans les creux et sur les bosses
        %   -> 4 couches (L=4)
        L=4;
        d_t(1:L)=d_AR;
        d_t(1)=d_arret;
        d_t(3)=d-d_arret-2*d_AR;
        F_t(1:L)=F;
        F_t(1)=1;
        E1_t(1:L)=E1;       E2_t(1:L)=E2;
        E1_t(2)=E_AR;       E2_t(1)=E_arret;
                            E2_t(4)=E_AR;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end

    case 7
        %7: profil rectangulaire avec couche d'arret continue et antireflet sur les bosses uniquement
        %   -> 3 couches (L=3) 
        L=3;
        d_t(1)=d_arret;
        d_t(2)=d-d_arret-d_AR;
        d_t(3)=d_AR;
        F_t(1:L)=F;
        F_t(1)=1;
        E1_t(1:L)=E1;       E2_t(1)=E_arret;
                            E2_t(2)=E2;
                            E2_t(3)=E_AR;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end

    case 8
        %8: profil rectangulaire avec couche d'arret continue et antireflet continue
        %   -> 3 couches (L=3) 
        L=3;
        d_t(1)=d_arret;
        d_t(2)=d-d_arret-d_AR;
        d_t(3)=d_AR;
        F_t(1:L)=1;
        F_t(2)=F;
        E1_t(1:L)=E1;       E2_t(1)=E_arret;
                            E2_t(2)=E2;
                            E2_t(3)=E_AR;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end

    case 9
        %9: profil rectangulaire simple avec couche d'arret discontinue
        %   -> 2 couches (L=2)
        Lb=2*Lb_tmp;    %période double ?????  => marche pas !! factorielle...
        L=2;
        d_t(1)=d_arret;
        d_t(2)=d-d_arret;
        F_t(1:L)=F;
        Etmp(1)=E_arret;    E1_t(1)=E_arret;       E2_t(1:L)=E2;
        Etmp(2)=E2;         E1_t(2)=E1;
        for l=1:L
            Eprof=layer9(F_t(l),Etmp(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end

    case 10
        %10: profil rectangulaire avec couche d'arret discontinue et antireflet dans les creux et sur les bosses
        %   -> 4 couches (L=4)
        % PAS REALISTE
        
    case 11
        %11: profil rectangulaire avec couche d'arret discontinue et antireflet sur les bosses uniquement
        %   -> 3 couches (L=3)
        % PAS REALISTE
        
    case 12
        %12: profil rectangulaire avec couche d'arret discontinue et antireflet continue
        %   -> 3 couches (L=3)
        Lb=2*Lb_tmp;    %période double ?????  => marche pas !! factorielle...
        L=3;
        d_t(1)=d_arret;
        d_t(2)=d-d_arret-d_AR;
        d_t(3)=d_AR;
        F_t(1:L)=F;
        F_t(3)=1;
        Etmp(1)=E_arret;    E1_t(1)=E_arret;       E2_t(1:L)=E2;
        Etmp(2)=E2;         E1_t(2)=E1;            E2_t(3)=E_AR;
        Etmp(3)=E_AR;       E1_t(3)=E_AR;
        for l=1:L
            Eprof=layer12(F_t(l),Etmp(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end

    case 13
        %13: profil trapézoidal simple
        %   -> L couches
        clear prof_inv prof
        global prof
        d_t(1:L)=d/L;
        for l=1:L
            y=(l-1/2)*d_t(l); 
            % a l'endroit
            F_t(l)=F+2*y*tan(pente)/Lb;
            % a l'envers
%            F_t(l)=F-2*y*tan(pente)/Lb;   % F obtenu au fond
%            F_t(l)=F-2*(y-d)*tan(pente)/Lb;   % F obtenu au sommet
%            F_t(l)=F-2*(y-d/2)*tan(pente)/Lb;   % F obtenu à d/2
        end
        E1_t(1:L)=E1;       E2_t(1:L)=E2;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
            prof_inv(l,:)=1./Eprof;
        end
          %figure; imagesc(prof)
          %break
        
        %14: profil trapézoidal avec couche antireflet dans les creux, sur les bosses et adhérent aux parois obliques
        %   -> L+2 couches
    case 15
        %15: profil trapézoidal avec couche antireflet sur les bosses uniquement
        %   -> L couches
        pente=rad(1);
        d_t(1:L-1)=(d-d_AR)/(L-1);
        d_t(L)=d_AR;
        for l=1:L-1
            y=(l-1/2)*d_t(l);
            F_t(l)=F-2*y*tan(pente)/Lb;   % F obtenu au fond
%            F_t(l)=F-2*(y-d)*tan(pente)/Lb;   % F obtenu au sommet
%            F_t(l)=F-2*(y-d/2)*tan(pente)/Lb;   % F obtenu à d/2
        end
        F_t(L)=F_t(L);
        E1_t(1:L)=E1;         E2_t(1:L-1)=E2;
                              E2_t(L)=E_AR;       
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
        end
        %figure; imagesc(prof)
        %break
        
        
        %17: profil trapézoidal simple avec couche d'arret continue
        %   -> L+1 couches
        %18: profil trapézoidal avec couche d'arret continue et antireflet dans les creux, sur les bosses et adhérent aux parois obliques
        %   -> L+2+1 couches
        %19: profil trapézoidal avec couche d'arret continue et antireflet sur les bosses uniquement
        %   -> L+1+1 couches
        %20: profil trapézoidal avec couche d'arret continue et antireflet continue
        %   -> L+1+1 couches
        
        
        %21: profil trapézoidal simple avec couche d'arret discontinue
        %   -> L+1 couches
        %22: profil trapézoidal avec couche d'arret discontinue et antireflet dans les creux, sur les bosses et adhérent aux parois obliques
        %   -> L+2+1 couches
        %23: profil trapézoidal avec couche d'arret discontinue et antireflet sur les bosses uniquement
        %   -> L+1+1 couches
        %24: profil trapézoidal avec couche d'arret discontinue et antireflet continue
        %   -> L+1+1 couches
        
        
    case 25 %profil LETI 1
        Lb=2*Lb_tmp;
        pente=1./180*pi;
        F=F_tmp;
        E1=E1_tmp;E2=E2_tmp;
        %L=15;
        L=L_tmp;
        for l=1:L
            if l==L
                d_t(l)=d_AR_tmp;
                for i=1:1024
                    prof_inv(l,i)=1/E_AR_tmp;
                    prof(l,i)=E_AR_tmp;
                end
            else
                d=d_tmp-d_AR_tmp;
                y(l)=(l)*d/(L-1);d_t(l)=d/(L-1);
                x1(l)=round(1024*((d-y(l))*tan(pente)+Lb/4-((1-F)*Lb_tmp)/2)/Lb);
                x2(l)=round(1024*((d-y(l))*tan(pente)+(1-F)*Lb_tmp+Lb/4-((1-F)*Lb_tmp)/2)/Lb);
                x3(l)=round(1024*((y(l)-d)*tan(pente)+Lb_tmp+Lb/4-((1-F)*Lb_tmp)/2)/Lb);
                x4(l)=round(1024*((y(l)-d)*tan(pente)+(2-F)*Lb_tmp+Lb/4-((1-F)*Lb_tmp)/2)/Lb);
                for i=1:1024
                    if i<=x1(l)
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    elseif i<=x2(l)
                        prof_inv(l,i)=1/E1;
                        prof(l,i)=E1;
                    elseif i<=x3(l)
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    elseif i<=x4(l)
                        prof_inv(l,i)=1/E1;
                        prof(l,i)=E1;
                    else
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    end
                end
            end
        end
%         imagesc(prof)
%         break
    
    case 26  
        Lb=2*Lb_tmp;
        pente=1/180*pi;
        F=F_tmp;
        E1=E1_tmp;E2=E2_tmp;
        %L=15;
        L=L_tmp;
        for l=1:L
            if l==1
                d=d_tmp-d_AR_tmp-d_arret_tmp;
                d_t(l)=d_arret_tmp/2;
                for i=1:1024
                    prof_inv(l,i)=1/E_arret_tmp;
                    prof(l,i)=E_arret_tmp;
                end
            elseif l==2
                d=d_tmp-d_AR_tmp-d_arret_tmp;
                d_t(l)=d_arret_tmp/2;y(l)=(l-1/2)*d/(L-3);
                x1(l)=round(1024*((d-y(l))*tan(pente))/Lb);
                x2(l)=round(1024*((d-y(l))*tan(pente)+(1-F)*Lb_tmp)/Lb);
                x3(l)=round(1024*((y(l)-d)*tan(pente)+Lb_tmp)/Lb);
                x4(l)=round(1024*((y(l)-d)*tan(pente)+(2-F)*Lb_tmp)/Lb);
                for i=1:1024
                    if i<=x1(l)
                        prof_inv(l,i)=1/E_arret_tmp;
                        prof(l,i)=E_arret_tmp;
                    elseif i<=x2(l)
                        prof_inv(l,i)=1/E_arret_tmp;
                        prof(l,i)=E_arret_tmp;
                    elseif i<=x3(l)
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    elseif i<=x4(l)
                        prof_inv(l,i)=1/E_arret_tmp;
                        prof(l,i)=E_arret_tmp;
                    else
                        prof_inv(l,i)=1/E_arret_tmp;
                        prof(l,i)=E_arret_tmp;
                    end
                end
                
            elseif l==L
                d_t(l)=d_AR_tmp;
                for i=1:1024,
                    prof_inv(l,i)=1/E_AR_tmp;
                    prof(l,i)=E_AR_tmp;
                end
            else
                
                y(l)=(l-1/2)*d/(L-3);d_t(l)=d/(L-3);
                x1(l)=round(1024*((d-y(l))*tan(pente))/Lb);
                x2(l)=round(1024*((d-y(l))*tan(pente)+(1-F)*Lb_tmp)/Lb);
                x3(l)=round(1024*((y(l)-d)*tan(pente)+Lb_tmp)/Lb);
                x4(l)=round(1024*((y(l)-d)*tan(pente)+(2-F)*Lb_tmp)/Lb);
                for i=1:1024
                    if i<=x1(l)
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    elseif i<=x2(l)
                        prof_inv(l,i)=1/E1;
                        prof(l,i)=E1;
                    elseif i<=x3(l)
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    elseif i<=x4(l)
                        prof_inv(l,i)=1/E1;
                        prof(l,i)=E1;
                    else
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    end
                end
            end
        end

        
    case 27
        Lb=Lb_tmp;F=F_tmp;
        E1=E1_tmp;
        E2=E2_tmp;
        %L=15;
        L=20+L_tmp;
        for l=1:L
            if l==L
                d_t(l)=d_tmp;
                Ftmp1=(1-F)/2*1024;
                Ftmp2=(1+F)/2*1024;
                for i=1:1024,
                    if i<=Ftmp1
                        prof_inv(l,i)=1/E1;   
                        prof(l,i)=E1;
                    elseif i<=Ftmp2
                        prof_inv(l,i)=1/E2;       
                        prof(l,i)=E2;               
                    else
                        prof_inv(l,i)=1/E1;    
                        prof(l,i)=E1;
                    end
                end
            elseif (mod(l,2))==0
                d_t(l)=d_AR2_tmp;
                for i=1:1024,
                    prof_inv(l,i)=1/E_AR2_tmp;   
                    prof(l,i)=E_AR2_tmp;  
                end
            elseif (mod(l,2))==1
                d_t(l)=d_AR1_tmp;
                for i=1:1024,
                    prof_inv(l,i)=1/E_AR1_tmp;   
                    prof(l,i)=E_AR1_tmp;  
                end
            end
        end
        
    case 28 %profil LETI 1
        Lb=2*Lb_tmp;
        pente=1/180*pi;
        F=F_tmp;
        E1=E1_tmp;E2=E2_tmp;
        %L=15;
        L=L_tmp;
        for l=1:L
            if l==L
                d_t(l)=d_AR_tmp;
                for i=1:1024
                    prof_inv(l,i)=1/E_AR_tmp;
                    prof(l,i)=E_AR_tmp;
                end
            else
                d=d_tmp-d_AR_tmp;
                y(l)=(l-1/2)*d/(L-1);d_t(l)=d/(L-1);
                x1(l)=round(1024*((d/2-y(l))*tan(pente)+Lb/4-((1-F)*Lb_tmp)/2)/Lb);
                x2(l)=round(1024*((d/2-y(l))*tan(pente)+(1-F)*Lb_tmp+Lb/4-((1-F)*Lb_tmp)/2)/Lb);
                x3(l)=round(1024*((y(l)-d/2)*tan(pente)+Lb_tmp+Lb/4-((1-F)*Lb_tmp)/2)/Lb);
                x4(l)=round(1024*((y(l)-d/2)*tan(pente)+(2-F)*Lb_tmp+Lb/4-((1-F)*Lb_tmp)/2)/Lb);
                for i=1:1024
                    if i<=x1(l)
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    elseif i<=x2(l)
                        prof_inv(l,i)=1/E1;
                        prof(l,i)=E1;
                    elseif i<=x3(l)
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    elseif i<=x4(l)
                        prof_inv(l,i)=1/E1;
                        prof(l,i)=E1;
                    else
                        prof_inv(l,i)=1/E2;
                        prof(l,i)=E2;
                    end
                end
            end
        end
        
    case 29 %supercellule bruitée
        ncell=64;
        ec=2048/ncell;        
        Lb=Lb_tmp*ncell;
        F=F_tmp;
        d=d_tmp;
        E1=E1_tmp;
        E2=E2_tmp;
        L=L_tmp;
        d_tmp=d/(L);
        for l=1:L
            d_t(l)=d_tmp;
            %         ec=1024/32;
            Ftmp=round(F*(ec));
            FFF=Ftmp/ec
            for i=1:ncell
                %                 tic=round(Ftmp+1*randn(1));
                tic=Ftmp;
                for j=1:ec
                    if j<=tic
                        prof_inv(l,(i-1)*ec+j)=1/E2(l);   
                        prof(l,(i-1)*ec+j)=E2(l);
                    else
                        prof_inv(l,(i-1)*ec+j)=1/E1(l);   
                        prof(l,(i-1)*ec+j)=E1(l);
                    end
                end
            end
        end

    case 69 %profil Fab (réseau blasé)
        %   -> L couches
        clear prof_inv prof
        global prof        
        pente=atan(d/Lb);
        d_t(1:L)=d/Lb;
        for l=1:L
            F_t(l)=l/(L+1);
        end
        E1_t(1:L)=E1;       E2_t(1:L)=E2;
        for l=1:L
            Eprof=layerFab(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
            prof_inv(l,:)=1./Eprof;            
        end
%         figure; imagesc(prof)
%         break        

    
    case 30
        %30: profil proto AGPM N-band (Uppsala 2009)
        %   -> L couches
        clear prof_inv prof
        global prof
        
%         clear
%         close all
%         E1 = 1;
%         E2 = 2;
%         L = 100;
%         d = 13;
%         F = 0.5;
%         Lb = 4.679;
        
        decal_F = F-0.2996;
        M = [0 0.5
            0.0001 0.4995
            0.02 0.38+decal_F/2
            0.5 0.26+decal_F/2
            0.65 0.23+decal_F/2
            0.9 0.15+decal_F/2
            0.97 0.09+decal_F/2
            0.9999 0.0005
            1 0];
% 
% 
%         % dessin
%         figure
%         hold on
%         %axis equal
%         t=size(M,1);
%         ytemp = [0:0.01:1];
%         ytemp = ytemp;%*d;
%         Mp(:,1) = M(:,1);%*d;
%         Mp(:,2) = M(:,2);%*Lb;
%         Pp = fit(Mp(:,1),Mp(:,2),'pchipinterp');
%         xtemp = -Pp(ytemp)+0.5;%*Lb;
%         plot(xtemp,ytemp,'linewidth',2)
%         xtemp =  Pp(ytemp)+0.5;%*Lb;
%         plot(xtemp,ytemp,'linewidth',2)
%         Mp(t+1:2*t-1,1)=Mp(1:t-1,1);
%         Mp(t+1:2*t-1,2)=Mp(1:t-1,2)+0.5;%*Lb;
%         Mp(1:t,2)=-Mp(1:t,2)+0.5;%*Lb;
%         Mp(t+1:2*t-1,:)=flipud(Mp(t+1:2*t-1,:));
%         plot(Mp(:,2),Mp(:,1),'ko')
%         % ----
        
        
        d_t(1:L)=d/L;
        P = fit(M(:,1),M(:,2),'pchipinterp');
        for l=1:L
            y=(l-1/2)/L;
            F_t(l) = 2*P(y);
        end
        E1_t(1:L)=E1;       E2_t(1:L)=E2;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
            prof_inv(l,:)=1./Eprof;
        end
        prof=flipdim(prof,1);
        d_t=flipdim(d_t,2);
        prof_inv=flipdim(prof_inv,1);
        
        
case 31
        %31: Pontus (trapeze + triangle on the top)
        %   -> L couches
        clear prof_inv prof
        global prof
        L1=round(L/6);
        L2=L-L1;
        d_t(1:L1)=d/6/L1;
        d_t(L1+1:L)=d*5/6/L2;
        for l=1:L1
            F_t(l)=(l-1/2)/L1*F;
        end
        for l=L1+1:L
            y=(l-L1-1/2)*d_t(l); 
            % a l'endroit
            F_t(l)=F+2*y*tan(pente)/Lb;
            % a l'envers
%            F_t(l)=F-2*y*tan(pente)/Lb;   % F obtenu au fond
%            F_t(l)=F-2*(y-d)*tan(pente)/Lb;   % F obtenu au sommet
%            F_t(l)=F-2*(y-d/2)*tan(pente)/Lb;   % F obtenu à d/2
        end
        E1_t(1:L)=E1;       E2_t(1:L)=E2;
        for l=1:L
            Eprof=layer(F_t(l),E1_t(l),E2_t(l));
            prof(l,:)=Eprof;
            prof_inv(l,:)=1./Eprof;
        end
%           figure; imagesc(prof)
%           break
end



%Renversement du profile (nécessaire pour l'algo Main)
%-----------------------------------------------------

d_t=flipdim(d_t,2);
prof=flipdim(prof,1);

%not needed: 
%prof_inv=flipdim(prof_inv,1);





