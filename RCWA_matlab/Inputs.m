%Configuration g�n�rale du probl�me
%----------------------------------

%Balayage
%--------

%Choix param�tre
%---------------

%REMARQUE: Unit�s en MICRONS !!!

%1: p�riode, pas du r�seau
%2: �paisseur
%3: facteur de remplissage Fa, grande base du trap�ze
%4: facteur de remplissage Fb, petite base du trap�ze
%5: angle d'incidence non conique
%6: angle d'incidence conique
%7: angle de polarisation
%8: longueur d'onde 
%9: �paisseur A/R
%10: NA

scan1=10;
mesh1=10;
    scan2=10;
    mesh2=20;
        scan3=10;
        mesh3=15;
            scan4=5;
            mesh4=20;
                scan5=8;
                mesh5=40;

%Param�tres de calculs
%---------------------

    %Troncature
    
        %En X (N donc 2N+1 ordres au total)  
    
        N=6;
              
%Param�tres du r�seau
%--------------------    

%0.7076    2.7115    0.7188    0.4040
   
                
    
%Pas du r�seau
             
    Lb=0.9561     ;
    Lbmin=0.45;Lbmax=0.5;
    
%Epaisseurs

    d=1.0023  ;
    dmin=0.75;dmax=1.5;
    
%Permittivit�s
%-------------

%Choix mat�riau
    
    %1: Vide/Air
    %2: CdTe
    %3: Diamant
    %4: Germanium
    %5: Silicium
    %6: ZnSe
    %7: YF3
    %8: Manuel
    E_man=1.57;
        %LaF2: E_man=1.57;
        %BaF2: E_man(l=1.97 microns)=1.4647
        %BaF2: E_man(l=2.3253 microns)=1.4635
        %BaF2: E_man(l=2.6738 microns)=1.4623
    
    %Permittivit�s milieu ext�rieur �mergent (!!!! INVERSION)
    
    EI_choice=2;
    
    %Permittivit�s milieu ext�rieur incident (!!!! INVERSION)
    
    EIII_choice=1;
    
    %Permittivit�s du r�seau
    %-----------------------
    
        %R�seau en relief de surface
        %---------------------------
        
            %Permittivit� milieu 1
            
            E1_choice=1;
            
            %Permittivit� milieu 2
            
            E2_choice=2;
                                   
            %Facteurs de remplissage
            
            Fa=0.15 ;% Grande Base trap�ze
            Famin=0.65;Famax=0.8;
            
            %Fb=0.3;
            Fb=Fa;% Petite Base trap�ze
            if Fb==Fa  
            Fbmin=Famin;Fbmax=Fbmin;
            else
            Fbmin=0.5;Fbmax=0.9;
            end;
            
        %Permittivit�s pseudo A/R
        %------------------------
           
            E_AR_choice=1;
            d_AR=0.413;
            d_AR_min=0.320;d_AR_max=0.42;
              
        %R�seau en volume
        %----------------
                
            %Permittivit� moyenne
            
            %epsm(l)=;
            
            %Modulation de permittivit�
            
            %modmin(l)=;
            %modmax(l)=;
            %mod(l)=;                       
            
%Nombre de couches
%-----------------

    L_wo_AR=1;%hors A/R
    
    %Sans Pseudo-A/R YF3 : 1
    
    ar=1;%par d�faut ar=1
       
    %Avec Pseudo-A/R : 2
    
    %Avec Pseudo-A/R sur r�seau trap�zoidal : 3
    
%     if ar==2
%         if Fb==Fa
%         ar=2;%3 couches
%         else
%         %TBD!!!
%         ar=3;%au moins 5 couches, pour la pr�cision de la discr�tisation
%         end;
%     end;
    
%Onde incidente 
%--------------

    %Longueur d'onde
    
    lb=4;
    lbmin=6;lbmax=18;%nlb=32;

    %Angle incidence non conique
    
    
    theta=0;
    %theta=39.5*pi/180;
    thetamin=1.1495-30/60/180*pi;thetamax=1.1495+30/60/180*pi;
   
    %Angle incidence conique

    phi=8/60/180*pi;
    phimin=0;phimax=0;
        
    %Angle polarisation
    
    psi=pi/4;
    psimin=0;psimax=0;
    
        
        
    

    

    
