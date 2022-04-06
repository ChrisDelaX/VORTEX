%Configuration générale du problème
%----------------------------------

%Balayage
%--------

%Choix paramètre
%---------------

%REMARQUE: Unités en MICRONS !!!

%1: période, pas du réseau
%2: épaisseur
%3: facteur de remplissage Fa, grande base du trapèze
%4: facteur de remplissage Fb, petite base du trapèze
%5: angle d'incidence non conique
%6: angle d'incidence conique
%7: angle de polarisation
%8: longueur d'onde 
%9: épaisseur A/R
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

%Paramètres de calculs
%---------------------

    %Troncature
    
        %En X (N donc 2N+1 ordres au total)  
    
        N=6;
              
%Paramètres du réseau
%--------------------    

%0.7076    2.7115    0.7188    0.4040
   
                
    
%Pas du réseau
             
    Lb=0.9561     ;
    Lbmin=0.45;Lbmax=0.5;
    
%Epaisseurs

    d=1.0023  ;
    dmin=0.75;dmax=1.5;
    
%Permittivités
%-------------

%Choix matériau
    
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
    
    %Permittivités milieu extérieur émergent (!!!! INVERSION)
    
    EI_choice=2;
    
    %Permittivités milieu extérieur incident (!!!! INVERSION)
    
    EIII_choice=1;
    
    %Permittivités du réseau
    %-----------------------
    
        %Réseau en relief de surface
        %---------------------------
        
            %Permittivité milieu 1
            
            E1_choice=1;
            
            %Permittivité milieu 2
            
            E2_choice=2;
                                   
            %Facteurs de remplissage
            
            Fa=0.15 ;% Grande Base trapèze
            Famin=0.65;Famax=0.8;
            
            %Fb=0.3;
            Fb=Fa;% Petite Base trapèze
            if Fb==Fa  
            Fbmin=Famin;Fbmax=Fbmin;
            else
            Fbmin=0.5;Fbmax=0.9;
            end;
            
        %Permittivités pseudo A/R
        %------------------------
           
            E_AR_choice=1;
            d_AR=0.413;
            d_AR_min=0.320;d_AR_max=0.42;
              
        %Réseau en volume
        %----------------
                
            %Permittivité moyenne
            
            %epsm(l)=;
            
            %Modulation de permittivité
            
            %modmin(l)=;
            %modmax(l)=;
            %mod(l)=;                       
            
%Nombre de couches
%-----------------

    L_wo_AR=1;%hors A/R
    
    %Sans Pseudo-A/R YF3 : 1
    
    ar=1;%par défaut ar=1
       
    %Avec Pseudo-A/R : 2
    
    %Avec Pseudo-A/R sur réseau trapézoidal : 3
    
%     if ar==2
%         if Fb==Fa
%         ar=2;%3 couches
%         else
%         %TBD!!!
%         ar=3;%au moins 5 couches, pour la précision de la discrétisation
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
    
        
        
    

    

    
