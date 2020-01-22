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
%9: NA

scan1=9;

    mesh1=100;

    scan2=9;
    
        mesh2=100;
        
            scan3=5;
    
            mesh3=100;
            
                scan4=8;
    
                mesh4=20;

%Paramètres de calculs
%---------------------

    %Troncature
    
        %En X (N donc 2N+1 ordres au total)  
    
        N=6;
              
%Paramètres du réseau
%--------------------    
    
%Nombre de couches

    L=5;

%Pas du réseau
    
    Lb=0.9617;
    Lbmin=;Lbmax=;
    
%Epaisseurs

    d=0.7647;
    dmin=;dmax=;
    
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
    
    %Permittivités milieu extérieur incident
    
    EI_choice=2;
    
    %Permittivités milieu extérieur émergent
    
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
            
            Fa=0.1979;% Grande Base trapèze
            Famin=;Famax;
            Fb=0.65;% Petite Base trapèze
            Fbmin=;Fbmax=;
            
        %Réseau en volume
        %----------------
                
            %Permittivité moyenne
            
            %epsm(l)=;
            
            %Modulation de permittivité
            
            %modmin(l)=;
            %modmax(l)=;
            %mod(l)=;
                                        
%Onde incidente 
%--------------

    %Longueur d'onde
    
    lb=0.633;
    lbmin=6;lbmax=11;

    %Angle incidence non conique
    
    theta=0.6985;%0*pi/180;
    thetamin=;thetamax=;
   
    %Angle incidence conique

    phi=0*pi/180;
    phimin=;phimax=;
        
    %Angle polarisation
    
    psi=pi/4;
    psimin=;psimax=;
    
        
        
    

    

    
