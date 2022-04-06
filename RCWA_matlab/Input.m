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
%9: NA

scan1=9;

    mesh1=100;

    scan2=9;
    
        mesh2=100;
        
            scan3=5;
    
            mesh3=100;
            
                scan4=8;
    
                mesh4=20;

%Param�tres de calculs
%---------------------

    %Troncature
    
        %En X (N donc 2N+1 ordres au total)  
    
        N=6;
              
%Param�tres du r�seau
%--------------------    
    
%Nombre de couches

    L=5;

%Pas du r�seau
    
    Lb=0.9617;
    Lbmin=;Lbmax=;
    
%Epaisseurs

    d=0.7647;
    dmin=;dmax=;
    
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
    
    %Permittivit�s milieu ext�rieur incident
    
    EI_choice=2;
    
    %Permittivit�s milieu ext�rieur �mergent
    
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
            
            Fa=0.1979;% Grande Base trap�ze
            Famin=;Famax;
            Fb=0.65;% Petite Base trap�ze
            Fbmin=;Fbmax=;
            
        %R�seau en volume
        %----------------
                
            %Permittivit� moyenne
            
            %epsm(l)=;
            
            %Modulation de permittivit�
            
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
    
        
        
    

    

    
