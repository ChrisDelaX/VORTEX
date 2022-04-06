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
            scan4=10;
            mesh4=10;
                scan5=8;
                mesh5=100;

%Param�tres de calculs
%---------------------

    %Troncature
    
        %En X (N donc 2N+1 ordres au total)  
    
        N=10;
              
%Param�tres du r�seau
%--------------------    
    
%Pas du r�seau
    
    Lb=1.1211     ;
    Lbmin=0.45;Lbmax=0.5;
    
%Epaisseurs

    d=0.7166        ;
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
    
    %Permittivit�s milieu ext�rieur incident
    
    EI_choice=6;
    
    %Permittivit�s milieu ext�rieur �mergent
    
    EIII_choice=1;
    
    %Permittivit�s du r�seau
    %-----------------------
    
        %R�seau en relief de surface
        %---------------------------
        
            %Permittivit� milieu 1
            
            E1_choice=1;
            
            %Permittivit� milieu 2
            
            E2_choice=6;
                                   
            %Facteurs de remplissage
            
            Fa=0.31 ;% Grande Base trap�ze
            Famin=0.65;Famax=0.8;
            
            Fb=0.25;% Petite Base trap�ze
                    if Fb==Fa  
                    Fbmin=Famin;Fbmax=Fbmin;
                    else
                    Fbmin=0.5;Fbmax=0.9;
                    end;
            
        %Permittivit�s pseudo A/R
        %------------------------
           
            E_AR_choice=7;
            
            %Epaisseur de la couche AR
            
            d_AR=0.4231;
            d_AR_min=0.320;d_AR_max=0.42;
            
        %Nano-roughness
        %--------------
        
            r=0.005;
            
            %Epaisseur en surface
            
            d_R=r;
             %d=d-2*d_R; 
        %R�seau en volume
        %----------------
                
            %Permittivit� moyenne
            
            %epsm(l)=;
            
            %Modulation de permittivit�
            
            %modmin(l)=;
            %modmax(l)=;
            %mod(l)=;                       
            
%Choix du type de r�seau et Nombre de couches de discr�tisation
%--------------------------------------------------------------
    
    %R�seau rectangulaire/trap�zoidal simple sans Pseudo-A/R : 1
    %R�seau rectangulaire/trap�zoidal simple avec Pseudo-A/R : 2
    %R�seau trap�zoidal avec Pseudo-A/R adh�rent au parois : 3
    %R�seau trap�zoidal avec nano-roughness : 4
    
    ar=1;%par d�faut ar=1
       
    %Nombre de couche hors A/R
    L_wo_AR=32;
    
    %Nombre de couche hors nano-roughness
    L_wo_R=round((d/r)/8);
    
    if ar==2
        if Fb==Fa
        ar=2;%3 couches
        else
        %TBD!!!
        ar=3;%au moins 5 couches, pour la pr�cision de la discr�tisation
        end;
    end;
    
%Onde incidente 
%--------------

    %Longueur d'onde
    
    lb=4;
    lbmin=6;lbmax=11;nlb=10;

    %Angle incidence non conique
      
    theta=0.6048;
    %theta=39.5*pi/180;
    thetamin=35*pi/180;thetamax=50*pi/180;
   
    %Angle incidence conique

    phi=0*pi/180;
    phimin=0;phimax=0;
        
    %Angle polarisation
    
    psi=pi/4;
    psimin=0;psimax=0;
    
        
        
    

    

    
