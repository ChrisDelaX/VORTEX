%Programme principal
%-------------------
%-------------------
                
%Vecteur de polarisation incident
%--------------------------------

ux=cos(psi)*cos(theta)*cos(phi)-sin(psi)*sin(phi);
uy=cos(psi)*cos(theta)*sin(phi)+sin(psi)*cos(phi);
uz=-cos(psi)*sin(theta);
    
%Calcul Vecteur d'onde
%---------------------

k=2*pi/lb;
    
%Boucle sur "j"
%--------------

%Construction ky
%---------------

    ky=k*real(sqrt(EI))*sin(theta)*sin(phi);

    for j=1:2*N+1,
        jbis=j-(N+1);
        %Construction kx_vec
        %-------------------
        kx_vec(j)=k*(real(sqrt(EI))*sin(theta)*cos(phi)-jbis*(lb/Lb));   
        %Construction kIz_vec 
        %--------------------
        kIz_vec(j)=sqrt(k^2*EI-(kx_vec(j))^2-(ky)^2);
        %Construction kIIIz_vec 
        %----------------------
        kIIIz_vec(j)=sqrt(k^2*EIII-(kx_vec(j))^2-(ky)^2);
    end;%Fin de la boucle sur "j"

%Matrice kx
%----------

kx_mat=diag(kx_vec)/k;
kxoo=kx_mat(N+1,N+1);

%Matrice ky
%----------

ky_vec=ky*ones(2*N+1,1);
ky_mat=diag(ky_vec)/k;
kyoo=ky_mat(N+1,N+1);

%Matrice kIz
%-----------

kIz_mat=diag(kIz_vec)/k;

%Matrice kIIIz
%-------------

kIIIz_mat=diag(kIIIz_vec)/k;   

%Matrice delta
%-------------

delta=eye(2*N+1);
delta(N+1,N+1)=1;

%Boucle "l", le nombre de couches
%--------------------------------
%--------------------------------

%Méthode des matrices S
%----------------------

%Initialisation
%--------------

Tuu_l=eye(2*(2*N+1));
Rud_l=zeros(2*(2*N+1));
Rdu_l=zeros(2*(2*N+1));
Tdd_l=eye(2*(2*N+1));

%Calcul de W_lp1 initial
%-----------------------

F_lp1=[eye(2*N+1)                         zeros(2*N+1)                       eye(2*N+1)                       zeros(2*N+1)                  
       zeros(2*N+1)                       eye(2*N+1)                         zeros(2*N+1)                     eye(2*N+1)
       kx_mat*ky_mat/kIIIz_mat            (ky_mat^2+kIIIz_mat^2)/kIIIz_mat   -kx_mat*ky_mat/kIIIz_mat         -(ky_mat^2+kIIIz_mat^2)/kIIIz_mat
       -(kx_mat^2+kIIIz_mat^2)/kIIIz_mat  -kx_mat*ky_mat/kIIIz_mat           (kx_mat^2+kIIIz_mat^2)/kIIIz_mat  kx_mat*ky_mat/kIIIz_mat ];
  
X_lp1=[eye(2*(2*N+1))];
   
for l=1:L,

%Valeur et vecteurs propres de la matrice Omega
%----------------------------------------------

if var(prof(l,:))<1e-10
    prof_tmp=prof(l,1);clear E;E=eye(2*N+1)*prof_tmp;clear A;A=eye(2*N+1)*1/prof_tmp;
%     kx_vec
%     ky_vec
    Sigmatmp=diag([((-prof_tmp+(kx_vec/k).^2+(ky_vec'/k).^2)) ((-prof_tmp+(kx_vec/k).^2+(ky_vec'/k).^2))]);
    W=eye(2*(2*N+1));   
    
else

E=E_mat(prof(l,:),N);

A=A_mat(prof_inv(l,:),N);

B=B_mat(kx_mat,E,N,A);

D=D_mat(ky_mat,E,N,A);

Omega=Omega_mat(kx_mat,ky_mat,E,A,B,D,N);

[W,Sigmatmp]=eig(Omega);

end;

    %Matrice F
    %---------
    
        %Matrices V
        %----------
            
            %Matrice Sigma - vecteur sigma
            %-----------------------------
            
            Sigmatmp=diag(sqrt(Sigmatmp));
            
            %Test sur les valeurs propres
            %----------------------------
            
            for j=1:2*N+1,
                if (real(Sigmatmp(j))+imag(Sigmatmp(j)))>0
                Sigmatmp(j)=Sigmatmp(j);
                elseif (real(Sigmatmp(j))+imag(Sigmatmp(j)))<0
                Sigmatmp(j)=-Sigmatmp(j);
                end;
            end;            
            sigma=Sigmatmp;%clear Sigmatmp;                        
            Sigma=diag(sigma);
            
            %Matrices Q
            %----------
            
            Q11=kx_mat*ky_mat;
            Q12=inv(A)-ky_mat^2;
            Q21=kx_mat^2-E;
            Q22=-kx_mat*ky_mat;
            
            %Matrices W
            %----------
            
            W1=W(1:2*N+1,1:2*(2*N+1));
            W2=W(2*N+1+1:2*(2*N+1),1:2*(2*N+1));
            
            %Matrices V
            %----------

invSigma=inv(Sigma);		
            
            V1=(Q11*W1+Q12*W2)*invSigma;
            V2=(Q21*W1+Q22*W2)*invSigma;
        
    F11=W2;
    F12=W2;
    F21=W1;
    F22=W1;
    F31=sqrt(-1)*V2;
    F32=-sqrt(-1)*V2;
    F41=sqrt(-1)*V1;
    F42=-sqrt(-1)*V1;
            
    %Matrices X
    %----------
    
    X_l=X_lp1;     
    X_lp1=diag(exp(-k*sigma*d_t(l)));
    
    %Matrices T et F
    %---------------
    
    F_l=F_lp1;
    
    F_lp1=[W2          W2
           W1          W1
           sqrt(-1)*V2 -sqrt(-1)*V2
           sqrt(-1)*V1 -sqrt(-1)*V1];clear W1 W2 V1 V2
                 
    %T=inv(F_lp1)*F_l;
    
T=F_lp1\F_l;

    t11=T(1:2*(2*N+1),1:2*(2*N+1));
    t12=T(1:2*(2*N+1),2*(2*N+1)+1:2*2*(2*N+1));
    t21=T(2*(2*N+1)+1:2*2*(2*N+1),1:2*(2*N+1));
    t22=T(2*(2*N+1)+1:2*2*(2*N+1),2*(2*N+1)+1:2*2*(2*N+1));
    clear T;
    
    invt22=inv(t22);
    
    %Matrice S
    %---------
    
        %Matrice s
        %---------
          
        s=[t11-t12*invt22*t21  t12*invt22
           -invt22*t21         invt22];clear t11 t12 t21 t22 invt22;
        
            %Matrice s chapeau
            %-----------------
    
            s_chap_g=[eye(2*(2*N+1))    zeros(2*(2*N+1)) 
                      zeros(2*(2*N+1))  X_l];   
                      
            s_chap_d=[X_l               zeros(2*(2*N+1)) 
                      zeros(2*(2*N+1))  eye(2*(2*N+1))];   
                  
            clear X_l;
                  
            s_chap=s_chap_g*s*s_chap_d;clear s;
            
            tuu=s_chap(1:2*(2*N+1),1:2*(2*N+1));
            rud=s_chap(1:2*(2*N+1),2*(2*N+1)+1:2*2*(2*N+1));
            rdu=s_chap(2*(2*N+1)+1:2*2*(2*N+1),1:2*(2*N+1));
            tdd=s_chap(2*(2*N+1)+1:2*2*(2*N+1),2*(2*N+1)+1:2*2*(2*N+1));
            clear s s_chap_g s_chap_d;
            
    Tuu_lm1=Tuu_l;
    Rud_lm1=Rud_l;
    Rdu_lm1=Rdu_l;
    Tdd_lm1=Tdd_l;
    
    Tuu_l=tuu*inv(eye(2*(2*N+1))-Rud_lm1*rdu)*Tuu_lm1;
    Rud_l=rud+tuu*Rud_lm1*inv(eye(2*(2*N+1))-rdu*Rud_lm1)*tdd;
    Rdu_l=Rdu_lm1+Tdd_lm1*rdu*inv(eye(2*(2*N+1))-Rud_lm1*rdu)*Tuu_lm1;       
    Tdd_l=Tdd_lm1*inv(eye(2*(2*N+1))-rdu*Rud_lm1)*tdd;         
    clear Tuu_lm1 Rud_lm1 Rdu_lm1 Tdd_lm1;
    clear tuu rud rdu tdd;
    
end;

%Fin de la boucle à rebourd sur "l"

%Calcul de Sn
%------------

F_l=F_lp1;

F_lp1=[eye(2*N+1)                     zeros(2*N+1)                       eye(2*N+1)                     zeros(2*N+1)                       
       zeros(2*N+1)                   eye(2*N+1)                         zeros(2*N+1)                   eye(2*N+1)  
       kx_mat*ky_mat/kIz_mat          (ky_mat^2+kIz_mat^2)/kIz_mat       -kx_mat*ky_mat/kIz_mat         -(ky_mat^2+kIz_mat^2)/kIz_mat                                 
       -(kx_mat^2+kIz_mat^2)/kIz_mat  -kx_mat*ky_mat/kIz_mat             (kx_mat^2+kIz_mat^2)/kIz_mat   kx_mat*ky_mat/kIz_mat];  
       
    %T=inv(F_lp1)*F_l;
    
T=F_lp1\F_l;
    
    t11=T(1:2*(2*N+1),1:2*(2*N+1));
    t12=T(1:2*(2*N+1),2*(2*N+1)+1:2*2*(2*N+1));
    t21=T(2*(2*N+1)+1:2*2*(2*N+1),1:2*(2*N+1));
    t22=T(2*(2*N+1)+1:2*2*(2*N+1),2*(2*N+1)+1:2*2*(2*N+1));
    clear T;

    %Matrice S
    %---------
    
        %Matrice s
        %---------
          
        invt22=inv(t22);
        
        s=[t11-t12*inv(t22)*t21  t12*invt22
           -invt22*t21           invt22];clear t11 t12 t21 t22 invt22;
        
            %Matrice s chapeau
            %-----------------
    
            s_chap_g=[eye(2*(2*N+1))    zeros(2*(2*N+1)) 
                      zeros(2*(2*N+1))  X_lp1];   
                      
            s_chap_d=[X_lp1             zeros(2*(2*N+1)) 
                      zeros(2*(2*N+1))  eye(2*(2*N+1))];   
                  
            clear X_lp1;   
                  
            s_chap=s_chap_g*s*s_chap_d;clear s s_chap_g s_chap_d;
            
            tuu=s_chap(1:2*(2*N+1),1:2*(2*N+1));
            rud=s_chap(1:2*(2*N+1),2*(2*N+1)+1:2*2*(2*N+1));
            rdu=s_chap(2*(2*N+1)+1:2*2*(2*N+1),1:2*(2*N+1));
            tdd=s_chap(2*(2*N+1)+1:2*2*(2*N+1),2*(2*N+1)+1:2*2*(2*N+1));
            clear s_chap;
            
    Tuu_lm1=Tuu_l;
    Rud_lm1=Rud_l;
    Rdu_lm1=Rdu_l;
    Tdd_lm1=Tdd_l;
    
    Tuu_l=tuu*inv(eye(2*(2*N+1))-Rud_lm1*rdu)*Tuu_lm1;
    Rud_l=rud+tuu*Rud_lm1*inv(eye(2*(2*N+1))-rdu*Rud_lm1)*tdd;
    Rdu_l=Rdu_lm1+Tdd_lm1*rdu*inv(eye(2*(2*N+1))-Rud_lm1*rdu)*Tuu_lm1;       
    Tdd_l=Tdd_lm1*inv(eye(2*(2*N+1))-rdu*Rud_lm1)*tdd;    
    clear Tuu_lm1 Rud_lm1 Rdu_lm1 Tdd_lm1 tuu rud rdu tdd;
    
%Equation finale
%---------------
%---------------

S=[Tuu_l Rud_l
   Rdu_l Tdd_l];

fintmp=zeros(2*N+1,1);
fintmp(N+1)=1;

fin1=[zeros(2*N+1,1)
      zeros(2*N+1,1)
      ux*fintmp
      uy*fintmp];

%Solutions
%---------
%---------
    
    solut=S*[fin1];
 
%Resultats
%---------
%---------

    %Amplitudes dans la base  xyz
    %----------------------------

    Rx=solut(1:2*N+1);
    Ry=solut((2*N+1)+1:2*(2*N+1));

    Tx=solut(2*(2*N+1)+1:3*(2*N+1));
    Ty=solut(3*(2*N+1)+1:4*(2*N+1));

    Rz=(-kx_vec'.*Rx-ky_vec.*Ry)./kIz_vec';
    Tz=(-kx_vec'.*Tx-ky_vec.*Ty)./kIIIz_vec';

    %Amplitudes dans la base s(TE) p(TM)
    %-----------------------------------

    for j=1:2*N+1,  
        if ky==0    
        Phi(j)=0;
        else
        Phi(j)=atan(ky./kx_vec(j));  
        end;
    end;

    Rs=cos(Phi').*Ry-sin(Phi').*Rx;
    Rp=(1/k).*(cos(Phi').*(kIz_vec'.*Rx-kx_vec'.*Rz)-sin(Phi').*(ky_vec.*Rz-kIz_vec'.*Ry));

    Ts=(cos(Phi').*Ty-sin(Phi').*Tx);
    Tp=(1/k)*(cos(Phi').*(kIIIz_vec'.*Tx-kx_vec'.*Tz)-sin(Phi').*(ky_vec.*Tz-kIIIz_vec'.*Ty));

    %Efficacités de diffraction
    %--------------------------
    
    nIc=k*sqrt(EI)*cos(theta);
    
    %nIIIc=k*sqrt(EIII)*cos(theta);
    
    DER_s = (Rs.*conj(Rs).*real(kIz_vec/nIc)')';
    DER_p = (Rp.*conj(Rp).*real((kIz_vec/EI)/nIc)')';
    
    DET_s = (Ts.*conj(Ts).*real(kIIIz_vec/nIc)')';
    DET_p = (Tp.*conj(Tp).*real((kIIIz_vec/EIII)/nIc)')';
 
    test=1-sum(DER_s+DER_p+DET_s+DET_p);
    
    %Test conservation énergie algo
    if abs(test) > 1e-10    
          test ; 
    end;
   
    
    phi_s_T=atan2(imag(Ts(N+1)),real(Ts(N+1)));
    phi_s_R=atan(imag(Rs(N+1))/real(Rs(N+1)));
    phi_p_T=atan2(imag(Tp(N+1)),real(Tp(N+1)));
    phi_p_R=atan(imag(Rp(N+1))/real(Rp(N+1)));
    
    dphi_sp_T=((atan2(imag(Tp(N+1)),real(Tp(N+1)))-atan2(imag(Ts(N+1)),real(Ts(N+1))))); 
    dphi_sp_R=(atan2(imag(Rp(N+1)),real(Rp(N+1)))-atan2(imag(Rs(N+1)),real(Rs(N+1))));
    
    %dphi_sp(i)=-(angle(Tp(N+1))-angle(Ts(N+1)));
    
    
%FIN PROGRAMME MAIN