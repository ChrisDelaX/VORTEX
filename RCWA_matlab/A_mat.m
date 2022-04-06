function[tp_inv_mat]=A_mat(prof_inv,N)

%Construction matrice eps inverse dite de TOEPLITZ INVERSE

%Passage de la permittivité inverse du réseau à sa transformée de Fourier

Etmp=ifft(prof_inv);%/(1024);

Etmp=fftshift(Etmp);

i=-N;

for k=1:2*N+1,
    
    p=-N;
    
    for j=1:2*N+1,
        
        tp_inv_mat(k,j)=Etmp(513+(i-p));
        
        p=p+1;
        
    end;
    
i=i+1;    
    
end;

