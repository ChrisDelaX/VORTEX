function[tp_mat]=E_mat(prof,N)

%Construction matrice eps dite de TOEPLITZ

%Passage de la permittivité du réseau à sa transformée de Fourier

Etmp=ifft(prof);
Etmp=fftshift(Etmp);

%Etmp2=fft(prof)./sum(prof(1,:));
%Etmp-Etmp2

i=-N;

for k=1:2*N+1
    
    p=-N;
    
    for j=1:2*N+1
        
        tp_mat(k,j)=Etmp(513+(i-p));
        
        p=p+1;
        
    end
    
i=i+1;    
    
end