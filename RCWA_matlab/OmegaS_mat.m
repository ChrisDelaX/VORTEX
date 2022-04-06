function[OmegaS]=OmegaS_mat(ky,A,B,N)

%Construction matrice Omega l

OmegaS=[ky^2*eye(2*N+1)+B*inv(A)];