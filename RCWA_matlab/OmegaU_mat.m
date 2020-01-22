function[OmegaU]=OmegaU_mat(kx_mat,ky,E,N)

%Construction matrice Omega l

OmegaU=[ky^2*eye(2*N+1)+kx_mat^2-E];