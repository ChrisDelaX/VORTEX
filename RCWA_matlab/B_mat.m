function[Btmp]=B_mat(kx_mat,E,N,A)

%Construction matrice B

Btmp=kx_mat*inv(E)*kx_mat-eye(2*N+1);
