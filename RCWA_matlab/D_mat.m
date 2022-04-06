function[Dtmp]=D_mat(ky_mat,E,N,A)

%Construction matrice D

Dtmp=ky_mat*inv(E)*ky_mat-eye(2*N+1);