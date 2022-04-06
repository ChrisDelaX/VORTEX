function[Omega]=Omega_mat(kx_mat,ky_mat,E,A,B,D,N)

%Construction matrice Omega l

%Omega=[ky_mat^2+B*inv(A)                     zeros(2*N+1)
%       zeros(2*N+1)   ky_mat^2+kx_mat^2-E ];

Omega=[kx_mat^2+D*E                      ky_mat*(inv(E)*kx_mat*inv(A)-kx_mat)
       kx_mat*(inv(E)*ky_mat*E-ky_mat)   ky_mat^2+B*inv(A) ];