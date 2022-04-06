%INPUT
%
%A: Sparse matrix that represent Lx or Ly to use for ADI
%A_RHS: Sparse matrix, to use for RHS of the ADI (explicit)
%h:  The space size, for cell centered grid
%k: time step size
%D: diffusion constant
%max_t: maximum time to run solver for
%ic: initial data in a grid.


D=0.1;
N=81;
h=1/(N);
k=h;
ic_v = @(X,Y) exp(-100*(X.^2+Y.^2));
[X,Y] = meshgrid(-1+h/2:2*h:1-h/2,-1+h/2:2*h:1-h/2);
ic = ic_v(X,Y);
[A,A_rhs]= ...
  nma_generate_A_and_ARHS_for_2D_diffusion_Neumman(N,D,k,h);
max_t=1;
[u,X,Y,~]=nma_solve_2D_diffusion_ADI(...
  A,A_rhs,h,k,D,max_t,ic,false)
mesh(X,Y,u);