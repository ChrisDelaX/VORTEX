function [ A,A_rhs] = nma_generate_A_and_ARHS_for_2D_diffusion_Neumman(n,D,k,h)
% generate the A and B matrices used by the ADI solver.

%nma_generate_A_and_ARHS_for_2D_diffusion_Neumman
% This function generate the A and B matrices to use
% for solving problem 2, HW 2.
% Used by the ADI solver.
%
% INPUT:
%   n size of matrix (1D size)
%   D diffusion constant
%   k time step
%   h space_step
%
% OUTPUT
%   A: represent the Lx or Ly operators, but applied to all the
%      rows or columns at once.
%   B: the RHS matrix used for the explicit part of the ADI solver. See
%      my HW report for more information.
%

% Generate the sparse A matrix

r  = D*k/(2*h^2);
e  = ones(n,1);
B  = [-r*e (1+2*r)*e -r*e];
Lx = spdiags(B,[-1 0 1],n,n);
Ix = speye(n);
A  = kron(Ix,Lx);

% adjust A, see HW report, some terms on the diagonal need to be
% 1+r, not 1+2*r

pos = 1:n:n^2;
for i = 1:length(pos)
   A(pos(i),pos(i))=1+r;
end

pos = n:n:n^2;

for i=1:length(pos)
   A(pos(i),pos(i))=1+r;
end

% Now generate the B matrix

e = ones(n^2,1);
B = [r*e (1-r)*e r*e];
A_rhs = spdiags(B,[-n 0 n],n^2,n^2);

%adjust matrix, see HW report
pos = n+1:n^2-n;
for i=1:length(pos)
   A_rhs(pos(i),pos(i))=1-2*r;
end

end
