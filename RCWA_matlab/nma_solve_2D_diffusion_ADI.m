function [u,X,Y,u_steady_state]=nma_solve_2D_diffusion_ADI(...
   A,A_RHS,h,k,D,max_t,ic,SHOW_SUM_AT_EACH_STEP)
%Solves 2D diffusion PDE $u_t = D  (u_xx+u_yy)$ on unit square Neumann BC using cell centered grid points and ADI schem

%nma_solve_2D_diffusion_ADI Solves the 2D diffusion PDE problem
%u_t = D * (u_xx+u_yy)  on unit square (0,1)x(0,1) assuming
%Neumann boundary conditions. Uses cell centered grid points
%and used ADI scheme.
%
%INPUT
%
%A: Sparse matrix that represent Lx or Ly to use for ADI
%A_RHS: Sparse matrix, to use for RHS of the ADI (explicit)
%h:  The space size, for cell centered grid
%k: time step size
%D: diffusion constant
%max_t: maximum time to run solver for
%ic: initial data in a grid.

%SHOW_SUM_AT_EACH_STEP:  % flag to tell if to show sum of solution
%
%OUTPUT
%u: The solution u
%X: X matrix from meshgrid
%Y: Y matrix from meshgrid
%u_steady_state: The steady state solution
%
%EXAMPLE CALL:
% D=0.1;
% N=81;
% h=1/N;
% k=h;
% ic_v = @(X,Y) exp(-100*(X.^2+Y.^2));
% [X,Y] = meshgrid(h/2:h:1-h/2,h/2:h:1-h/2);
% ic = ic_v(X,Y);
% [A,A_rhs]= ...
%   nma_generate_A_and_ARHS_for_2D_diffusion_Neumman(N,D,k,h);
% max_t=1;
% [u,X,Y,~]=nma_solve_2D_diffusion_ADI(...
%   A,A_rhs,h,k,D,max_t,ic,false)
% mesh(X,Y,u);
%
%Nasser M. Abbasi
%Math 228B, UC Davis, Winter 2011


[X,Y]  = meshgrid(h/2:h:1-(h/2),h/2:h:1-(h/2)); % coordinates
u_mean = mean(mean(ic)); %quad2d(ic,h/2,1-(h/2),h/2,1-(h/2));
grid_size = size(ic,1); %cell centered grid
u              = ic;  %(X,Y); % initial U
u_steady_state = zeros(size(u));
u_steady_state(:,:) = u_mean;

%for plotting only, to set z limit
u_max  = max(max(u));
u_min  = min(min(u));

% initialize the loop
current_time_step = 0;
number_of_time_steps = round(max_t/k);
display_plots = false; %flag to display or not each frame

if display_plots
   fig_handle = figure();
end

if SHOW_SUM_AT_EACH_STEP
   fprintf('\tcurrent_time\t\tsum(U(current_time)\n');
end

while current_time_step < number_of_time_steps
   % solve for U*
   tmp = reshape(flipud(u)',grid_size^2,1);
   sol = A\(A_RHS*tmp);
   u   = flipud(reshape(sol,grid_size,grid_size)');
   
   % solve for U_n+1
   tmp = reshape(flipud(u),grid_size^2,1);
   sol = A\(A_RHS*tmp);
   u   = flipud(reshape(sol,grid_size,grid_size));
   
   update_sum(SHOW_SUM_AT_EACH_STEP,current_time_step,k,u);
   if display_plots
      update_plot(X,Y,u,fig_handle,u_min,u_max,current_time_step,k,D,h,...
         u_mean,grid_size,u_steady_state);
   end
   
   %update current time and check if reached end of time
   current_time_step = current_time_step + 1;
end

if display_plots
   close(fig_handle);
end

end

%-----------------------------
function update_plot(X,Y,u,fig_handle,u_min,u_max,current_time_step,k,D,h,...
   u_mean,grid_size,u_steady_state)


set(0,'CurrentFigure',fig_handle);
surf(X,Y,u);
colormap cool;
title(sprintf('solution at time = %1.3f, D=%.3f, N=%d\nsteady state=%1.4f, h=%1.5f',...
   current_time_step*k,D,grid_size,u_mean,h));
zlim([u_min u_max]);

hold on;
mesh(X,Y,u_steady_state);
zlim([u_min u_max]);
drawnow();
hold off;

end
%---------------------------------
function update_sum(SHOW_SUM_AT_EACH_STEP,current_time_step,k,u)

% display the sum if asked, this is to verify discrete
% conservation, since Neumann boundary conditions
if SHOW_SUM_AT_EACH_STEP
   fprintf('\t%.5f\t\t\t%5.5f\n',current_time_step*k,sum(sum(u)));
end

end