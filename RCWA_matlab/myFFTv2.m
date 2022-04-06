clear all
close all

% Input paramters
% ---------------
N = 2^9;%10;%
%N2 = 2*N+1;
Rmax = 4;%2;%
D = 1; % lambda/D
lim = 8;%15;
mask = 1;  % 0=FQPM ; 1=AGPM
Lyot_rat = 0.9; % Lyot-stop ratio
switch mask
    % FQPM
    case 0
        maskstr = 'FQPM';
    % AGPM 
    case 1
        maskstr = 'AGPM';
end

% Entrance pupil
% --------------
P = zeros(N);
xy = linspace(-Rmax,Rmax,N);
[X,Y] = meshgrid(xy);
[TH,R] = cart2pol(X,Y);
P(X.^2 + Y.^2 <= (D/2)^2) = 1;
% newFig
% imagesc(xy,xy,abs(P))
% colormap('gray')
% hC=colorbar;
% axis image
% title('{\bf Entrance pupil}')
% xlabel('$D/2$')
% ylabel('$D/2$')
% tick2latex
% print('-depsc2',sprintf('fft_%s_N=%d_r=%d_entrance.eps',maskstr,N,Rmax), '-r300');


% % Tip-tilt
% % --------
% Tx=0.5;%0;%
% Ty=0;
% P=P.*exp(-1i*pi*(X.*Tx+Y.*Ty));

% % Companion
% % ---------
% Rc=1;%3/(D/2); 
% THc=20 ; Xc=Rc*cos(deg2rad(THc)); Yc=Rc*sin(deg2rad(THc));
% DRmag=15;
% FR=10^(-DRmag/2.5);
% fprintf('\nCompanion Separation = %3.2f lbd/D \n',Rc)
% fprintf('\nCompanion Flux = %3.2s \n',FR)
% Pc=sqrt(FR).*P.*exp(-1i*pi*(X*Xc+Y*Yc));
% P=P+Pc;


% Focal plane
% -----------
F = ((ifftshift(P)));
%F = fftshift(fft2(P));
F=F./sqrt((max(max(abs(F).^2))));
Xmax = (N/2) /(2*Rmax);     % lambda/D = 1/(2*Rmax) = 2*Xmax/N
xyF = linspace(-Xmax,Xmax,N);
% newFig
% imagesc(xyF,xyF,abs(F))
% colormap(pink)
% hC=colorbar;
% axis ([-lim lim -lim lim])
% title('{\bf Focal plane before mask}')
% xlabel('$\lambda/D$')
% ylabel('$\lambda/D$')
% tick2latex
% print('-depsc2',sprintf('fft_%s_N=%d_r=%d_before.eps',maskstr,N,Rmax), '-r300');


% Phase mask
% ----------
M = ones(N);
switch mask
    % FQPM
    case 0
        M(X>0 & Y>0)=exp(1i*pi);
        M(X<0 & Y<0)=exp(1i*pi);
        %M(N+1,:)=0;
        %M(:,N+1)=0;
        %M(N:N+2,N:N+2)=0;
    % AGPM
    case 1
        M=exp(2*1i*TH);
end
newFig
imagesc(xyF,xyF,real(M))
colormap('gray')
hC=colorbar;
axis image
title('{\bf Phase mask}')
xlabel('$\lambda/D$')
ylabel('$\lambda/D$')
tick2latex
print('-depsc2',sprintf('fft_%s_N=%d_r=%d_PM.eps',maskstr,N,Rmax), '-r300');


% Focal plane after the mask
% --------------------------
F2 = F.*M;
Phase = angle(F2);
newFig
imagesc(xyF,xyF,Phase)
colormap('gray')
hC=colorbar;
axis ([-lim lim -lim lim])
title('{\bf Phase after the mask}')
xlabel('$\lambda/D$')
ylabel('$\lambda/D$')
tick2latex
print('-depsc2',sprintf('fft_%s_N=%d_r=%d_Phase_after.eps',maskstr,N,Rmax), '-r300');

break


% Relayed pupil plane
% -------------------
P2 = (fft2(F2))./(N^2);
% newFig
% imagesc(xy,xy,abs(P2))
% colormap('gray')
% hC=colorbar;
% axis image
% title('{\bf Relayed pupil plane}')
% xlabel('$D/2$')
% ylabel('$D/2$')
% tick2latex
% print('-depsc2',sprintf('fft_%s_N=%d_r=%d_Relayed.eps',maskstr,N,Rmax), '-r300');


% Lyot stop
% ---------
L=ones(N);
L(R>Lyot_rat*D/2)=0;
P3=P2.*L;
% newFig
% imagesc(xy,xy,abs(P3))
% colormap('gray')
% hC=colorbar;
% axis image
% title('{\bf Lyot stop}')
% xlabel('$D/2$')
% ylabel('$D/2$')
% tick2latex
% print('-depsc2',sprintf('fft_%s_N=%d_r=%d_Lyot.eps',maskstr,N,Rmax), '-r300');


% Intensity on the detector
% -------------------------
I = abs((fft2(P3))).^2;
newFig
imagesc(xyF,xyF,I)
colormap(pink)
hC=colorbar;
axis ([-lim lim -lim lim])
title('{\bf Coronagraphic image}')
xlabel('$\lambda/D$')
ylabel('$\lambda/D$')
tick2latex
print('-depsc2',sprintf('fft_%s_N=%d_r=%d_corona.eps',maskstr,N,Rmax), '-r300');
break
% log scale
newFig
imagesc(xyF,xyF,log10(I))
colormap(pink)
hC=colorbar;
axis ([-lim lim -lim lim])
title('{\bf Coronagraphic image (log$_{10}$)}')
xlabel('$\lambda/D$')
ylabel('$\lambda/D$')
tick2latex
print('-depsc2',sprintf('fft_%s_N=%d_r=%d_corona_log.eps',maskstr,N,Rmax), '-r300');



