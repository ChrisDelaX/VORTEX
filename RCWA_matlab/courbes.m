close all
figure; plot(Xparam,nulldpt.^-1)
%title('Diamant/SiO2 (N=8, d_{AR}=0.5µm) pente 1° (L=16)')
title('Diamant/SiO2 (N=8, d_{AR}=0.5µm) pente 1° (L=16)')
xlabel('Grating depth (microns)');
ylabel('Rejection ratio');
figure; plot(Xparam,Yparam(:,1))
title('Diamant/SiO2 (N=8, d_{AR}=0.5µm) pente 1° (L=16)')
xlabel('Grating depth (microns)');
ylabel('Grating period (microns)');
figure; plot(Xparam,Yparam(:,2))
title('Diamant/SiO2 (N=8, d_{AR}=0.5µm) pente 1° (L=16)')
xlabel('Grating depth (microns)');
ylabel('Grating filling factor');