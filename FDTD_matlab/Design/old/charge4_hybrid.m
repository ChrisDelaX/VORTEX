clear all
close all

periodMax = 1.42;%1;%
alph = 3;%   %slope in degrees
depth = 5.6;% %grating depth
F_max = 0.51;
Lb_min = 2*depth*tan(deg2rad(alph))/(1-F_max);
%Lb_min=periodMax*1e-4;%sqrt(.5);%0.8;

Fmoy=.45;%.48;%.1;%
periodRat = Lb_min/periodMax;%.5;%0.9;%.5;%   % ratio Lb_min/periodMax

rMax=40;%10;%5;%1000;%
nbPts=20;  % for one period Lb
nbStep = 3.5;  % offset at each level

% Figures
% -------
newFig
grid off
set(gca,'visible','off')
axis equal

tic
depart_tmp=now;
depart=datestr(depart_tmp)

i=0;
j=0;
periodStep = periodMax;
rExt=periodStep*(1+Fmoy)/4;
rInt=periodStep*(1-Fmoy)/4;
rMin = 0;
sinMax = 1;
sinMin = periodRat^(1/2);
aMax0 = 2*asin(sinMax);
aMin0 = 2*asin(sinMin);
while 2*rInt*sinMin < rMax
    
    while 2*rInt*sinMin < rMax
        
        if 2*rExt*sinMax > rMin

            % rExt - aMax
            % -----------
            if 2*rExt*sinMax <= rMax
                aMax = aMax0;
            elseif 2*rExt*sinMin < rMax
                aMax = 2*asin(rMax/(2*rExt));
            else
                aMax = aMin0;
                rExt = rMax/(2*sinMin);
            end
            % rExt - aMin
            % -----------
            if 2*rExt*sinMin >= rMin
                aMin = aMin0;
            else
                aMin = 2*asin(rMin/(2*rInt));
            end
            nbTh = ceil((aMax-aMin)*rExt/periodMax*nbPts+1);
            TH1 = linspace(aMin,aMax,nbTh);
            R1 = rExt*ones(1,nbTh);
            [X1,Y1] = pol2cart(TH1,R1);
            xx = -X1+rExt;
            yy = Y1;
            if 2*rInt*sinMax >= rMax %%%%%%
                aMax = 2*asin(rMax/(2*rInt));
            else

                aMax = aMax0;
            end

            % rInt - aMax
            % -----------
            if 2*rInt*sinMax <= rMax
                aMax = aMax0;
                if 2*rExt*sinMax > rMax
                    xx(end+1) = rMax*sinMax;
                    yy(end+1) = rMax*cos(asin(sinMax));
                end
            else
                aMax = 2*asin(rMax/(2*rInt));%%
            end
            
            % rInt - aMin
            % -----------
            if 2*rInt*sinMin >= rMin
                aMin = aMin0;
            elseif 2*rInt*sinMax > rMin
                aMin = 2*asin(rMin/(2*rInt));
            else
                aMin = aMax0;
                rInt = rMin/(2*sinMax);
            end
            nbTh = ceil((aMax-aMin)*rInt/periodMax*nbPts+1);
            TH2 = linspace(aMax,aMin,nbTh);
            R2 = rInt*ones(1,nbTh);
            [X2,Y2] = pol2cart(TH2,R2);
            %X2 = -X2+rInt;
            xx(end+1:end+nbTh) = -X2+rInt;
            yy(end+1:end+nbTh) = Y2;
%             if 2*rExt*sinMin >= rMin && 2*rInt*sinMin < rMin
%                 %aMax = 2*asin(rMax/(2*rInt));
%                 x(end+1)=rMin*sinMin;
%                 y(end+1)=rMin*cos(asin(sinMin));
%             end
            patch(xx,yy,mygrey,'edgecolor',mygrey);
            patch(-xx,yy,mygrey,'edgecolor',mygrey);
            patch(xx,-yy,mygrey,'edgecolor',mygrey);
            patch(-xx,-yy,mygrey,'edgecolor',mygrey);
        end
        
        i=i+1;
        rExt = periodStep*(1+Fmoy+2*i)/4;
        rInt = periodStep*(1-Fmoy+2*i)/4;
        
    end
    %increments
    j=j+1;
    i=0;%5*j;%
    periodStep = periodMax/periodRat^(j);
    rExt = periodStep*(1+Fmoy+2*i)/4;
    rInt = periodStep*(1-Fmoy+2*i)/4;
    rMin = periodStep*(1-Fmoy+2*nbStep*j)/4;
    sinMax = periodRat^(j/2);
    sinMin = periodRat^((j+1)/2);
    aMax0 = 2*asin(sinMax);
    aMin0 = 2*asin(sinMin);

    
end
nb_dc=j-1  % number of discontinuities


%depart
%fin
toc
%datestr(now)


% Graphics
la_tmp=.05;%
line([-rMax rMax]*(1+la_tmp),[0 0],'Color',mygreen,'LineStyle','--')
text(rMax*(1+1.5*la_tmp),0,'$0^{\circ}$')
text(-rMax*(1+5*la_tmp),0,'$180^{\circ}$')
line([-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),[-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(rMax/sqrt(2)*(1+la_tmp),rMax/sqrt(2)*(1+2*la_tmp),'$45^{\circ}$')
text(-rMax/sqrt(2)*(1+4*la_tmp),-rMax/sqrt(2)*(1+3*la_tmp),'$225^{\circ}$')
line([0 0],[-rMax rMax]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(-rMax*la_tmp,rMax*(1+2*la_tmp),'$90^{\circ}$')
text(-rMax*la_tmp,-rMax*(1+2*la_tmp),'$270^{\circ}$')
line([rMax/sqrt(2) -rMax/sqrt(2)]*(1+la_tmp),[-rMax/sqrt(2) rMax/sqrt(2)]*(1+la_tmp),'Color',mygreen,'LineStyle','--')
text(-rMax/sqrt(2)*(1+4*la_tmp),rMax/sqrt(2)*(1+2.5*la_tmp),'$135^{\circ}$')
text(rMax/sqrt(2)*(1+la_tmp),-rMax/sqrt(2)*(1+2.5*la_tmp),'$315^{\circ}$')
axis equal
xlabel('$\mu m$')
ylabel('$\mu m$')
set(gca,'xLim',[-rMax rMax]*(1+8*la_tmp))
set(gca,'ylim',[-rMax rMax]*(1+6*la_tmp))
str = sprintf('L band $$\\rightarrow \\Lambda_{max}= %3.3g \\mu m, \\:\\: h = %3.3g \\mu m \\rightarrow \\frac{\\Lambda_{min}}{\\Lambda_{max}}= %.2g \\rightarrow F = %.2g $$',periodMax,depth,periodRat,Fmoy);
title(str,'Fontname',fnz,'FontSize',fsz,'FontWeight','bold')
%tick2latex

%save
print('-depsc2',sprintf('curves_sketch_Lbrat=%3.2f_F=%3.2f_r=%d.eps',periodRat,Fmoy,rMax), '-r300');








