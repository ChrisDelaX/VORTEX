close all
f=figure('FileName','Diamant');
set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
set(f,'color',[1 1 1])
hold on
%grid on
xlabel('F_{top}','FontSize',20,'FontWeight','bold')
ylabel('h (µm)','FontSize',20,'FontWeight','bold')
%set(gca,'YDir','reverse')
set(gca,'FontSize',20,'FontWeight','bold')

% données
a = 2.75;
tan_a = tan(a/180*pi);
Lb = 4.6;
h = 10:.1:15;

% optimum F
hmoy = 13.6875;
Feqmoy = 0.4570;
Fmoy = Feqmoy / (1+tan_a/Lb*hmoy);

% droite pour Feq
b1 = 28.9017;
b0 = hmoy-Feqmoy*b1;

Feq = (h-b0)./b1;
%plot(Feq,h,'k','Linewidth',2)

F = Feq ./ (1 + tan_a/Lb.*h);
plot(F,h,':r','Linewidth',2)

% droite pour F
b1 = 37.594;
b0 = hmoy-Fmoy*b1;

Flin = (h-b0)./b1;
plot(Flin,h,'b','Linewidth',1.5)

legend('h = f(F)','linear regression')

% axes
xmin = 0.3815;
xmax = 0.4214;
ymin = 13;
ymax = 14.5;
axis([xmin xmax ymin ymax])
