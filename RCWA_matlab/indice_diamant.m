
clear all
nlb=1e3;
for ib=1:nlb
    lb_min=2;%10.5;%
    lb_max=2.1;%.65;%
    lb(ib)=lb_min+(lb_max-lb_min)*(ib-1)/(nlb-1);
    A_d=1 ; B_d=0.3306 ; C_d=175.0 ; D_d= 4.3356 ; E_d=106.0 ;
    E2=sellmeier_diamant(lb(ib),A_d,B_d,C_d,D_d,E_d);
    n2(ib)=sqrt(E2);
end

figure
hold on
plot(lb,n2)

% plot(lb,polyval(polyfit(lb,n2,100),lb))
% 
% p = polyfit(lb,n2,2);
% lbi = .633
% nlbi = polyval(p,lbi)