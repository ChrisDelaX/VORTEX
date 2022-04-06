sansBAL

texto=strcat('pente_',num2str(pente/pi*180));
fn=strcat('fig_',texto,'.ps')
path='C:\Documents and Settings\Christianek\Bureau\RCWA_2008\data\2009-06-28\';
fn_new=[path fn];
fn1=strcat('fig_',texto,'.jpg')
fn1_new=[path fn1];

%print -deps -r300 
print(gcf, '-dps','-r300', fn_new )
print(gcf, '-djpeg','-r300', fn1_new )