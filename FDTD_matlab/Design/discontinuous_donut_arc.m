function discontinuous_donut_arc(x0,y0,r1,r2,width1,width2,name)
%writes to a text file an AutoCAD command line script that will draw
%circular polyline shapes similar to a long rectangle twisted around the
%center until the ends meet. x0,y0: center coordinates. r1,r2: radius of the
%center circle and the cut off radius. width1,width2: distance between
%bands and band thickness. name: optional argument to add a custom name to
%the file. All lengths and coordinates in µm.

%change µm to mm
%x0=x0/1000;
%y0=y0/1000;
%width1=width1/1000;
%width2=width2/1000;
%r1=r1/1000;
%r2=r2/1000;

r=r1+width1; %starting radius

if nargin==7
    fid=fopen(['discontinuous_donut_' name '.scr'],'w');
else
    fid=fopen('discontinuous_donut.scr','w');
end
fprintf(fid,'-osnap\noff\n');

%draw center circle
fprintf(fid,'pline\n');
fprintf(fid,[num2str(x0+r1,8) ',' num2str(y0,8) '\n' 'a\nce\n' num2str(x0) ',' num2str(y0) '\na\n180\nclose\n']);

%draw discontinuous donuts
tic
while r<r2 %loop until cut off radius
    fprintf(fid,'pline\n'); %start polyline object

    fprintf(fid,[num2str(x0+r,8) ',' num2str(y0,8) '\n' 'a\nce\n' num2str(x0) ',' num2str(y0) '\na\n180\n' num2str(x0+r,8) ',' num2str(y0,8) '\nl\n']);
    r=r+width2; %increase to outer radius
    fprintf(fid,[num2str(x0+r,8) ',' num2str(y0,8) '\n' 'a\nce\n' num2str(x0) ',' num2str(y0) '\na\n-180\n' num2str(x0+r,8) ',' num2str(y0,8) '\nl\nclose\n']);
    r=r+width1; %increase radius for next interation
end
fprintf(fid,'-osnap\nend,cen,int,ext\n'); %reset object snap

fclose(fid); %close file
toc
disp('end')
disp(num2str(datestr(now)))
