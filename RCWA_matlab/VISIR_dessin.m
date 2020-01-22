
[vi,i]=min(nulldpt,[],1);
[min_nulldpt,j]=min(vi);

f=figure('FileName','Diamant');
set(f,'color',[1 1 1])
hold on
%grid on


hS = pcolor(Xparam,Yparam,log10(nulldpt'));
shading interp
hC = colorbar;
set(hC,'FontSize',16,'FontWeight','bold')
%set(hC,'YScale','log')



[b]=plot(Xparam(i(j)),Yparam(j),'Marker','+');
set(b,'Linewidth',2,'color',[1 0 0])
xlabel('F','FontSize',16,'FontWeight','bold')
ylabel('d','FontSize',16,'FontWeight','bold')
%title('Mean Null Depth (11-13.2µm)  --  \Lambda = 4.6 µm / \alpha = 4°','FontSize',16,'FontWeight','bold')
set(gca,'YDir','reverse')
set(gca,'FontSize',16)
set(gca,'FontWeight','bold')


