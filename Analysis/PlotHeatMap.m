function PlotHeatMap(contacts,RFrange,title,save)
% Function to plot heat map. Inputs are contact matrix and title of
% figure

map=[ones(1000,1), linspace(0,1,1000)',linspace(0,1,1000)'];
RFlist=(RFrange(1):RFrange(2));
logcont=log(contacts./max(max(contacts)));
revlog=max(max(logcont))-logcont;
if length(revlog)>length(RFlist)
    RFlist=[RFlist (RFrange(2)+1)];
    RFrange(end)=RFrange(end)+1;
end

figure
pcolor(RFlist,RFlist,revlog);shading flat;axis image;colormap bone; axis ij;
colormap(map);
ticklocs=linspace(min(min(revlog)),max(max(revlog(revlog~=Inf))),5);
c=colorbar('Direction','reverse','Ticks',ticklocs,'TicksMode','manual',...
    'TickLabels',floor(linspace(0,min(min(logcont(logcont~=-Inf))),5)));
c.Label.String='ln($IF_{ij}$)';
c.Label.Interpreter='latex';
c.Label.FontSize=15;
xlabel('RF label','Interpreter','latex','FontSize',20);
ylabel('RF label','Interpreter','latex','FontSize',20);
if save
    saveas(gcf,title,'png')
end

end