function plotIFTriangles(serIF,parIF,rangeRF,title,saveornot)
% Parallel IF on top
pIF=parIF./max(max(parIF));
sIF=serIF./max(max(serIF));
totsize=length(sIF);
newMap=zeros(totsize);
for i=1:totsize
    for j=1:totsize
        if i>=j
            newMap(i,j)=sIF(i,j);
        else
            newMap(i,j)=pIF(i,j);
        end
    end
end

PlotHeatMap(newMap,rangeRF,title,saveornot);
end