function [runsum, rhos, rhoktemp, weightktemp]=CompareIFs(filename1,filename2,outputfileprefix)
% filename1 has IF outputs of serial run
% filename2 has IF outputs of parallel name
% outputfileprefix for where to save the comparison metrics and figure
% normfactor=1e-5; % numbers in files are actually IF/1e-5

% open and read serial IF output
fid=fopen(filename1,'r');
% First line has the starting and ending rows/columns
firstline=fgetl(fid);
bounds=str2num(firstline(3:end));
minRF=min(bounds);maxRF=max(bounds);
numRF=maxRF-minRF+1;
serialIF=zeros(numRF); % empty matrix to hold serialIF
parIF=zeros(numRF); % empty matrix to hold parallelIF

% collect rest of data
TempDat=textscan(fid,'%f %s %d %s %d');
fclose(fid);
IFDat=TempDat{1};
RF1=TempDat{3}-minRF+1;RF2=TempDat{5}-minRF+1;
% Populate serial matrix
sersum=sum(IFDat);
for i=1:length(IFDat)
    serialIF(RF1(i),RF2(i))=IFDat(i);
    serialIF(RF2(i),RF1(i))=IFDat(i);
end

% Handle parallel output
fid2=fopen(filename2,'r');
firstline2=fgetl(fid);
bounds2=str2num(firstline2(3:end));
minRF2=min(bounds2);maxRF2=max(bounds2);
if (minRF2~=minRF) || (maxRF~=maxRF2)
    error("TSV files do not have same RF bounds\n");
end
% collect rest of data
TempDat2=textscan(fid,'%f %s %d %s %d');
fclose(fid2);
IFDat2=TempDat2{1};
RF12=TempDat2{3}-minRF2+1;RF22=TempDat2{5}-minRF2+1;
% Populate parallel matrix
parsum=sum(IFDat2);
for i=1:length(IFDat2)
    parIF(RF12(i),RF22(i))=IFDat2(i);
    parIF(RF22(i),RF12(i))=IFDat2(i);
end

% squared residuals
ressquared=((parIF./parsum-serialIF./sersum)).^2;
% sum of squared residuals for upper triangular
runsum=0;
for i=1:numRF
    for j=i:numRF
        runsum=runsum+ressquared(i,j);
    end
end

% Compare with SCC, no smoothing
[rhos, rhoktemp, weightktemp]=CompSCC(serialIF,parIF,numRF-1,0);

% plot both 
plotIFTriangles(serialIF,parIF,[minRF maxRF],[outputfileprefix '_TriPlot'],1);

fid3=fopen([outputfileprefix '_Comparison.txt'],'w');
fprintf(fid3,'Sum of squared residuals (re-normalized)clc: %4.3e \n',runsum);
fprintf(fid3,'Stratified correlation coefficient: %4.3e \n',rhos);
fprintf(fid3,'Stratified correlation coefficient per diagonal:\n');
for i=1:numRF-1
    fprintf(fid3,'%4.3e \n',rhoktemp(i));
end
fprintf(fid3,'Weight of stratified correlation coefficient per diagonal:\n');
for i=1:numRF-1
    fprintf(fid3,'%4.3e \n',weightktemp(i));
end
fclose(fid3);

end

