function [rhos, rhoktemp, weightktemp]=CompSCC(map1,map2,dmax,smoothh)
if any(any(map1>1))
    Xold=map1./(max(max(map1))); % normalize maps
else
    Xold=map1;
end
if any(any(map2>1))
    Yold=map2./(max(max(map2)));
else
    Yold=map2;
end

if smoothh
    X=smoothmap(Xold,smoothh);
    Y=smoothmap(Yold,smoothh);
    
else
    X=Xold;
    Y=Yold;
end

NB=length(Xold); % Number of RFs

% need to optimize h
rhoktemp=zeros(dmax,1);
weightktemp=zeros(dmax,1);
r1k=zeros(dmax,1);
r2k=zeros(dmax,1);
r2kV=zeros(dmax,1);
for d=1:dmax % for every distance
    X1=[];Y1=[];
    for i=1:(NB-d)
        X1(i)=X(i,i+d);
        Y1(i)=Y(i,i+d);
    end
    if numel(unique(X1)) && numel(unique(Y1)) >1
        
        r1k(d)=mean(X1.*Y1)-mean(X1).*mean(Y1);
        r2k(d)=sqrt( var(X1,1).*var(Y1,1));
        tempCorr=corrcoef(X1,Y1);
        rhoktemp(d)=tempCorr(2);
        
        [~,rankX]=ismember(X1,sort(X1,'descend'));
        [~,rankY]=ismember(Y1,sort(Y1,'descend'));
        sortX=rankX./max(rankX);
        sortY=rankY./max(rankY);
        
        nXY=length(X1);
        r2kV(d)=sqrt( var(sortX,1).*var(sortY,1));
        weightktemp(d)=nXY.*r2k(d);
    else
        rhoktemp(d)=NaN;
        weightktemp(d)=NaN;
    end
end

rhok=rhoktemp(~(isnan(rhoktemp) | isnan(weightktemp)));
weightk=weightktemp(~(isnan(rhoktemp) | isnan(weightktemp)));

rhos=sum(rhok.*weightk)./sum(weightk);

%% functions
    function newX=smoothmap(X,h)
        U=length(X);
        newX=zeros(U);
        for l=1:U
            for r=1:U
                m=max([1, l-h]);
                n=max([1, r-h]);
                o=min([U, l+h]);
                p=min([U, r+h]);
                subMat=X(m:o,n:p);
                newX(l,r)=sum(sum(subMat))./((1+2.*h).^2);
            end
        end
    end
end