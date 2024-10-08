function xMat=zScoreX(xMat, replaceNans)

if nargin<2
    replaceNans=false;
end




if nanstd(xMat(:,1))==0
    xMat=xMat(:,2:end);
    
    
    xMat=(xMat-repmat(nanmean(xMat, 1), size(xMat,1), 1))./...
        repmat(nanstd(xMat), size(xMat,1), 1);
    
    xMat=[ones(size(xMat,1),1) xMat];
else
    
    xMat=(xMat-repmat(nanmean(xMat, 1), size(xMat,1), 1))./...
        repmat(nanstd(xMat), size(xMat,1), 1);
    
        
end

if replaceNans
   xMat(~isfinite(xMat))=0;
end
