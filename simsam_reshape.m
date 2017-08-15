function rdat = simsam_reshape(dat,inbrain)

if isvector(dat)
direction = 'vec2map';
else 
direction = 'map2vec';
end

nr = size(inbrain,1);
nc = size(inbrain,2);
ibvec = reshape(inbrain,1,nr*nc);

switch direction
    case 'vec2map'
        rdat = zeros(1,nr*nc);
        rdat(ibvec)=dat;
        rdat = reshape(rdat,nr,nc);
    case 'map2vec'
        rdat = reshape(dat,1,nr*nc);
        rdat = rdat(ibvec);
end

return
