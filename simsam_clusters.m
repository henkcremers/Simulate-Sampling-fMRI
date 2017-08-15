function [coord,clustermap] = simsam_clusters(map)
% [coord,clustermap] = simsam_clusters(map)
nr = size(map,1);
nc = size(map,2);

rcoor = repmat([1:nr]',1,nc);
ccoor = repmat([1:nc],nr,1);

% create vector 
vec = reshape(map,1,nr*nc);
coord(1,:) = reshape(ccoor,1,nr*nc);
coord(2,:) = reshape(rcoor,1,nr*nc);
% only for the inbrainvoxel
coord = coord(:,vec==1);

clusters = spm_clusters([coord; ones(1,size(coord,2))],18);
clustervec = zeros(1,nr*nc);
clustervec(vec) = clusters;
clustermap = reshape(clustervec,nr,nc);

end

