function [setup] = simsam_setup_data(setupinfo)
% ========================================================
% set-up the brain-behaviour correlations and main effects. 
% ========================================================

%% set up the BB effects
[activation rrange] = create_clusters(setupinfo.BBclusters,setupinfo.BBrange,setupinfo.BBsmooth,[110 110]);

% plot & crop the data
[brainmap,contour,underlay] = simsam_plotslice(activation,setupinfo.BBrange(1),'or',setupinfo.orientation,'slicenum',setupinfo.slice,'name','Brain-Behavior Effects');

% add the info
setup.BBmap = brainmap;
setup.BBcontour = contour;
setup.inbrain = underlay.inbrain;
setup.BBvec = simsam_reshape(setup.BBmap,setup.inbrain);
setup.underlay = underlay;
setup.BBrange = rrange;

% cluster info BB
[coord,clustermap] = simsam_clusters(abs(brainmap)>0);
setup.BBclustermap = clustermap;

% coordinates of all in-brain voxels
[coord,clustermap] = simsam_clusters(setup.inbrain);
setup.coord = coord;

%% set up the main effects
[activation range] = create_clusters(setupinfo.MEclusters,setupinfo.MErange,setupinfo.MEsmooth,[110 110]);
[brainmap,contour,underlay] = simsam_plotslice(activation,setupinfo.MErange(1),setup.underlay,'name','MainEffects');
setup.MEmap = brainmap;
setup.MEcontour = contour;
[coord,clustermap] = simsam_clusters(abs(brainmap)>0);
setup.MEclustermap = clustermap;
setup.MErange = range;
return

function [activation rrange] = create_clusters(clusters,range,smooth,dim);
nclus = size(clusters,1);
activation = zeros(dim(2),dim(1));

for cl = 1:nclus;
    x = clusters(cl,1);
    y = clusters(cl,2);
    r = clusters(cl,3);
    [rr cc] = meshgrid(1:dim);
    blob = sqrt((rr-x).^2+(cc-y).^2)<=r;
    clusterloc(cl).dat = blob;
    loc = blob==1;
    blob = blob .*clusters(cl,4);
    activation(loc) = blob(loc);
end

% smooth
if smooth>0
    activation = simsam_smooth(activation,smooth);
end

% scale
activation = ((activation./max(max(activation)).*range(2)));

% Cut off values below minimum
cutoffloc = abs(activation)<range(1) & abs(activation)>0;
activation(cutoffloc)=0;
actabs = abs(activation);
rrange(1) = min(actabs(actabs>0));
rrange(2) = max(actabs(actabs>0));
return