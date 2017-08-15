function [simsam] = simsam_generate_data(setup,varargin)

% Defaults
% ------------
mbehav = 30;  % mean behavioural values
n = 10000;    % full sample size
stdbehav = 3; % standard deviation behavioural data
stdbrain = 3; % standard deviation brain data
FWHM = 6;     % FWHM

% Input.
% -------------
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'mbehav'
                mbehav = varargin{i+1};
            case 'n'
                n = varargin{i+1};
            case 'stdbehav'
                stdbehav = varargin{i+1};
            case 'stdbrain'
                stdbrain = varargin{i+1};
            case 'FWHM'
                FWHM = varargin{i+1};
        end
    end
end

simsam.FWHM = FWHM;
simsam.inbrain = setup.inbrain;
nr = size(simsam.inbrain,1);
nc = size(simsam.inbrain,2);
nn = sum(sum(simsam.inbrain));

%% =========================================================================
% Set-up the covariance
% =========================================================================

% set-up some variables.
bbVox = simsam_reshape(setup.BBmap,setup.inbrain);
simsam.coord = setup.coord;
distvox = dist(simsam.coord); simsam.distvox = distvox;
nV = length(bbVox);
% mu = [mbehav zeros(1,nV)];
meVox = simsam_reshape(setup.MEmap,setup.inbrain);

% main effects
mu = [mbehav meVox.*stdbrain];
s  = [stdbehav^2 stdbrain^2*ones(1,nV)];
Sigma = zeros(nV+1,nV+1);

% Smooth
sig = FWHM./(2*(sqrt(2*log(2))));
smoothmap = gaussmf(distvox,[sig 0]);

% Difference in Correlation among voxels related to behavioural variable:
dMat = abs(repmat(bbVox',1,nn)) - abs(repmat(bbVox,nn,1));
pMat = bbVox'*bbVox;
pMat = sign(pMat).*((abs(pMat)).^0.5);

% set-up the cormat:
cormat = Sigma;
cormat(1,:) = [1 bbVox];
cormat(:,1) = [1 bbVox];

% Covariance - compromise between distance (smooth) and correlation with B
sMat = smoothmap;
% copy the smoothmap, weigh with the discrepenacy..
cMat = sMat.*(1-dMat);      % - 20170812 needed???
% + "hard coded" overrule of smoothness
locP = abs(pMat)>sMat & abs(pMat)>0.6;
cMat(locP) = pMat(locP);
cormat(2:end,2:end) = cMat; %update the correlation matrix

% Set-up the covariance
Sigma(1,1) = stdbehav.^2;
Sigma(2:end,:) = stdbehav*stdbrain;
Sigma(:,2:end) = stdbehav*stdbrain;
Sigma(2:end,2:end) = stdbrain.^2;
Sigma = cormat.*Sigma;
Sigma(logical(eye(nV+1,nV+1))) = s;
Sigma = (Sigma + Sigma.') / 2;

% Create SPD covariance matrix.
disp('set-up covariance..')
tic; SigmaSPD = nearestSPD(Sigma); toc

% True Effects
simsam.Sigma = SigmaSPD;
simsam.TBBvec = simsam.Sigma(1,2:end)./(stdbrain*stdbehav);
simsam.TBBmap = simsam_reshape(simsam.TBBvec,simsam.inbrain);

%%  create the Full sample data
% ======================
disp('create the data..')
tic; dat = mvnrnd(mu, SigmaSPD, n); toc

simsam.behav = dat(:,1);
simsam.data = dat(:,2:end);

[c p] = corrcoef(dat);

simsam.FSBBvec = c(1,2:end);
simsam.FSBBp   = p(1,2:end);
simsam.FSMEvec = mean(simsam.data)./std(simsam.data);
clear dat;

% Maps of the full sample data
simsam.FSBBmap = simsam_reshape(simsam.FSBBvec,simsam.inbrain);
simsam.FSMEmap = simsam_reshape(simsam.FSMEvec,simsam.inbrain);

%% add some info to the structure
% -------------------------------

FSBBabs = abs(simsam.FSBBvec);
FSBBact = FSBBabs>=0.05;
minact = min(FSBBabs(FSBBact));
maxnull = max(FSBBabs(~FSBBact));

[coord,clustermap] = simsam_clusters(abs(simsam.FSBBmap)>=minact);
simsam.BBclustermap = clustermap;
simsam.BBclustervec = simsam_reshape(clustermap,simsam.inbrain);

actsum = sum(FSBBabs(FSBBact));
actperc = (sum(FSBBact)./length(FSBBact))*100;
actmean = mean(FSBBabs(FSBBact));
[maxact maxloc] = max(FSBBabs);

simsam.FSBBeffsum = actsum;
simsam.FSBBmean = actmean;
simsam.FSBBabs = FSBBabs;
simsam.FSBBact = FSBBact;
simsam.FSBBactperc = actperc;
simsam.FSBBmax = maxact;
simsam.FSBBmaxloc = maxloc;
simsam.FSBBmin = minact;
simsam.FSBBmaxnull = maxnull;

return

