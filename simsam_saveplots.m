function [] = simsam_saveplots(simsam,sampledata,setup,varargin)

rep = 4;
samplesize = [30 150];
thresholds = {0.05 0.01 'FDR' 'BONF'};

% input.
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'rep'
                rep = varargin{i+1};
            case 'samplesize'
                samplesize = varargin{i+1};
            case 'thresholds'
                treshold = varargin{i+1};
        end
    end
end

% create a new dir

mkdir('figures')
cd('figures')

%% plot the figures for the full sample

% plot the slice
maxvox = simsam.coord(:,simsam.FSBBmaxloc);
[data,contour] = simsam_plotslice(simsam.FSBBmap,simsam.FSBBmin,'underlay',setup.underlay,'plotcir',[maxvox(1) maxvox(2)],'saveplot','fullsample');
nvox = length(simsam.coord);

% % plot the scatter
maxdat = [simsam.behav simsam.data(:,simsam.FSBBmaxloc)];
plotname = 'fullsample_corrmaxplot';
simsam_plotcor(maxdat,'saveplot',plotname);

%% save some subsamples..

ns = length(sampledata.nsamp);
for j=1:ns; for i=1:length(samplesize); if samplesize(i) == sampledata.nsamp(j).n; nvec(i) = j;end; end; end

for k = 1:rep;
    for nloc = nvec;
        n = sampledata.nsamp(nloc).n;
        rmap = sampledata.nsamp(nloc).repd(:,k);
        [maxval maxloc] = max(abs(rmap));
        rmap = simsam_reshape(rmap,simsam.inbrain);
        pmap = sampledata.nsamp(nloc).repp(:,k);
        pmap = simsam_reshape(pmap,simsam.inbrain);
        
        %%
        close all;
        for t = 1:length(thresholds);
            
            if isnumeric(thresholds{t})
                thr =  thresholds{t};
                loc =  pmap<thr;
            elseif strmatch(thresholds{t},'FDR')
                loc = sampledata.nsamp(nloc).fdrloc(:,rep);
                loc = simsam_reshape(loc,simsam.inbrain);
            elseif strmatch(thresholds{t},'BONF')
                loc = pmap<(0.05/nvox);
            else
                error('unknown threshold')
            end
            if isnumeric(thr); thr = num2str(thr); else thr = thresholds{t}; end
            % treshold the map
            rmap(~loc) = 0;
            
            imname = ['n_' num2str(n) '_rep' num2str(k) '_thr' thr];
            if sum(sum(loc))>0;
                simsam_plotslice(rmap,0.05,'contour',contour,'plotcir',simsam.coord(:,maxloc)','saveplot',imname);
            else
                simsam_plotslice(rmap,0.05,'contour',contour,'saveplot',imname);
            end
            close all
            
        end
        
        %% save the scatter plot
        sample = sampledata.nsamp(nloc).sample(:,k);
        sampledat = [simsam.behav(sample) simsam.data(sample,maxloc)];
        plotname = ['n_' num2str(n) '_rep' num2str(k) '_scatter'];
        simsam_plotcor(sampledat,'saveplot',plotname);
    end
end


return

