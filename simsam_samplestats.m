function [samplestats,sampledata] = simsam_samplestats(simsam,varargin)

% basic info
% ---------------------------
np = length(simsam.behav);    % full sample size
nr = size(simsam.inbrain,1);
nc = size(simsam.inbrain,2);
nvox = length(simsam.coord);  % number of inbrain voxels

FSval = simsam.FSBBvec;       % full sample values
FSact = simsam.FSBBact;       % full sample logical vector

% defaults.
rep = 2000;                            % number of samples in total
nfigdat = 4;                           % number of samples to save sampledata for figure
nsamp = 10:5:150;                      % vector with the sample sizes
thresholds = {0.05 0.01 'FDR' 'BONF'}; % thresholds numeric (0-1),'FDR' or 'BONF')

% input.
for i = 1:length(varargin)
    arg = varargin{i};
    if ischar(arg)
        switch lower(arg)
            case 'rep'
                rep = varargin{i+1};
            case 'nsamp'
                nsamp = varargin{i+1};
            case 'thresholds'
                thresholds = varargin{i+1};
            case 'nfigdat'
                nfigdat = varargin{i+1};
        end
    end
end
tic
%% start
% ------
for s = 1:length(nsamp)
    
    disp(['Sample Size: ' num2str(nsamp(s))])
    for k = 1:rep
        
        fullsample = 1:np;
        sample = randsample(fullsample,nsamp(s));
        psample = simsam.behav(sample,1);
        bsample = simsam.data(sample,:);
        
        % BB effects
        [dat,pval] = corr(psample,bsample);
        
        % main effects
        % [h,p,CI,STATS] = ttest(bsample);
        % dat = STATS.tstat;
        % pval = p;
        
        repdat(:,k) = dat;
        reppval(:,k) = pval;
        
        % fix zero p values
        pval(pval<1e-20)=1e-20;
        
        % Get the FDR value
        thrFDR = FDR(pval,0.05);
        if ~isempty(thrFDR); thrFDRloc = pval<thrFDR; else thrFDRloc = zeros(1,length(pval)); thrFDR = NaN; end
        repfdr.val(1,k) = thrFDR;
        repfdr.loc(:,k) = thrFDRloc;
        
        % save the first samples for figures
        if k<=nfigdat
            sampledata.nsamp(s,1).n = nsamp(s);
            sampledata.nsamp(s,1).repd(:,1:k) = repdat(:,1:k);
            sampledata.nsamp(s,1).repp(:,1:k) = reppval(:,1:k);
            sampledata.nsamp(s,1).fdrloc(:,k) = thrFDRloc;
            sampledata.nsamp(s,1).fdrval(1,k) = thrFDR;
            sampledata.nsamp(s,1).sample(:,k) = sample;
        end
    end
    
    %% Set up the sample data per treshold
    for t = 1:length(thresholds)
        if isnumeric(thresholds{t})
            disp(['  Threshold: ' num2str(thresholds{t})])
            thr =  thresholds{t};
            samact = reppval<thr;
            samval = repdat;
        elseif strmatch(thresholds{t},'FDR')
            disp(['  Threshold: ' thresholds{t}])
            % filter out extremely low FDR values - ? 
            mfdr = nanmedian(repfdr.val);
            filter = repfdr.val>(mfdr/100) | isnan(repfdr.val);
            %loc = repfdr.loc;
            samval = repdat(:,filter);
            samact = repfdr.loc(:,filter);
        elseif strmatch(thresholds{t},'BONF')
            disp(['  Threshold: ' thresholds{t}])
            samact = reppval<(0.05/nvox);
            samval = repdat;
        else
            error('unknown threshold')
        end
        
        % calculate the stats & save in samplestats structure
        samplestats.thr(t).level{1} = thresholds{t};
        samplestats.thr(t).n = nsamp;
        [stats] = simsam_stats(samval,samact,simsam);
        
        samplestats.thr(t).spec(s,1) = stats.spec;
        samplestats.thr(t).sens(s,1) = stats.sens;
        samplestats.thr(t).sensCI95(s,1:2) = stats.sensCI95;
        
        samplestats.thr(t).ACC(s,1) = stats.ACC;
        samplestats.thr(t).FPR(s,1) = stats.FPR;
        samplestats.thr(t).FDR(s,1) = stats.FDR;
        samplestats.thr(t).sensANY(s,1) = stats.sensANY;
        samplestats.thr(t).FWE(s,1) = stats.FWE;
        samplestats.thr(t).F1(s,1) = stats.F1;
        samplestats.thr(t).prec(s,1) = stats.prec;
        samplestats.thr(t).LRP(s,1) = stats.LRP;
        samplestats.thr(t).LRN(s,1) = stats.LRN;
        samplestats.thr(t).ER(s,1) = stats.ER;
        
        % core measures
        samplestats.thr(t).meaneff{s,1} = stats.meaneff;
        samplestats.thr(t).meaneffM(s,1) = stats.meaneffM;
        samplestats.thr(t).meaneffCI95(s,1:2) = stats.meaneffCI95;
        
        samplestats.thr(t).maxeff{s,1} = stats.maxeff;
        samplestats.thr(t).maxeffM(s,1) = stats.maxeffM;
        samplestats.thr(t).maxeffCI95(s,1:2) = stats.maxeffCI95;
        
        samplestats.thr(t).actper{s,1} = stats.mactper;
        samplestats.thr(t).actperM(s,1) = stats.mactperM;
        samplestats.thr(t).actperCI95(s,1:2) = stats.mactperCI;
        
        samplestats.thr(t).diffrep{s,1} = stats.diffrep;
        samplestats.thr(t).diffrepM(s,1) = stats.diffrepM;
        samplestats.thr(t).diffrepCI95(s,1:2) = stats.diffrepCI95;
        
        samplestats.thr(t).corrrep{s,1} = stats.correp;
        samplestats.thr(t).corrrepM(s,1) = stats.correpM;
        samplestats.thr(t).corrrepCI95(s,1:2) = stats.correpCI95;
        
        samplestats.thr(t).DC{s,1} = stats.dicerep;
        samplestats.thr(t).DCM(s,1) = stats.dicerepM;
        samplestats.thr(t).DCCI95(s,1:2) = stats.dicerepCI95;
        
        
        
    end
end
toc
end
