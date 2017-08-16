function [stats] = simsam_stats(samval,samact,simsam)

%% reshape data 
rep = size(samval,2);
nvox = size(samval,1);
FSval = simsam.FSBBvec;
FSact = simsam.FSBBact;

actvox = sum(FSact);
FSvalmat = repmat(FSval',1,rep);
FSactmat = repmat(FSact',1,rep);

% Effect Size
% ------------

% mean 
samsig = samval.*double(samact);
samsigm = sum(abs(samsig))./sum(samact);
mEffsam = nanmean(samsigm);
mEffsamCI95 = prctile(samsigm(~isnan(samsigm)),[2.5 97.5]);
stats.meaneff = samsigm;
stats.meaneffM = mEffsam;
stats.meaneffCI95 = mEffsamCI95;

% max 
sammax = max(samval);
stats.maxeff = sammax;
stats.maxeffM = mean(sammax);
stats.maxeffCI95 = prctile(sammax,[2.5 97.5]);

% percentage of significant voxels
% -------------------------------
actper = (sum(samact)./nvox).*100;
mActper  = mean(actper);
mActperCI = prctile(actper,[2.5 97.5]);
stats.mactper = actper;
stats.mactperM = mActper;
stats.mactperCI = mActperCI;

% Differences Subsequent replications
% -----------------------------
s1 = samval(:,1:end-1);
s2 = samval(:,2:end);
d = s1-s2; d = abs(d);
dm = mean(d);
stats.diffrep  = dm;
stats.diffrepM = mean(dm);
stats.diffrepCI95 = prctile(dm,[2.5 97.5]);

% Dice Coefficient
% ---------------
act1 = samact(:,1:end-1);
act2 = samact(:,2:end);
sameact = act1==1 & act2==1;
sumsa = sum(sameact);
sumtot= sum([act1;act2]);
DCvec = (2.*sumsa)./sumtot;
% set NaN to zero - no overlap. 
DCvec(isnan(DCvec)) = 0;
stats.dicerep = DCvec;
stats.dicerepM = mean(DCvec); 
stats.dicerepCI95 = prctile(DCvec,[2.5 97.5]); 

% Spatial Correlation
% -----------------------
for j = 1:(rep-1)
    if sum(sameact(:,j))>1;
        sigcorr= corr(s1(sameact(:,j),j),s2(sameact(:,j),j),'type','Spearman');
        spatcor(j) = sigcorr;
    else
        spatcor(j) = 0;
    end
end
stats.correp  = spatcor;
stats.correpM = mean(spatcor);
stats.correpCI95 = prctile(spatcor,[2.5 97.5]);

%% sensitivity / specificity 
% --------------------------
true_pos = (samact==1 & FSactmat==1);
TPvec = sum(true_pos);
TP = mean(TPvec);

false_pos = (samact==1 & FSactmat==0);
FPvec = sum(false_pos);
FP = mean(FPvec);

true_neg = (samact==0 & FSactmat==0);
TNvec = sum(true_neg);
TN = mean(TNvec);

false_neg = (samact==0 & FSactmat==1);
FNvec = sum(false_neg);
FN = mean(FNvec);

tot = TP + FP + TN + FN;

stats.tp = TP;
stats.fp = FP;
stats.tn = TN;
stats.fn = FN;

sens = TP/(TP + FN);
sensvec = TPvec./(TPvec+FNvec);
sensCI95 = prctile(sensvec,[2.5 97.5]);
spec = TN/(TN + FP);

%FWE
FWE = sum(FPvec>0)/rep;

sensANY = sum(TPvec>0)/rep;

%FDR
FDROVER = FP/(FP+TP);

%Accuracy
ACC = (TP + TN)/(nvox);

%false positive rate.
FPR = 1-spec;

LRP = sens/(1-spec);
LRN = (1-sens)/spec;
ER = (1-sens)/(1-spec);

%precision
prec = TP/(TP+FP);

%F1
F1 = 2*TP/(2*TP+FP+FN);

% add to structure
stats.spec = spec;
stats.sens = sens;
stats.sensCI95 = sensCI95;
stats.ACC  = ACC;
stats.FPR  = FPR;
stats.FDR  = FDROVER;
stats.sensANY = sensANY;
stats.FWE  = FWE;
stats.F1   = F1;
stats.prec = prec;
stats.LRN = LRN;
stats.LRP = LRP;
stats.ER  = ER;

end
     

