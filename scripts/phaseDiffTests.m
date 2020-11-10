% compute whether there is significant differences in phase preferences of
% different cell types

function [mu1,mu2,muDiff,pval_ww,pval_perm] = phaseDiffTests(results,simDat,cellType1,cellType2)
nPerms = 1e4;

[phaseCount,~,~] = phasePref(results, simDat,0);
% phases = results.LFPbands.phaseLFP;
if isnumeric(cellType1) && isnumeric(cellType2)
    pd1 = phaseCount{cellType1}; pd2 = phaseCount{cellType2};
else
    pd1 = phaseCount{strcmp(simDat.cellTypeNames,cellType1)};
    pd2 = phaseCount{strcmp(simDat.cellTypeNames,cellType2)};
end
% watson-williams test
[pval_ww, table]= circ_wwtest(pd1,pd2);

% calculate mean angle of circular distribution
mu1 = circ_mean(pd1'); mu2 = circ_mean(pd2');
% difference between means
muDiff = circ_dist(mu1,mu2);

id1 = zeros(length(pd1),1); id2 = ones(length(pd2),1);
idC = [id1; id2]; pdC = [pd1'; pd2'];
N = length(idC);
nullDist = NaN(nPerms,1);
% permutation test
for p = 1:nPerms

% shuffle assignment of cell type
    order = randperm(N);
    tmpID = idC(order);
    
    tmpPd1 = pdC(tmpID==0); tmpPd2 = pdC(tmpID==1);
    
    % compute difference in circular mean phase
    tmp_mu1= circ_mean(tmpPd1); tmp_mu2 = circ_mean(tmpPd2);
    nullDist(p) = circ_dist(tmp_mu1,tmp_mu2);
end


% get pvalue for permutation test
pval_perm = sum(muDiff<nullDist)/nPerms;