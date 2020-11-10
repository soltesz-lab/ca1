%% [phaseCount, mvl, pVals] = phasePref(results,simDat,varargin)
% calculate preferred phase of multi-unit spiking activity within ripple
% frequency band (130-200). I have yet to generalize this to do any phase
% 
% inputs:
%   results: described in spectralAnalyses_subsample
%   simDat: described in createChronuxFiles
%   [runPerm]: if set to 1, run permutations test for significance of phase
%   preference, if 0, don't run permuations test
% 
% outputs:
%   phaseCounts: histogram of spikes in each phase bin
%   mvl: mean vector length of population spiking
%   pVals: permutation test achieved p-value for phase preference



function [phaseCount,mvl,pVals] = phasePref(results, simDat,varargin)
%%
ripI = 4; % index in results to get ripple data
% start and stop of ripple - code needs to be improved
[ripStart,envStart, ripStop, envStop] = findRippleTimes(results,simDat);
eventTimes = envStart:envStop;

if ~isempty(varargin)
    runPerm = varargin{1};
else
    runPerm=1;
end


phaseDat = results.LFPbands.phaseLFP;
instPhase = phaseDat(eventTimes,ripI); % time series of instantaneous phase
phaseCount = cell(length(simDat.cellTypeNames),1);
pVals = NaN(length(simDat.cellTypeNames),1);
mvl = NaN(length(simDat.cellTypeNames),1);

for c = 1:length(simDat.cellTypeNames)
  
     if ~ismember(simDat.cellTypeNames{c},{'eccell','ca3cell','ca3ripcell'})
         tmpCount = [];
         mua = results.MUA(c).spikeCount(eventTimes);
         
         spikeInds = mua>0; % find indices of spikes and make a mask
         
         phase = instPhase(spikeInds); 
         mua = mua(spikeInds);
         
         % get counts of phase values
         for s = 1:length(phase);
             tmpCount = [tmpCount repmat(phase(s),1,mua(s))];
         end
         if length(tmpCount)==0
             tmpCount = zeros(1,length(phase));
         end
         
         phaseCount{c} = tmpCount;
         
         if runPerm
            [mvl(c),p,mvlDist] = phasePermTest(tmpCount,instPhase);
            pVals(c) = p;
         else
             [histS,edges] =  histcounts(tmpCount,[-180:20:180]*pi/180,'normalization','probability');
             mvl = norm(histS,2);
         end
     end
end
end

function [mvl,p,nullDist] = phasePermTest(phaseDist,instPhase)
    nPerms = 1e4; % 10000 permutations
    totalSpikes = length(phaseDist);
    nullSpikeInds = randi(length(instPhase),[nPerms totalSpikes]); % random data indexing
    
    % get histogram of spike phases and calculate mean vector length
    [histS,edges] =  histcounts(phaseDist,[-180:20:180]*pi/180,'normalization','probability');
    mvl = norm(histS,2);
    
    
    nullDist = zeros(nPerms,1);
    for i = 1:nPerms % scramble data and repeat
        [histS,edges] =  histcounts(instPhase(nullSpikeInds(i,:)),[-180:20:180]*pi/180,'normalization','probability');
        nullDist(i) = norm(histS,2);
    end
    
    % calculate p value
    p=sum(nullDist>mvl)/nPerms;
end
