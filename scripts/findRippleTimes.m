% define ripple event using some "objective" criteria and calculate the
% coherence and phase-preference within this event only

function [ripStart,envStart, ripStop, envStop] = findRippleTimes(results,simDat)

% beginning and end of ripple is onset of first CA3 spike and offset of
% last CA3 spike
ripSpike = sort(simDat.rasterPlotCell{12}(:,1));
ripStart = floor(ripSpike(1)*10); ripStop = floor(ripSpike(end)*10);

% beginning of ripple is when envelope of filtered LFP rises above 0+1 std 
envAll = results.LFPbands.envLFP; envStd = std(envAll(:,4)); % standardDeviation
envStart = find(envAll(:,4)>envStd,1); envStop = find(envAll(envStart:end,4)<envStd,1)+envStart;
% cutoff for end of ripple is when envelope of filtered LFP falls below 0+1 std 

