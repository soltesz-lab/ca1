%%  plotSpectralAnalyses(simDat,results)
% plot results from spectralAnalyses_subsample
% summary results for LFP as well as cell type specific plots are created
% inputs: simDat described in createChronuxFiles.m
%         results described in spectralAnalyses_subsample
% 
% Mark Plitt July 2016

function plotSpectralAnalyses(simDat,results)

% LFP
plot(simDat.LFP(:,1),simDat.LFP(:,2)); xlabel('time (ms)'); ylabel('Potential (mV)');
title('LFP');

% spike-raster
spikeRaster(simDat);

% spectrum
figure;
plot_vector(results.spectrum.S,results.spectrum.f*1000,'n'); ylabel('Power (dB)');


% spectrogram
figure;
plot_matrix(results.spectrogram.S,results.spectrogram.t,results.spectrogram.f*1000,'n'); % chronux function for spectrogram plots


% find the start and end of the ripple 
[ripStart,envStart, ripStop, envStop] = findRippleTimes(results,simDat); % this code is a little buggy
phaseCount = phasePref(results,simDat,0); % determine phase of spiking for each cell type
for ct = 1:length(simDat.cellTypeNames) % for each cell type
    
    if ~ismember(simDat.cellTypeNames{ct},{'eccell','ca3cell','ca3ripcell'}) % skip stimulus cells
        % plot LFP with superimposed population activity
        figure;subplot(4,2,1:2); 
        yyaxis left 
        l = plot(simDat.LFP(envStart-500:envStop+500,1),simDat.LFP(envStart-500:envStop+500,2));
        yyaxis right
        m = bar(simDat.LFP(envStart-500:envStop+500,1),results.MUA(ct).spikeCount(envStart-500:envStop+500),'r');
        
        title(simDat.cellTypeNames{ct});
        % coherence
        subplot(4,2,3);
        plot_vector(results.MUA(ct).coherence.C, results.MUA(ct).coherence.f,'n'); 
        ylabel('Coherence'); xlabel('f');

        % coherogram
        subplot(4,2,4);
        plot_matrix(results.MUA(ct).coherogram.C, results.MUA(ct).coherogram.t,...
            results.MUA(ct).coherogram.f,'n'); 
        ylabel('f'); xlabel('time'); title('');
        
        % ripple-filtered data
        subplot(4,2,5:6)
        filtLFP = results.LFPbands.filtLFP; 
        phaseLFP = results.LFPbands.phaseLFP;
        plot(simDat.LFP(envStart-500:envStop+500,1),filtLFP(envStart-500:envStop+500,4)); hold on;
%         plot(simDat.LFP(envStart-500:envStop+500,1),phaseLFP(envStart-500:envStop+500,4));
        
        % phase-pref
        subplot(4,2,7)
        yyaxis left
        histogram(phaseCount{ct}); 
        yyaxis right 
        plot(-pi:.01:pi,cos(-pi:.01:pi)); xlim([-pi pi]);
        
        % rose plot of phase-pref
        subplot(4,2,8);
        rose(phaseCount{ct})
        
    end
end