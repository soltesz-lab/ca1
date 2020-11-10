%% results = spectralAnalyses_subsample(simDat)
% perform some basic signal processing on the outputs of the simulation
% using Chronux. For single-cell analyses, just takes a random subsample of
% cells
% 
% inputs: 
%   simDat: described in createChronuxFiles.m. A structure with
%   Chronux-friendly formatting of simulation results
% 
% outputs:
%   results: structure containing some basic signal processing results
%       fields:
%           spectrum: spectrum of LFP, fields {S (power), f (freq)} 
%           spectrogram: spectrogram of LFP, fiels {S,f,t (time)}
%           LFPbands: 1x(# of frequency bands) with fields
%                   bandOrder - name of frequency band
%                   filtLFP - filtered LFP
%                   phaseLFP - instantaneous phase of filtered LFP
%                   envLFP - envelope of filtered LFP
% Mark Plitt - July 2016

function results = spectralAnalyses_subsample(simDat)
%% plot initial results before analyzing spike times

%plot LFP
figure;
plot(simDat.LFP(:,1),simDat.LFP(:,2));
xlabel('time (ms)'); ylabel('V (mV)');
title('LFP');

% plot spike raster
spikeRaster(simDat) % 

% plot spectrum
[S,f] = mtspectrumc(simDat.LFP(:,2),simDat.params); % calculate spectrum for continuous timeseries
figure;
plot_vector(S,f*1000,'n'); ylabel('Power (dB)'); % chronux function for plotting spectrum
spectrum.S = S; spectrum.f = f;
results.spectrum = spectrum;

% plot spectrogram
[S,t,f] = mtspecgramc(simDat.LFP(:,2),[100 1],simDat.params); % calculates spectrogram
figure;
plot_matrix(S,t,f*1000,'n'); % chronux function for spectrogram plots
spectrogram.S = S; spectrogram.t = t; spectrogram.f = f;
results.spectrogram = spectrogram;

%% do some basic filtering, coherence, and phase preference calculations

% filter LFP for different frequency bands
bandOrder = {'theta', 'epsilon','gamma', 'ripple'};
fpass = [5 10; 25 90; 90 130; 130 200]; % hard-coded frequency bands

Fs = f; %10000; %sampling rate (Hz)
% Fs_new = 1000; % downsample to this rate
n = 4; %order of filter

% downsample 
% LFP_ds = resample(simDat.LFP(:,2),1,10);
% LFP_t_ds = resample(simDat.LFP(:,1),1,10);
% results.LFP_ds = LFP_ds; results.LFP_t_ds = LFP_t_ds;

filtLFP = zeros(length(simDat.LFP(:,2)),length(bandOrder));
phaseLFP = zeros(length(simDat.LFP(:,2)),length(bandOrder));
envLFP = zeros(length(simDat.LFP(:,2)),length(bandOrder));

for band = 1:length(bandOrder) % for each frequency band
   
%     if band ~=1
%         Hd = designfilt('bandpassfir','FilterOrder', 20, 'CutoffFrequency1',fpass(band,1),...
%             'CutoffFrequency2',fpass(band,2), 'SampleRate',1000);
%         filtDat = filtfilt(Hd,LFP_ds);
%     else
%         
%         [b,a] = butter(n,fpass(band,:).*2/Fs_new);
%         filtDat = filtfilt(b,a,LFP_ds);
%     end

    [b,a] = butter(n,fpass(band,:).*2/Fs); % use basic butterworth filter 
    filtDat = filtfilt(b,a,simDat.LFP(:,2));
    filtDat_a = hilbert(filtDat); % create analytical signal

%     unrolled_phase = unwrap(angle(filtDat_a));
    inst_phase = angle(filtDat_a); % instantaneous phase of signal
    env = abs(filtDat_a); % calculate envelope of signal

    filtLFP(:,band) = filtDat;
    phaseLFP(:,band) = inst_phase;
    envLFP(:,band) = env;
end
results.LFPbands = struct('bandOrder',bandOrder,'filtLFP',filtLFP,...
    'phaseLFP',phaseLFP,'envLFP',envLFP);


fprintf('plotted basic frequency analysis of LFP. \n');
pause;
close all;



%% coherence and phase-locking

% for each cell type
for ct = 1:length(simDat.cellTypeNames)
    
    % single-cell analysis
    
    % subsample nCells random cells of each type
    nCells = 100;
    order = randperm(length(simDat.cellSpikes(ct).dat));
    order = order(1:nCells);
    
    if ~ismember(simDat.cellTypeNames{ct},{'eccell','ca3cell','ca3ripcell'})
        
        LFP = repmat(simDat.LFP(:,2),[1,nCells]);
        
        [C,phi,S12,S1,S2,f]=coherencycpt(LFP,... % calculates coherence between continuous and point process
            simDat.cellSpikes(ct).dat(order),simDat.params);
        
        [C_g,phi_g,S12,S1,S2,t_g,f_g]=cohgramcpt(LFP,... % same as above but time-binned 
            simDat.cellSpikes(ct).dat(order),[100 .1],simDat.params);
        
        coherence(ct).C = C; % coherence
        coherence(ct).phi = phi; % phase
        coherence(ct).f = f; % frequency
        
        coherogram(ct).C = C_g;
        coherogram(ct).phi = phi_g;
        coherogram(ct).t = t_g; % time
        coherogram(ct).f =f_g;
        
        fprintf('%s coherence estimate is done \n',simDat.cellTypeNames{ct});
    else
        fprintf('skipping coherence estimate for %s cell type \n',simDat.cellTypeNames{ct});
        coherence(ct).C = 'skipping cell type';
        coherence(ct).phi = 'skipping cell type';
        coherogram(ct).C = 'skipping cell type';
        coherogram(ct).phi = 'skipping cell type';
        
    end
    
    
    % calculate multi-unit coherence and phase-locking for each cell type
    
    %     spike count
    if ~ismember(simDat.cellTypeNames{ct},{'eccell','ca3cell','ca3ripcell'})
        sc = zeros(size(simDat.LFP,1),1);
        t = 0:.1:599.9;
        for t_i = 1:length(t)-1;
            sc(t_i) = sum((simDat.rasterPlotCell{ct}(:,1)>=t(t_i) & simDat.rasterPlotCell{ct}(:,1)<t(t_i+1)));
        end
        MUA(ct).spikeCount = sc;
        
        % coherence between continuous signal and probability distribution
        [C,phi,S12,S1,S2,f]=coherencycpb(simDat.LFP(:,2),sc,simDat.params); 
        
        [C_g,phi_g,S12,S1,S2,t_g,f_g]=cohgramcpb(simDat.LFP(:,2),sc,[100 .1],simDat.params);
        
        

        MUA(ct).coherence = struct('C',C,'phi',phi,'f',f);

        MUA(ct).coherogram = struct('C',C_g,'phi',phi_g,'f',f_g,'t',t_g);
        

         fprintf('%s multi-unit coherence estimate is done \n',simDat.cellTypeNames{ct});
    else
        MUA(ct).coherence = 'skipping cell type';
        MUA(ct).coherogram = 'skipping cell type';
    end
    
    
end
results.coherence = coherence;
results.coherogram = coherogram;
results.MUA = MUA;

save -v7.3 results.mat results
