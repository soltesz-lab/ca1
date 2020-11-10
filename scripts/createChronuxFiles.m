%% simDat = createChronuxFils(resultsDir)
% For a simulation (specified by resultsDir), creates a structure
% (simDat)that contains all the spike times and LFP in a Chronux parasble
% manner
% 
% inputs:
%   resultsDir: string that specifies the directory in which of the files
%   from the run are saved
% outputs:
%   simDat: structure with LFP and cell spike data
%       fields:
%           LFP - (# of time steps)x2 array where the first column is the t
%           (ms) and the second column is the value of the local field
%           potential
%           params - parameters used by filtering functions in Chronux
%           library. There are hard-coded parameters here that may need to
%           be fixed. 
%           cellTypeNames - cell array of names of cell types used by
%           Bezaire & Raikov's code
%           cellTypeInds - (# of cell types)x3 array. 1st column is the
%           index of the cell type used by the simulation source code. 2nd
%           column is the first index (0-indexing) for the cells of this
%           cell type. 3rd column is the last index (0-indexing) for the 
%           cells of this cell type
%           rasterPlotCell - (# of cell types)x1 cell array. Each cell
%           contains an array of spike times (1st column). The index of the
%           cell that spiked is found in the 2nd column
%           cellSpikes - Chronux requires each cell to have its own
%           structure of spike times. The format to acheive this is as
%           follows:
%               cellSpike: (# of cell types)x1 struct array, each with
%               field "dat"
%                   dat: 1x(# of cells) structure
%                       dat(i).times: times of spikes for ith cell of jth
%                       cell type
%           
%           
% Mark Plitt - July 2016

function simDat = createChronuxFiles(resultsDir)
%% set things up
% hardcode base directory for results
dataDir = fullfile('/home/igr/src/model/ca1/results/',resultsDir);

cd(dataDir);

% Make sure LFP was saved in this results directory. Heuristic to ensure
% run finished on cluster
if ~exist('lfp.dat','file');
    disp('No LFP data found. Check to make sure run finished!');
    ME = MException('MyComponent:noSuchFile', 'LFP data is not found');
else
    disp('Found the LFP data! :D');
    load('lfp.dat');
    simDat.LFP = lfp(2:end,:);
end


%% set up chronux parameter structure
params.Fs = length(simDat.LFP(2:end,1))/simDat.LFP(end,1); % sampling frequency of LFP(kHz)  
params.tapers = [3 5]; % [TW (time-bandwidth product), K (# of tapers)] 
params.fpass= [0 .3]; % look at frequencies between 0 and 300 Hz
% params.err = [2 .05]; % jacknife (LOOCV) error calculation at 95% CI

simDat.params = params;


%% create structures for cell spike times for use in chronux
% get info on cell types
cellTypes = importdata('celltype.dat','\t'); 
cellTypeInds = cellTypes.data; %columns: type index, range start, range end
cellTypeNames = cellTypes.textdata(2:end,1); 
simDat.cellTypeNames = cellTypeNames; % save in structure
simDat.cellTypeInds = cellTypeInds;

% concatenate spike raster files
files = dir('spikeraster.dat');
files = {files(:).name};

% dump spike rasters into one large array
spikes = importdata(files{1});
for i = 2:length(files)
    tmp = importdata(files{i});
    spikes = [spikes; tmp];
end

% sort raster by cell index
[~,I] = sort(spikes(:,2));
spikes_s = spikes(I,:);

disp('spike rasters concatenated and sorted by cell');

% unique cell numbers
cellNums = unique(spikes_s(:,2));

% put spike times for each cell type into their own arrays and save in a cell
% array
spikeDatCell = cell(size(cellTypeInds,1),1);
for type = 1:size(cellTypeInds,1)
    mask = (spikes_s(:,2)>=cellTypeInds(type,2) & spikes_s(:,2)<=cellTypeInds(type,3));
    spikeDatCell{type} = spikes_s(mask,:);
%     figure; scatter(spikeDatCell{type}(:,1),spikeDatCell{type}(:,2),'.k'); xlim([0 600]);
end
simDat.rasterPlotCell = spikeDatCell;
disp('spike raster for plotting code created');


% create structure for each cell - required by Chronux
% cellType x 1 structure
%   each entry is a N cells x 1 structure with fieldname times
disp('creating structure of spike times for each cell');
simDat.cellSpikes = struct([]);
for type = 1:size(cellTypeInds,1)
    cellN_vec = cellTypeInds(type,2):cellTypeInds(type,3);
    cellDat = struct([]);
%     parfor cellI = 1:length(cellN_vec)
    for cellI = 1:length(cellN_vec)
        cellN = cellN_vec(cellI)
        cellDat(cellI).times = spikes_s(spikes_s(:,2)==cellN,1);
%         cellDat(cellI).cellN = cellN;    
    end
    
    simDat.cellSpikes(type,1).dat = cellDat;
    disp([cellTypeNames{type} ' is done...']);
end

%save('simulationData.mat','simDat');

end
