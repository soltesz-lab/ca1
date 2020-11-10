%% Plot subset of intracellular traces for CA1 ripple simulations
% This script was used for debugging when I was having numerical stability
% problems

% gather filenames
cd 'C:\cygwin64\home\mplitt\ca1\results\ca1_ripplestim_1x_44_dtChange';

cellTypes = {'axoaxoniccell','bistratifiedcell','cckcell','ivycell',...
    'ngfcell','olmcell','pyramidalcell','pvbasketcell','scacell'};

names = cell(length(cellTypes),1);
for c = 1:length(cellTypes)
    cellStruc = dir(strcat('trace_',cellTypes{c},'*.dat'));
    cellStruc = {cellStruc(:).name};
    names{c} = cellStruc;
end


for n = 1:2:length(names{1})
%     figure;
%     for c = 1:length(cellTypes)
%         subplot(length(cellTypes),1,c);
%         cellDat = importdata(names{c}{n},'\t');
%         plot(cellDat.data(2500:5000,1),cellDat.data(2500:5000,2));
%         title(cellTypes{c});    
%     end
    
    for c=[7 8]
        cellDat = importdata(names{c}{n},'\t');
        figure; plot(cellDat.data(:,1),cellDat.data(:,2));
        title(cellTypes{c});
    end
    
    pause
end
        