% Preliminary script used to visualize spatial organization of firing.
% This could be used to look at a ripple oscillation as a "traveling wave"

load simulationData.mat
% results directory
posFile = 'position.dat';
posStruct = importdata(posFile);

% plot x,y position and each time frame will be a scatter plot
time = simDat.LFP(:,1);
if ~exist('2DPlotMat.mat','file')
    % get max an min values for x and y
    range = [min(posStruct.data(:,2:3)); max(posStruct.data(:,2:3))];
    % bin space for plotting
    xBins = linspace(range(1,1),range(2,1),21); yBins = linspace(range(1,2),range(2,2),21);
    
    % allocate for speed
    plotDat = zeros(length(xBins)-1,length(yBins)-1,length(time));
    
    cellIDs = cell(length(xBins)-1,length(yBins)-1); % cell array of cellIDs for each bin
    % Put spike data into bins
    for x = 1:length(xBins)-1
        if x ~= length(xBins)-1
            xRows = (posStruct.data(:,2)>=xBins(x) & posStruct.data(:,2) < xBins(x+1));
        else
            xRows = (posStruct.data(:,2)>=xBins(x) & posStruct.data(:,2) <= xBins(x+1));
        end
        
        for y = 1:length(yBins)-1
            if y ~= length(yBins)-1
                yRows = (posStruct.data(:,3)>=yBins(y) & posStruct.data(:,3) < yBins(y+1));
            else
                yRows = (posStruct.data(:,3)>=yBins(y) & posStruct.data(:,3) <= yBins(y+1));
            end
            rows = xRows & yRows;
            
            cellIDs{x,y} = posStruct.data(rows,1);
        end
    end
    
    allSpikes = cell2mat(simDat.rasterPlotCell);
    
    for x = 1:size(cellIDs,1);
        for y = 1:size(cellIDs,2);
            tmpSpikes = [];
            for id = 1:length(cellIDs{x,y})
                tmpSpikes = [tmpSpikes; allSpikes(allSpikes(:,2)==cellIDs{x,y}(id),1)];
            end
            
            for t = 1:length(time)
                plotDat(x,y,t) = sum(tmpSpikes==time(t));
            end
        end
    end
    
    % save data
    save 2DPlotMat.mat plotDat
else
    
    load 2DPlotMat.mat
end

% plot binned data as heat map movie 
subplot(5,1,1);
plot(simDat.LFP(:,1),simDat.LFP(:,2)); hold on;
plot([3000 3000], [-15, 15]);
ax1 = gca;
ax1.NextPlot = 'replaceChildren';

subplot(5,1,2:5)
imagesc(plotDat(:,:,3000),[0 20]); 
ax2 = gca;
ax2.NextPlot = 'replaceChildren';

loops = 1001; size(plotDat,3);
F(loops) = struct('cdata',[],'colormap',[]);
for j = 1:loops
    
    subplot(5,1,1);
    plot(simDat.LFP(:,1),simDat.LFP(:,2),'b'); hold on;
%     hard coded start and end times
    plot([simDat.LFP(j+3000-1,1) simDat.LFP(j+3000-1,1)], [-15, 15],'r'); hold off;
    
    subplot(5,1,2:5)
    imagesc(plotDat(:,:,j+3000-1),[0 20]);
    drawnow;
    F(j) = getframe(gcf);
end
save -v7.3 spikeRasterMovie.mat F

