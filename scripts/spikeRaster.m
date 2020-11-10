%% spikeRaster(simDat)
% simply plots the spike raster for the entire population broken up by cell
% type. As you can see I used the brute force method (ctrl+c & ctrl+v), and
% hardcoded everything. Marianne's code is much better, but did not work on
% my computer. 
% Mark Plitt - July 2016
function spikeRaster(simDat)

spikeDatCell = simDat.rasterPlotCell;
figure;
% play with proportions
subplot(800,1,1:13);
ymin = min(spikeDatCell{1}(:,2)); ymax = max(spikeDatCell{1}(:,2));
scatter(spikeDatCell{1}(:,1),spikeDatCell{1}(:,2),'.k'); xlim([0 600]); ylim([ymin ymax]);
ylabel('AA'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,15:35);
ymin = min(spikeDatCell{2}(:,2)); ymax = max(spikeDatCell{2}(:,2));
scatter(spikeDatCell{2}(:,1),spikeDatCell{2}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('BS');set(gca,'YTick',[],'XTick',[]);

subplot(800,1,37:71);
ymin = min(spikeDatCell{3}(:,2)); ymax = max(spikeDatCell{3}(:,2));
scatter(spikeDatCell{3}(:,1),spikeDatCell{3}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('CCK');set(gca,'YTick',[],'XTick',[]);

subplot(800,1,73:159);
ymin = min(spikeDatCell{4}(:,2)); ymax = max(spikeDatCell{4}(:,2));
scatter(spikeDatCell{4}(:,1),spikeDatCell{4}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('Ivy'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,161:195);
ymin = min(spikeDatCell{5}(:,2)); ymax = max(spikeDatCell{5}(:,2));
scatter(spikeDatCell{5}(:,1),spikeDatCell{5}(:,2),'.k'); xlim([0 600]); ylim([ymin ymax]);
ylabel('NGF'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,197:212);
scatter(spikeDatCell{6}(:,1),spikeDatCell{6}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('OLM'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,214:380);
ymin = min(spikeDatCell{7}(:,2)); ymax = max(spikeDatCell{7}(:,2));
scatter(spikeDatCell{7}(:,1),spikeDatCell{7}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('Pyr'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,382:420);
ymin = min(spikeDatCell{8}(:,2)); ymax = max(spikeDatCell{8}(:,2));
scatter(spikeDatCell{8}(:,1),spikeDatCell{8}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('PV'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,422:450);
ymin = min(spikeDatCell{9}(:,2)); ymax = max(spikeDatCell{9}(:,2));
scatter(spikeDatCell{9}(:,1),spikeDatCell{9}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('SCA'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,452:550);
ymin = min(spikeDatCell{10}(:,2)); ymax = max(spikeDatCell{10}(:,2));
scatter(spikeDatCell{10}(:,1),spikeDatCell{10}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('EC'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,552:670);
ymin = min(spikeDatCell{11}(:,2)); ymax = max(spikeDatCell{11}(:,2));
scatter(spikeDatCell{11}(:,1),spikeDatCell{11}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('CA3'); set(gca,'YTick',[],'XTick',[]);

subplot(800,1,672:800);
ymin = min(spikeDatCell{12}(:,2)); ymax = max(spikeDatCell{12}(:,2));
scatter(spikeDatCell{12}(:,1),spikeDatCell{12}(:,2),'.k'); xlim([0 600]);ylim([ymin ymax]);
ylabel('CA3Rip'); set(gca,'YTick',[]);


