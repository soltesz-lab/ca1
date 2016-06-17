function Plot_Spike_Pies(~,~,handles)
global RunArray

if exist('data\RastArray.mat','file')==0

    idx=searchRuns('ExecutionDate','',0,'~');
    
    for m=1:length(idx)
        ind = idx(m);

        handles.curses=[];
        handles.curses.ind = ind;
        guidata(handles.btn_generate,handles);
        
        spikeraster(handles.btn_generate,guidata(handles.btn_generate))
        handles=guidata(handles.btn_generate);

        getcelltypes(handles.btn_generate,guidata(handles.btn_generate))
        handles=guidata(handles.btn_generate);

        RastArray(m).RunName = RunArray(ind).RunName;
        RastArray(m).spikerast = handles.curses.spikerast;
        RastArray(m).cells = handles.curses.cells;

        for x=1:length(RastArray(m).cells)
            RastArray(m).CellbySpike(x).Name = RastArray(m).cells(x).name(1:3);
            RastArray(m).CellbySpike(x).NumSpikes = length(find(RastArray(m).spikerast(:,2)>=RastArray(m).cells(x).range_st & RastArray(m).spikerast(:,2)<=RastArray(m).cells(x).range_en));
            RastArray(m).CellbySpike(x).SpikeRate = RastArray(m).CellbySpike(x).NumSpikes/RastArray(m).cells(x).numcells;
        end

        RastArray(m).DegreeStim = RunArray(ind).DegreeStim;
        RastArray(m).SimDuration = RunArray(ind).SimDuration;
        RastArray(m).NumData = RunArray(ind).NumData;
        RastArray(m).ConnData = RunArray(ind).ConnData;
        RastArray(m).SynData = RunArray(ind).SynData;
        RastArray(m).ind = ind;
    end
    save('data\RastArray.mat','RastArray','-v7.3');
else
    load('data\RastArray.mat');
end

colorvec={'b','c','g','y','m','r'};

idxb=find([RastArray(:).ConnData]>306);
useme='SpikeRate'; % NumSpikes SpikeRate

for h=307:315
    bh(h-306).f = figure('Color','w','Name',['ConnData ' num2str(h)]);
    bh(h-306).i = 1;
    bh(h-306).runs = find([RastArray(:).ConnData]==h);
    bh(h-306).stims = [RastArray(bh(h-306).runs).DegreeStim];
    [~, tmpi] = sort(bh(h-306).stims);
    bh(h-306).runs = bh(h-306).runs(tmpi);
    for z=1:length(bh(h-306).runs)
        m = bh(h-306).runs(z);
        figure(bh(RastArray(m).ConnData-306).f);
        subplot(3,5,bh(RastArray(m).ConnData-306).i);
        bh(RastArray(m).ConnData-306).i = bh(RastArray(m).ConnData-306).i + 1;
        try
            pie([RastArray(m).CellbySpike([RastArray(m).CellbySpike(:).(useme)]>0).(useme)],{RastArray(m).CellbySpike([RastArray(m).CellbySpike(:).(useme)]>0).Name})
        catch
            m
        end
        title([RastArray(m).RunName ': '  num2str(RastArray(m).DegreeStim)]) % num2str(RastArray(m).ConnData) ', '
    end
end

% 
% for z=1:length(idxb)
%     m=idxb(z);
% %     if bh(RastArray(m).ConnData-306).i>15
% %         try
% %             figure(bh(RastArray(m).ConnData-306+20).f);
% %         catch
% %             bh(RastArray(m).ConnData-306+20).f = figure('Color','w','Name',['Next ConnData ' num2str(RastArray(m).ConnData)]);
% %         end
% %     else
%         figure(bh(RastArray(m).ConnData-306).f);
% %     end
%     subplot(3,5,bh(RastArray(m).ConnData-306).i);
%     bh(RastArray(m).ConnData-306).i = bh(RastArray(m).ConnData-306).i + 1;
%     try
%         pie([RastArray(m).CellbySpike([RastArray(m).CellbySpike(:).(useme)]>0).(useme)],{RastArray(m).CellbySpike([RastArray(m).CellbySpike(:).(useme)]>0).Name})
%     catch
%         m
%     end
%     title([RastArray(m).RunName ': '  num2str(RastArray(m).DegreeStim)]) % num2str(RastArray(m).ConnData) ', '
% end

% for x=1:length(RastArray(1).cells)
%     figure('Color','w');
%     for z=1:length(idxb)
%         m=idxb(z);
%         if RastArray(m).ConnData<304
%             mycol = 'k';
%         else
%             mycol = colorvec{floor((RastArray(m).ConnData-304)/2)+1};
%         end
%         if RastArray(m).DegreeStim<2
%             plot(RastArray(m).DegreeStim, RastArray(m).CellbySpike(x).NumSpikes/RastArray(m).SimDuration,'MarkerSize',15,'Marker','.','Color',mycol);
%         end
%         hold on
%     end
%     xlabel('DegreeStim')
%     ylabel('Spikes/time')
%     title(RastArray(1).CellbySpike(x).Name)
% end
%     
