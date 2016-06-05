function Theta_Power_PS_Old(hObject,handles)
global RunArray 

idx1=searchRuns('ExecutionDate','',0,'~');
idx2=searchRuns('Scale',1,1,'=');

idx3=intersect(idx1,idx2);
idx4=searchRuns('SimDuration',250,1,'>');

idx5=intersect(idx3,idx4);
idx6=searchRuns('Stimulation','spontaneous',0,'=');

idx=intersect(idx5,idx6);


rez=1;
croptime = 100;
endtime = 1000;

if exist('myRuns.mat','file')
    load('myRuns.mat','myRuns');
    rs = length(myRuns)+1;
    chkflag=1;
else
    rs=1;
    chkflag=0;
end

for r=1:length(idx)
    if chkflag==1 && ~isempty(strmatch(RunArray(idx(r)).RunName,{myRuns(:).Name},'exact'))
        continue
    end
    myRuns(rs).Name =  RunArray(idx(r)).RunName;
    myRuns(rs).Ind =  idx(r);
    myRuns(rs).SimDuration =  RunArray(idx(r)).SimDuration;
    myRuns(rs).NumData =  RunArray(idx(r)).NumData;
    myRuns(rs).ConnData =  RunArray(idx(r)).ConnData;
    myRuns(rs).SynData =  RunArray(idx(r)).SynData;
    myRuns(rs).DegreeStim =  RunArray(idx(r)).DegreeStim;

    handles.curses = [];
    handles.curses.ind = idx(r);
    guidata(handles.btn_generate,handles)
    spikeraster(handles.btn_generate,guidata(handles.btn_generate))
    handles=guidata(handles.btn_generate);
    getcelltypes(handles.btn_generate,guidata(handles.btn_generate))
    handles=guidata(handles.btn_generate);
    numcons(handles.btn_generate,guidata(handles.btn_generate));
    handles=guidata(handles.btn_generate);
    handles.curses.spikerast = addtype2raster(handles.curses.cells,handles.curses.spikerast,3);
    guidata(handles.btn_generate, handles)
    
    timeidx = find(handles.curses.spikerast(:,1)>croptime & handles.curses.spikerast(:,1)<endtime);
    handles.curses.spikerast = handles.curses.spikerast(timeidx,:);
    
    flagr = 1;
    for c=1:length(handles.curses.cells)
        spidx = find(handles.curses.spikerast(:,2)>=handles.curses.cells(c).range_st & handles.curses.spikerast(:,2)<=handles.curses.cells(c).range_en);
        idt = find(handles.curses.spikerast(spidx,1)>croptime & handles.curses.spikerast(spidx,1)<endtime);

        [f, fft_results]=myfft(handles.curses.spikerast(spidx(idt),1),RunArray(idx(r)).SimDuration,rez);
        theta_range=find(f(:)>=4 & f(:)<=12);
        [~, theta_idx] = max(fft_results(theta_range));
        gamma_range=find(f(:)>=25 & f(:)<=100);
        [~, gamma_idx] = max(fft_results(gamma_range));
        
        handles.curses.cells(c).theta.freq = f(theta_range(theta_idx));
        handles.curses.cells(c).theta.power = fft_results(theta_range(theta_idx));
        handles.curses.cells(c).theta.norm = fft_results(theta_range(theta_idx))/max(fft_results);
        
        handles.curses.cells(c).gamma.freq = f(gamma_range(gamma_idx));
        handles.curses.cells(c).gamma.power = fft_results(gamma_range(gamma_idx));
        handles.curses.cells(c).gamma.norm = fft_results(gamma_range(gamma_idx))/max(fft_results);

        y=histc(handles.curses.spikerast(spidx(idt),2),handles.curses.cells(c).range_st:handles.curses.cells(c).range_en);
        z=hist(y,[0:max(y)]);
        handles.curses.cells(c).spikehist.numcells = z;
        handles.curses.cells(c).spikehist.numspikes = [0:max(y)];
        handles.curses.cells(c).poprate = mean(y/((endtime-croptime)/1000));
        handles.curses.cells(c).subrate = mean(y(y~=0)/((endtime-croptime)/1000));
        handles.curses.cells(c).totspikes = sum(y);
        
        if strcmp(handles.curses.cells(c).name,'pyramidalcell')==1
            flagr = c;
        end
    end
    
    % Theta Period
if rs==40
    rs
end
    thetaper = 1000/handles.curses.cells(flagr).theta.freq;

    pyrangle = 0;
    for c=1:length(handles.curses.cells)
        spidx = find(handles.curses.spikerast(:,2)>=handles.curses.cells(c).range_st & handles.curses.spikerast(:,2)<=handles.curses.cells(c).range_en);
        idt = find(handles.curses.spikerast(spidx,1)>croptime & handles.curses.spikerast(spidx,1)<endtime);
        spiketimes = handles.curses.spikerast(spidx(idt),1);

        n=length(spiketimes);
        modspiketimes = mod(spiketimes, thetaper);

        xbar = 1/n*sum(sin(modspiketimes*pi/(thetaper/2)));
        ybar = 1/n*sum(cos(modspiketimes*pi/(thetaper/2)));

        magnitude(c)=sqrt(xbar^2+ybar^2);
        if xbar>0
            angle = acos(ybar/magnitude(c));
        else
            angle = 2*pi - acos(ybar/magnitude(c));
        end

        rdir(c) = angle; % * pi/180;
        if strcmp(handles.curses.cells(c).name,'pyramidalcell')==1
            pyrangle=angle;
        end
    end
    
    refphase = 20;
    refangle=refphase/180*pi;
    pyrangle_shift=pyrangle-refangle;
    
    for c=1:length(handles.curses.cells)
        handles.curses.cells(c).theta.phase = mod((rdir(c)-pyrangle_shift)*180/pi+360,360);
        handles.curses.cells(c).theta.mag = magnitude(c);
        handles.curses.cells(c).theta.pyrangle = pyrangle*180/pi;
    end
    
     
    % Gamma Period
    
    thetaper = 1000/handles.curses.cells(flagr).gamma.freq;

    pyrangle = 0;
    for c=1:length(handles.curses.cells)
        spidx = find(handles.curses.spikerast(:,2)>=handles.curses.cells(c).range_st & handles.curses.spikerast(:,2)<=handles.curses.cells(c).range_en);
        idt = find(handles.curses.spikerast(spidx,1)>croptime);
        spiketimes = handles.curses.spikerast(spidx(idt),1);

        n=length(spiketimes);
        modspiketimes = mod(spiketimes, thetaper);

        xbar = 1/n*sum(sin(modspiketimes*pi/(thetaper/2)));
        ybar = 1/n*sum(cos(modspiketimes*pi/(thetaper/2)));

        magnitude(c)=sqrt(xbar^2+ybar^2);
        if xbar>0
            angle = acos(ybar/magnitude(c));
        else
            angle = 2*pi - acos(ybar/magnitude(c));
        end

        rdir(c) = angle; % * pi/180;
        if strcmp(handles.curses.cells(c).name,'pyramidalcell')==1
            pyrangle=angle;
        end
    end
    
    refphase = 0;
    refangle=refphase/180*pi;
    pyrangle_shift=pyrangle-refangle;
    
    for c=1:length(handles.curses.cells)
        handles.curses.cells(c).gamma.phase = (rdir(c)-pyrangle_shift)*180/pi;
        handles.curses.cells(c).gamma.mag = magnitude(c);
        handles.curses.cells(c).gamma.pyrangle = pyrangle*180/pi;
    end
       
    myRuns(rs).Theta = handles.curses.cells(flagr).theta;
    myRuns(rs).Gamma = handles.curses.cells(flagr).gamma;
    myRuns(rs).Cells = handles.curses.cells;
    
    myRuns(rs).Spikes.Real = 0;
    myRuns(rs).Spikes.Pyr = 0;
    myRuns(rs).Spikes.Inrn = 0;
    
    for c=1:length(myRuns(rs).Cells)
        if strcmp(myRuns(rs).Cells(c).techname(1:2),'pp')==0
            myRuns(rs).Spikes.Real = myRuns(rs).Spikes.Real + myRuns(rs).Cells(c).totspikes;
            if strcmp(myRuns(rs).Cells(c).name,'pyramidalcell')==1
                myRuns(rs).Spikes.Pyr = myRuns(rs).Spikes.Pyr + myRuns(rs).Cells(c).totspikes;
            else
                myRuns(rs).Spikes.Inrn = myRuns(rs).Spikes.Inrn + myRuns(rs).Cells(c).totspikes;
            end
        end
    end
    
    myRuns(rs).Spikes.PyrFrac = myRuns(rs).Spikes.Pyr/myRuns(rs).Spikes.Real;
    
    rs = rs + 1;
    if mod(r,10)==0
        save('myRuns.mat','myRuns');
        disp(['Just finished run ' num2str(r) ' out of ' num2str(length(idx)) ' runs.'])
    end
end

save('myRuns.mat','myRuns');

% idx2use = 1:length(myRuns);

idx2useA = find([myRuns(:).NumData]==106 | [myRuns(:).NumData]==101);
idx2useB = find([myRuns(:).SimDuration]>=600);

idx2use=intersect(idx2useA,idx2useB);

pidx=[1     2     4     6    38    43    39    40    41    42    45 51   148   149   156];

bidx = setxor(idx2use,pidx);
tidx = pidx;

subsint(1).idx=bidx;
subsint(1).col='k';
subsint(1).desc='All other runs';

subsint(2).idx=[1     2     4     6];
subsint(2).col='m';
subsint(2).desc='MayFlower runs';

subsint(3).idx=[38    43    39    40    41    42    45 51];
subsint(3).col='g';
subsint(3).desc='Basic runs';
% 
subsint(4).idx=[148   149   156];
subsint(4).col='c';
subsint(4).desc='Cut runs';


% tmpidx = find(arrayfun(@(x) x.Theta.norm, myRuns(idx2use))==1);
% focusidx = idx2use(tmpidx);
% 
% tmpidx = find(arrayfun(@(x) x.Theta.norm, myRuns(idx2use))<1);
% bidx = idx2use(tmpidx);
% 
% subsint(1).idx=bidx;
% subsint(1).col='k';
% subsint(1).desc='Cut <1 theta norm';
% 
% subsint(2).idx=focusidx;
% subsint(2).col='g';
% subsint(2).desc='Cut max theta norm';
% 
% 
% idx2useA = find([myRuns(:).NumData]==101);
% idx2useB = find([myRuns(:).SimDuration]>350);
% 
% idx2use=intersect(idx2useA,idx2useB);
% 
% tmpidx = find(arrayfun(@(x) x.Theta.norm, myRuns(idx2use))==1);
% focusidx = idx2use(tmpidx);
% 
% tmpidx = find(arrayfun(@(x) x.Theta.norm, myRuns(idx2use))<1);
% bidx = idx2use(tmpidx);
% 
% subsint(3).idx=bidx;
% subsint(3).col='b';
% subsint(3).desc='Pool <1 theta norm';
% 
% subsint(4).idx=focusidx;
% subsint(4).col='m';
% subsint(4).desc='pool max theta norm';
if ispc
    sl='\';
else
    sl='/';
end

repos=RunArray(1).ModelDirectory;
dc=dir([repos sl 'datasets' sl 'conndata_*.dat']);
for d=1:length(dc)
    myconns(d).num = str2num(dc(d).name(10:end-4));
    fid = fopen([repos sl 'datasets' sl dc(d).name],'r');                
    numlines = fscanf(fid,'%d\n',1) ;
    filedata = textscan(fid,'%s %s %f %f %f\n') ;
    fclose(fid);
    for r=1:numlines
        myconns(d).data.(filedata{1}{r}).(filedata{2}{r}).weight = filedata{3}(r);
        myconns(d).data.(filedata{1}{r}).(filedata{2}{r}).numcons = filedata{4}(r);
        myconns(d).data.(filedata{1}{r}).(filedata{2}{r}).syns = filedata{5}(r);
    end
    try
        myconns(d).ec2pyr=myconns(d).data.eccell.pyramidalcell.weight;
    catch
        myconns(d).ec2pyr=0;
    end
    try
        myconns(d).ca32pyr=myconns(d).data.ca3cell.pyramidalcell.weight;
    catch
        myconns(d).ca32pyr=0;
    end
end

for m=1:length(myRuns)
    c=find([myconns.num]==myRuns(m).ConnData);
    myRuns(m).conns = myconns(c);
end


figure('Color','w','Name','Pyramidal Fraction of Spikes')
for s=1:length(subsint)
    plot([myRuns(subsint(s).idx).DegreeStim],arrayfun(@(x) x.Spikes.PyrFrac, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
% for t=1:length(subsint(2).idx)
%     text(myRuns(subsint(2).idx(t)).DegreeStim,myRuns(subsint(2).idx(t)).Spikes.PyrFrac,myRuns(subsint(2).idx(t)).Name)
% end
legend({subsint.desc})
xlabel('Stimulation (Hz)')
ylabel('Pyramidal Spike Fraction')



colorvec=[.0 .0 .6;
          .0 .75 .65;
          1 .75 .0;
          1 .5 .3;
          1 .0 .0;
          .6 .4 .1;
          .5 .0 .6;
          .6 .6 .6;
          1 .1 1;
          1 0 0;
          0 0 1;];
      
mh = figure('Color','w','Name',['Mayflower Theta Freq']);
[x,y] = pol2cart(0,1);
h_fake1=compass(x,y);
hold on

bh = figure('Color','w','Name',['Basic Theta Freq']);
[x,y] = pol2cart(0,1);
h_fake2=compass(x,y);
hold on

ch = figure('Color','w','Name',['Cut Theta Freq']);
[x,y] = pol2cart(0,1);
h_fake3=compass(x,y);
hold on

for z=1:length(myRuns(1).Cells)
    figure(mh)
    for s=1:4
        try
        [x,y] = pol2cart(arrayfun(@(x) x.Cells(z).theta.phase, myRuns(subsint(2).idx(s))),arrayfun(@(x) x.Cells(z).theta.mag, myRuns(subsint(2).idx(s))))
        k=compass(x,y);
        z = strmatch(myRuns(1).Cells(z).name,{'pyramidalcell','pvbasketcell','cckcell','scacell','axoaxoniccell','bistratifiedcell','olmcell','ivycell','ngfcell','supercell','deepcell'});
        if isempty(z)
        set(k,'Color','k','LineWidth',4);
        else
        set(k,'Color',colorvec(z,:),'LineWidth',4);
        end
        hold on
        end
    end
    
    figure(bh)
    for s=5:12
        try
        [x,y] = pol2cart(arrayfun(@(x) x.Cells(z).theta.phase/180*pi, myRuns(subsint(2).idx(s))),arrayfun(@(x) x.Cells(z).theta.mag, myRuns(subsint(2).idx(s))))
        k=compass(x,y);
        z = strmatch(myRuns(1).Cells(z).name,{'pyramidalcell','pvbasketcell','cckcell','scacell','axoaxoniccell','bistratifiedcell','olmcell','ivycell','ngfcell','supercell','deepcell'});
        if isempty(z)
        set(k,'Color','k','LineWidth',4);
        else
        set(k,'Color',colorvec(z,:),'LineWidth',4);
        end
        hold on
        end
    end
    
    figure(ch)
    for s=13:15
        try
        [x,y] = pol2cart(arrayfun(@(x) x.Cells(z).theta.phase/180*pi, myRuns(subsint(2).idx(s))),arrayfun(@(x) x.Cells(z).theta.mag, myRuns(subsint(2).idx(s))))
        k=compass(x,y);
        z = strmatch(myRuns(1).Cells(z).name,{'pyramidalcell','pvbasketcell','cckcell','scacell','axoaxoniccell','bistratifiedcell','olmcell','ivycell','ngfcell','supercell','deepcell'});
        if isempty(z)
        set(k,'Color','k','LineWidth',4);
        else
        set(k,'Color',colorvec(z,:),'LineWidth',4);
        end
        hold on
        end
    end
end

set(h_fake1,'Visible','off')
set(h_fake2,'Visible','off')
set(h_fake3,'Visible','off')


figure('Color','w','Name','Pyramidal Contribution')
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Spikes.Real, myRuns(subsint(s).idx)),arrayfun(@(x) x.Spikes.PyrFrac, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
for t=1:length(subsint(2).idx)
    text(myRuns(subsint(2).idx(t)).DegreeStim,myRuns(subsint(2).idx(t)).Spikes.PyrFrac,myRuns(subsint(2).idx(t)).Name)
end
legend({subsint.desc})
xlabel('# Spikes')
ylabel('Pyramidal Spike Fraction')


figure('Color','w','Name','# Spikes v. Stimulation')
for s=1:length(subsint)
    plot(arrayfun(@(x) x.DegreeStim, myRuns(subsint(s).idx)),arrayfun(@(x) x.Spikes.Real, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
for t=1:length(subsint(2).idx)
    text(myRuns(subsint(2).idx(t)).DegreeStim,myRuns(subsint(2).idx(t)).Spikes.PyrFrac,myRuns(subsint(2).idx(t)).Name)
end
legend({subsint.desc})
xlabel('Stimulation')
ylabel('# Spikes')


figure('Color','w','Name','stim v firing rate')
for c=1:length(subsint)
    plot(arrayfun(@(x) x.DegreeStim, myRuns([subsint.idx])),arrayfun(@(x) x.Cells(c).subrate, myRuns([subsint.idx])),'Color',colorvec(c,:),'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
xlabel('Stimulation')
ylabel('rate')




figure('Color','w','Name','Spike Frac v Theta Freq')
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Spikes.PyrFrac, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
for t=1:length(subsint(2).idx)
    text(myRuns(subsint(2).idx(t)).Spikes.PyrFrac,myRuns(subsint(2).idx(t)).Theta.freq,myRuns(subsint(2).idx(t)).Name)
end
legend({subsint.desc})
xlabel('Pyramidal Spike Fraction')
ylabel('Theta Freq. (Hz)')

figure('Color','w','Name','Spike Frac v Gamma Freq')
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Spikes.PyrFrac, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
for t=1:length(subsint(2).idx)
    text(myRuns(subsint(2).idx(t)).Spikes.PyrFrac,myRuns(subsint(2).idx(t)).Gamma.freq,myRuns(subsint(2).idx(t)).Name)
end
legend({subsint.desc})
xlabel('Pyramidal Spike Fraction')
ylabel('Gamma Freq. (Hz)')



% return

figure('Color','w','Name','Stimulation v CA3->Pyr')
%subplot(1,2,1)
for s=1:length(subsint)
    plot([myRuns(subsint(s).idx).DegreeStim],arrayfun(@(x) x.conns.ca32pyr, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Stimulation (Hz)')
ylabel('CA3->Pyr')
title('Conn')

figure('Color','w','Name','Stimulation v Theta')
subplot(1,2,1)
for s=1:length(subsint)
    plot([myRuns(subsint(s).idx).DegreeStim],arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Stimulation (Hz)')
ylabel('Theta Freq. (Hz)')
title('Theta Frequency')

subplot(1,2,2)
for s=1:length(subsint)
    plot([myRuns(subsint(s).idx).DegreeStim],arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Stimulation (Hz)')
ylabel('Theta Power(|Y|)')
title('Theta Power')

figure('Color','w','Name','Stimulation v Gamma')
subplot(1,2,1)
for s=1:length(subsint)
    plot([myRuns(subsint(s).idx).DegreeStim],arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Stimulation (Hz)')
ylabel('Gamma Freq. (Hz)')
title('Gamma Frequency')

subplot(1,2,2)
for s=1:length(subsint)
    plot([myRuns(subsint(s).idx).DegreeStim],arrayfun(@(x) x.Gamma.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Stimulation (Hz)')
ylabel('Gamma Power(|Y|)')
title('Gamma Power')


figure('Color','w','Name','Freq v. Power')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Theta Freq. (Hz)')
ylabel('Theta Power(|Y|)')
title('Theta Freq. v Power')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Gamma Freq. (Hz)')
ylabel('Gamma Power(|Y|)')
title('Gamma Freq. v Power')


figure('Color','w','Name','Freq v. Norm Power')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Theta Freq. (Hz)')
ylabel('Theta Norm Power(|Y|)')
title('Theta Freq. v Norm Power')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Gamma Freq. (Hz)')
ylabel('Gamma Norm Power(|Y|)')
title('Gamma Freq. v Norm Power')


figure('Color','w','Name','Freq v. Other Power')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Theta Freq. (Hz)')
ylabel('Gamma Norm Power(|Y|)')
title('Theta Freq. v Gamma Power')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Gamma Freq. (Hz)')
ylabel('Theta Norm Power(|Y|)')
title('Gamma Freq. v Theta Power')


figure('Color','w','Name','Freq v. Other Norm Power')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Theta Freq. (Hz)')
ylabel('Gamma Power(|Y|)')
title('Theta Freq. v Gamma Power')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Gamma Freq. (Hz)')
ylabel('Theta Power(|Y|)')
title('Gamma Freq. v Theta Power')


figure('Color','w','Name','CA3->Pyr v. Freq')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.conns.ca32pyr, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('CA3->Pyr')
ylabel('Theta Freq. (Hz)')
title('Theta Frequency')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.conns.ca32pyr, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('CA3->Pyr')
ylabel('Theta Power(|Y|)')
title('Theta Power')

figure('Color','w','Name','CA3->Pyr v Gamma')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.conns.ca32pyr, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('CA3->Pyr')
ylabel('Gamma Freq. (Hz)')
title('Gamma Frequency')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.conns.ca32pyr, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('CA3->Pyr')
ylabel('Gamma Power(|Y|)')
title('Gamma Power')


figure('Color','w','Name','Theta v Gamma')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Theta Freq. (Hz)')
ylabel('Gamma Freq. (Hz)')
title('Frequency')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Theta Power(|Y|)')
ylabel('Gamma Power(|Y|)')
title('Power')



figure('Color','w','Name','Stimulation v Norm Osc Power')
subplot(1,3,1)
for s=1:length(subsint)
    plot([myRuns(subsint(s).idx).DegreeStim],arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Stimulation (Hz)')
ylabel('Theta Norm Power(|Y|)')
title('Theta Norm Power')

subplot(1,3,2)
for s=1:length(subsint)
    plot([myRuns(subsint(s).idx).DegreeStim],arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Stimulation (Hz)')
ylabel('Gamma Norm Power(|Y|)')
title('Gamma Norm Power')

subplot(1,3,3)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Theta Norm Power(|Y|)')
ylabel('Gamma Norm Power(|Y|)')
title('Norm Power')








figure('Color','w','Name','Spikes v Power')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Spikes.Real, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('# Spikes')
ylabel('Theta Power |Y|')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Spikes.Real, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.power, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('# Spikes')
ylabel('Gamma Power |Y|')

figure('Color','w','Name','Spikes v Norm Power')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Spikes.Real, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('# Spikes')
ylabel('Theta Norm Power |Y|')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Spikes.Real, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('# Spikes')
ylabel('Gamma Norm Power |Y|')




for s=1:length(subsint)
    figure('Color','w','Name',['Histograms: ' subsint(s).desc])
    subplot(2,2,1)
    hist(arrayfun(@(x) x.Spikes.PyrFrac, myRuns(subsint(s).idx)))
    title('Pyramidal Spike Fraction')

    subplot(2,2,2)
    hist(arrayfun(@(x) x.Spikes.Real, myRuns(subsint(s).idx)))
    title('# Spikes')

    subplot(2,2,3)
    hist(arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)))
    title('Theta Freq.')

    subplot(2,2,4)
    hist(arrayfun(@(x) x.Gamma.freq, myRuns(subsint(s).idx)))
    title('Gamma Freq.')
end

for s=1:length(subsint)
    figure('Color','w','Name',['Power Histograms: ' subsint(s).desc])
    subplot(2,2,1)
    hist(arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)))
    title('Theta Power')

    subplot(2,2,2)
    hist(arrayfun(@(x) x.Gamma.power, myRuns(subsint(s).idx)))
    title('Gamma Power')

    subplot(2,2,3)
    hist(arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)))
    title('Theta Norm.')

    subplot(2,2,4)
    hist(arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx)))
    title('Gamma Norm.')
end

% 
% figure('Color','w','Name','Freq. v. Power')
% subplot(1,2,1)
% plot(arrayfun(@(x) x.Theta.freq, myRuns(idx2use)),arrayfun(@(x) x.Theta.power, myRuns(idx2use)),'LineStyle','none','Marker','.','MarkerSize',10)
% xlabel('Theta Freq. (Hz)')
% ylabel('Theta Power(|Y|)')
% title('Theta')
% 
% subplot(1,2,2)
% plot(arrayfun(@(x) x.Gamma.freq, myRuns(idx2use)),arrayfun(@(x) x.Gamma.power, myRuns(idx2use)),'LineStyle','none','Marker','.','MarkerSize',10)
% xlabel('Gamma Freq. (Hz)')
% ylabel('Gamma Power(|Y|)')
% title('Gamma')
% 
% 
% 
% figure('Color','w','Name','Freq. v. Norm Power')
% subplot(1,2,1)
% plot(arrayfun(@(x) x.Theta.freq, myRuns(idx2use)),arrayfun(@(x) x.Theta.norm, myRuns(idx2use)),'LineStyle','none','Marker','.','MarkerSize',10)
% xlabel('Theta Freq. (Hz)')
% ylabel('Theta Power(|Y|)')
% title('Theta')
% 
% subplot(1,2,2)
% plot(arrayfun(@(x) x.Gamma.freq, myRuns(idx2use)),arrayfun(@(x) x.Gamma.norm, myRuns(idx2use)),'LineStyle','none','Marker','.','MarkerSize',10)
% xlabel('Gamma Freq. (Hz)')
% ylabel('Gamma Power(|Y|)')
% title('Gamma')


figure('Color','w','Name','Power v. Norm Power')
subplot(1,2,1)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Theta Power(|Y|)')
ylabel('Theta Norm Power(|Y|)')
title('Theta')

subplot(1,2,2)
for s=1:length(subsint)
    plot(arrayfun(@(x) x.Gamma.power, myRuns(subsint(s).idx)),arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('Gamma Power(|Y|)')
ylabel('Gamma Norm Power(|Y|)')
title('Gamma')
return

for s=1:length(subsint)
    figure('Color','w','Name',['Hist. of Norm. Powers: ' subsint(s).desc]);
    hist(arrayfun(@(x) x.Gamma.norm, myRuns(subsint(s).idx))+arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)))
    xlabel('Theta + Gamma Norm Power(|Y|)')
    title('Norm Power Sum')
end


for z=1:length(myRuns(1).Cells)
    figure('Color','w','Name',['Theta Norm Power v. Cell Phase:' myRuns(1).Cells(z).name]);
    for s=1:length(subsint)
        plot(arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)),arrayfun(@(x) x.Cells(z).theta.phase, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
        hold on
    end
    legend({subsint.desc})
    xlabel('Theta Norm Power(|Y|)')
    ylabel([myRuns(1).Cells(z).name ' Theta Phase'])
end

for z=1:length(myRuns(1).Cells)
    figure('Color','w','Name',['Theta Norm Power v. Cell Mag:' myRuns(1).Cells(z).name]);
    for s=1:length(subsint)
        plot(arrayfun(@(x) x.Theta.norm, myRuns(subsint(s).idx)),arrayfun(@(x) x.Cells(z).theta.mag, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
        hold on
    end
    legend({subsint.desc})
    xlabel('Theta Norm Power(|Y|)')
    ylabel([myRuns(1).Cells(z).name ' Theta Mag'])
end

for z=1:length(myRuns(1).Cells)
    figure('Color','w','Name',['Theta Power v. Cell Phase:' myRuns(1).Cells(z).name]);
    for s=1:length(subsint)
        plot(arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),arrayfun(@(x) x.Cells(z).theta.phase, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
        hold on
    end
    legend({subsint.desc})
    xlabel('Theta Power(|Y|)')
    ylabel([myRuns(1).Cells(z).name ' Theta Phase'])
end

for z=1:length(myRuns(1).Cells)
    figure('Color','w','Name',['Theta Power v. Cell Mag:' myRuns(1).Cells(z).name]);
    for s=1:length(subsint)
        plot(arrayfun(@(x) x.Theta.power, myRuns(subsint(s).idx)),arrayfun(@(x) x.Cells(z).theta.mag, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
        hold on
    end
    legend({subsint.desc})
    xlabel('Theta Power(|Y|)')
    ylabel([myRuns(1).Cells(z).name ' Theta Mag'])
end


for z=1:length(myRuns(1).Cells)
    figure('Color','w','Name',['Theta Freq v. Cell Phase:' myRuns(1).Cells(z).name]);
    for s=1:length(subsint)
        plot(arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Cells(z).theta.phase, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
        hold on
    end
    legend({subsint.desc})
    xlabel('Theta Freq (Hz)')
    ylabel([myRuns(1).Cells(z).name ' Theta Phase'])
end

for z=1:length(myRuns(1).Cells)
    figure('Color','w','Name',['Theta Freq v. Cell Mag:' myRuns(1).Cells(z).name]);
    for s=1:length(subsint)
        plot(arrayfun(@(x) x.Theta.freq, myRuns(subsint(s).idx)),arrayfun(@(x) x.Cells(z).theta.mag, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
        hold on
    end
    legend({subsint.desc})
    xlabel('Theta Freq (Hz)')
    ylabel([myRuns(1).Cells(z).name ' Theta Mag'])
end

typevec=[1 2 3 4 5 6 8 9];
figure('Color','w','Name','Phase v Phase');
for z=1:length(typevec)
    for w=1:length(typevec)
        subplot(length(typevec),length(typevec),(z-1)*length(typevec)+w)
        for s=1:length(subsint)
            plot(arrayfun(@(x) x.Cells(typevec(z)).theta.phase, myRuns(subsint(s).idx)),arrayfun(@(x) x.Cells(typevec(w)).theta.phase, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
            hold on
        end
        xlim([0 360])
        ylim([0 360])
        xlabel(myRuns(1).Cells(typevec(z)).name)
        ylabel(myRuns(1).Cells(typevec(w)).name)
    end
end

for s=1:length(subsint)
    figure('Color','w','Name','Phase v Phase');
    for z=1:length(typevec)
        subplot(2,4,z)
        hist(arrayfun(@(x) x.Cells(typevec(z)).theta.phase, myRuns(subsint(s).idx)),[0:20:360])
        title(myRuns(1).Cells(typevec(z)).name)
    end
end

% 
% county = 0;
% for rs=1:length(idx2use)
%     if myRuns(idx2use(rs)).Theta.norm>=.5 && myRuns(idx2use(rs)).Gamma.norm>=.5 % && myRuns(idx2use(rs)).Spikes.PyrFrac<.4
%         disp(['rs: ' num2str(idx2use(rs)) ', Run: ', myRuns(idx2use(rs)).Name]);
%         county = county + 1;
%     end
% end
% county

figure;
for s=1:length(subsint)
    plot(arrayfun(@(x) x.conns.ca32pyr, myRuns(subsint(s).idx)), arrayfun(@(x) x.Spikes.PyrFrac, myRuns(subsint(s).idx)),'Color',subsint(s).col,'LineStyle','none','Marker','.','MarkerSize',10)
    hold on
end
legend({subsint.desc})
xlabel('CA3->Pyr')
ylabel('Pyramidal Spike Fraction')
