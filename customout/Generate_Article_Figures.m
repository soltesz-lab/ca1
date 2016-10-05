function Generate_Article_Figures(~,~,handles)
global mypath myFontSize sl colorvec myFontWeight myFontName myhandles repodir savepath printflag draftqual printtable fid


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Editable parameters here
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist([mypath sl 'data' sl 'myrepos.mat'],'file')
    load([mypath sl 'data' sl 'myrepos.mat'], 'myrepos')
end

q=getcurrepos(handles); %#ok<NASGU>

repodir=myrepos(q).dir;
tmp=findstr(repodir,sl);
savepath=[repodir(1:tmp(end)) 'papers' repodir(tmp(end):end) ];
if exist(savepath,'dir')==0
    mkdir(savepath)
end

ctrlstr='ca1_centerlfp_long_exc_065_01';%'ca1_nlfp_long_exc_065_01_02';

myElectrode=[2000 500];
maxdist=500;

whichFigs2Print=[6]; % Only figs 3-8 are implemented right now
% 1: Anatomical constraints 
% 2: Electrophysiological constraints
% 3: Spectrograms & Heat Maps
% 4: Spike raster, LFP and MP traces
% 5: Theta firing phases & modulation
% 6: Network manipulations including GABAB, also supplemental fig 6 freqs
% 7: Network Clamp
% 8: Model performance

% Supplemental figures
% 11: S1 - axonal distributions -TBA
% 12: S2 - SDF, LFP explanation
% 13: S3 - Fig 6 frequencies - not in operation. will print with argument
% of 6
% 14: S4 - Scaled up raster
% 15: S5 - GABA inset -TBA

% Placeholder figures
% 20: supplemental full raster -TBA
% 21: network clamp MP peak hist -TBA
% 22: edge effects -TBA

% If you want eps files printed into the 'savepath' directory, set 1. 
% Otherwise 0 for nothing printed.
% If you set the printflag, some figures will be consolidated into a single
% larger figure prior to printing to ensure the subfigures line up in the image.

printflag=0;
draftqual=1;
printtable=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('FinalFig6Data.mat','file')
    load('FinalFig6Data.mat','FinalData','cellnumbers','LayerHeights','LongitudinalLength','TransverseLength','historical','BenchData');
else
    disp('Put the FinalFig6Data.mat file in your SimTracker directory and try again.');
    return
end

updatectrlflag=1;%0;
if isempty(strmatch(FinalData(1).Name,ctrlstr))
    FinalData(1).Name=ctrlstr;
    FinalData(1).Comments='Control';
    FinalData(1).Idx=0;
    updatectrlflag=1;
end

usenorm=0;
myszht = 1.1;

myFontSize=8;
myFontWeight='normal';
myFontName = 'Arial';

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
              

            colorstruct.pyramidal.color=[.0 .0 .6];
            colorstruct.pvbasket.color=[.0 .75 .65];
            colorstruct.cck.color=[1 .75 .0];
            colorstruct.sca.color=[1 .5 .3];
            colorstruct.axoaxonic.color=[1 .0 .0];
            colorstruct.bistratified.color=[.6 .4 .1];
            colorstruct.olm.color=[.5 .0 .6];
            colorstruct.ivy.color=[.6 .6 .6];
            colorstruct.ngf.color=[1 .1 1];
            colorstruct.pyramidal.pos=1;
            colorstruct.pvbasket.pos=2;
            colorstruct.cck.pos=3;
            colorstruct.sca.pos=4;
            colorstruct.axoaxonic.pos=5;
            colorstruct.bistratified.pos=6;
            colorstruct.olm.pos=7;
            colorstruct.ivy.pos=8;
            colorstruct.ngf.pos=9;

mytypes=fieldnames(colorstruct);            
            
if ispc
    sl='\';
else
    sl='/';
end

load(['data' sl 'MyOrganizer.mat'],'general')


addpath outputtypes
addpath tools
rrrange=[];
 
missingflag=0;
for rr=1:length(FinalData) %#ok<NODEF>
    if isempty(FinalData(rr).Idx) || isempty(FinalData(rr).Pos)
        updatectrlflag=1;
        rrrange=[rrrange rr];
    end
    FinalData(rr).Idx=searchRuns('RunName',FinalData(rr).Name,0);
    if isempty(FinalData(rr).Idx)
        disp(['There is no run named ' FinalData(rr).Name ' in this SimTracker.'])
        missingflag=1;
    elseif length(FinalData(rr).Idx)>1
        disp(['Somehow there are multiple runs named ' FinalData(rr).Name ' and this is a problem but we will try to use the first one.'])
        FinalData(rr).Idx=FinalData(rr).Idx(1);
    end
end

if missingflag
    disp('Some runs are missing from your SimTracker. Import them and try again.')
    disp('Alternatively, you may be in a different repository. Switch back to ca1 and try again.')
    return
end


if isequal(historical.myElectrode,myElectrode) && isequal(historical.maxdist,maxdist)
    eval(['global ' ctrlstr])
    evalin('caller',['global ' ctrlstr])
    evalin('base',['global ' ctrlstr])
    if exist(ctrlstr,'var')==0 || eval(['isempty(' ctrlstr ')']) || eval(['isstruct(' ctrlstr ')==0'])
        % export again
        handles.curses=[];
        dd=strmatch(ctrlstr,{FinalData.Name},'exact');
        if isempty(dd)
            disp(['Cancelling...Could not find a run entitled "' ctrlstr '" within the FinalData list:'])
            disp({FinalData.Name}')
            return
        elseif dd~=1
            disp('Attempting to Continue, but Problem: control run is not first in FinalData list...')
        end
        handles.curses.ind=FinalData(dd).Idx;
        guidata(handles.btn_generate, handles)
        h=load_results([],handles);
    end
    if updatectrlflag
        rrrange=[rrrange 1];
    end
else
        rrrange=1:length(FinalData);
end
% load control data & analyze data
for rr=rrrange %#ok<NODEF>
    eval(['global ' FinalData(rr).Name])
    evalin('caller',['global ' FinalData(rr).Name])
    evalin('base',['global ' FinalData(rr).Name])
    
    if exist(FinalData(rr).Name,'var')==0 || eval(['isempty(' FinalData(rr).Name ')']) || eval(['isstruct(' FinalData(rr).Name ')==0'])
        % export again
        handles.curses=[];
        handles.curses.ind=FinalData(rr).Idx;
        guidata(handles.btn_generate, handles)
        h=load_results([],handles);
    end

    eval(['myresultsstruct = ' FinalData(rr).Name ';']);
    myresultsstruct.general=general;
    myresultsstruct.btn_generate = handles.btn_generate;
    if strcmp(ctrlstr,FinalData(rr).Name)==1
        myresultsstruct.optarg='pwelch,sdf,type,table';
    else
        myresultsstruct.optarg='pwelch,sdf,pyr,table';
    end

    if exist('myhandles','var')==0 || isempty(myhandles) || isstruct(myhandles)==0
        myhandles = CustomNearElectrode(myresultsstruct,cellnumbers,LayerHeights,LongitudinalLength,TransverseLength,myElectrode,maxdist);
    end
    myresultsstruct.curses.cells=myhandles.curses.cells;

    [~, tbldata]=plot_spectral(myresultsstruct);

    for t=1:size(tbldata,1)
        FinalData(rr).Pos.Freqs(t).Name=tbldata{t,1};
        FinalData(rr).Pos.Freqs(t).Theta=tbldata{t,2};
        FinalData(rr).Pos.Freqs(t).ThetaP=tbldata{t,3};
        FinalData(rr).Pos.Freqs(t).Gamma=tbldata{t,4};
        FinalData(rr).Pos.Freqs(t).GammaP=tbldata{t,5};
        FinalData(rr).Pos.Freqs(t).Overall=tbldata{t,6};
        FinalData(rr).Pos.Freqs(t).OverallP=tbldata{t,7};
    end

    FinalData(rr).Pyr=strmatch('Pyr.',tbldata(:,1));
end
historical.myElectrode=myElectrode;
historical.maxdist=maxdist;

%if ~isempty(rrrange)
    save('FinalFig6Data.mat','FinalData','cellnumbers','LayerHeights','LongitudinalLength','TransverseLength','historical','BenchData')
%end

dd=strmatch(ctrlstr,{FinalData.Name},'exact');
disp(['mycontrolstruct = ' ctrlstr ';'])
eval(['mycontrolstruct = ' ctrlstr ';'])
mycontrolstruct.general=general;
mycontrolstruct.timevec=FinalData(dd).timevec;

if exist('myhandles','var')==0 || isempty(myhandles) || isstruct(myhandles)==0
    myhandles = CustomNearElectrode(mycontrolstruct,cellnumbers,LayerHeights,LongitudinalLength,TransverseLength,myElectrode,maxdist);
end
mycontrolstruct.curses.cells=myhandles.curses.cells;
mycontrolstruct.btn_generate = handles.btn_generate;

if exist('mycontrolstruct','var')==0
    disp('No control results loaded')
    return
end

if sum(ismember(whichFigs2Print,[3:5]))>0
    otherCtrlFigs(FinalData,mycontrolstruct,colorstruct,whichFigs2Print)
end

if ismember(7,whichFigs2Print)
    [anglebase, anglepv, anglesom, h, mydata]=mynewphaseshiftLFP();    
    if printflag
        for t=1:length(h)
            printeps(h(t),[savepath sl get(h(t),'Name')])
        end
    end
end


if ismember(8,whichFigs2Print)
    performancefig(BenchData);    
end

if ismember(2,whichFigs2Print)
    h=ephysFig2();
    if printflag
        for t=1:length(h)
            printeps(h(t),[savepath sl get(h(t),'Name')])
        end
    end
end

if ismember(12,whichFigs2Print)
    h=SDFwalkthru(mycontrolstruct);
    if printflag
        for t=1:length(h)
            printeps(h(t),[savepath sl get(h(t),'Name')])
        end
    end
end


if sum(ismember(6,whichFigs2Print))==0
    return;
end

if printtable
    fid=fopen([savepath '..' sl 'Tables' sl 'Table_OscPower.tex'],'w');
    fprintf(fid,'\\begin{table}[position htb]\n');
    fprintf(fid,'\\centering\n');
    fprintf(fid,'\\begin{tabular}{|p{1.5in}|R{.85in}R{.85in}|R{.85in}R{.85in}|R{.85in}R{.85in}|} \n');
    fprintf(fid,'\\hline\n');
    fprintf(fid,' & \\multicolumn{2}{c|}{\\textbf{Theta}}& \\multicolumn{2}{c|}{\\textbf{Gamma}}& \\multicolumn{2}{c|}{\\textbf{Overall}} \\\\ \n');
    fprintf(fid,'\\hline\n');
    fprintf(fid,'\\multicolumn{1}{|c|}{\\textbf{Condition}} & \\multicolumn{1}{c}{\\textbf{Frequency}} & \\multicolumn{1}{c|}{\\textbf{Power}} & \\multicolumn{1}{c}{\\textbf{Frequency}} & \\multicolumn{1}{c|}{\\textbf{Power}} & \\multicolumn{1}{c}{\\textbf{Frequency}} & \\multicolumn{1}{c|}{\\textbf{Power}}  \\\\ \n');
    fprintf(fid,'\\hline\n');
end

GetAllSDFs(FinalData,handles)
%%% Everything below here is for generating figure 6:

GABAFig(FinalData,handles);
% Figure 6: Network Manipulations (Excitation, Reduced, Muted, Converge,

% Excitation Level
Runs2Use=[];
for r=1:length(FinalData)
    if ~isempty(strfind(FinalData(r).Comments,'Excitation')) || (~isempty(strfind(FinalData(r).Comments,'Control')) && isempty(strfind(FinalData(r).Comments,'GrrControl'))) %~isempty(find(myfigs(r).figs==2))
        Runs2Use=[Runs2Use r];
    end
end
exclevel=[FinalData(Runs2Use).DegreeStim];
[exclevel, sorti]=sort(exclevel);
Runs2Use=Runs2Use(sorti);
mdx=1;

thetafreq=[];
thetapow=[];
thetanorm=[];
for bb=1:length(Runs2Use)
    thetafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Theta;
    thetapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP;
    thetanorm(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP/FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
end

myname='DiffExcLevel';
mytitle='Tonic Excitation Level (Hz)';
myapp='';

if usenorm==1
    myplotvals=thetanorm;
else
    myplotvals=thetapow;
end
mysize=[length(Runs2Use) myszht];
xL={};
for r=1:length(Runs2Use)
    xL{r}=sprintf('%.2f',exclevel(r));
end

if usenorm==1
    yL={'Norm. Theta Power','Relative to Control'};
else
    yL='Theta Power';
end

mycolors=[];
myfreqs=[];
freqstruct=[];
for r=1:length(Runs2Use)
        mycolors(r,:)=[0 0 0];
        myfreqs(r,:)=[FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Theta FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).ThetaP FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Overall FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).OverallP  FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Gamma FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).GammaP];
        freqstruct(r).freqs= FinalData(Runs2Use(r)).Pos.Freqs(1).Gamma;%unique([FinalData(mdx).tbldata{:,6}]);
        freqstruct(r).pows= FinalData(Runs2Use(r)).Pos.Freqs(1).GammaP;%unique([FinalData(mdx).tbldata{:,6}]);
end
myplotvals(isnan(myplotvals))=0;
gy=plotmebaronly(myname, mycolors, myplotvals,mysize,xL,yL,myfreqs,freqstruct,mytitle,myapp);
if printflag
printeps(gy(1),[savepath sl myname])
printeps(gy(2),[savepath sl myname 'Freq'])
end



% Ephys same
Runs2Use=[];
for r=1:length(FinalData)
    if ~isempty(strfind(FinalData(r).Comments,'Reduced')) || (~isempty(strfind(FinalData(r).Comments,'Control')) && isempty(strfind(FinalData(r).Comments,'GrrControl'))) %~isempty(find(myfigs(r).figs==2))
        Runs2Use=[Runs2Use r];
    end
end

thetafreq=[];
thetapow=[];
thetanorm=[];
for bb=1:length(Runs2Use)
    thetafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Theta;
    thetapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP;
    thetanorm(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP/FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
end
    
myname='DiffCellEphys';
mytitle='Single Interneuron E''phys. Profile';
myapp='';

if usenorm==1
    myplotvals=thetanorm;
else
    myplotvals=thetapow;
end

mysize=[length(Runs2Use) myszht];
xL={};
ttL={};
for r=1:length(Runs2Use)
    if ~isempty(strfind(FinalData(Runs2Use(r)).Comments,'Control')) && isempty(strfind(FinalData(Runs2Use(r)).Comments,'GrrControl'))
        xL{r}='Ctrl';
        ttL{r}='Ctrl';
    else
        ttmp=FinalData(Runs2Use(r)).Name(9:end-3);
        ttL{r}=ttmp;
        rr=strmatch(ttmp,{'OLM','CCK','PV','NGF'});
        nicett={'O-LM','CCK+B','PV+B','NGF'};
        xL{r}=nicett{rr};
    end
end

if usenorm==1
    yL={'Norm. Theta Power','Relative to Control'};
else
    yL='Theta Power';
end

mycolors=[];
myfreqs=[];
freqstruct=[];
for r=1:length(Runs2Use)
    tt=strmatch(lower(ttL{r}),mytypes);
    if isempty(tt)
        mycolors(r,:)=[0 0 0];
    else
        mycolors(r,:)=colorstruct.(mytypes{tt}).color;
    end
        myfreqs(r,:)=[FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Theta FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).ThetaP FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Overall FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).OverallP  FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Gamma FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).GammaP];
        freqstruct(r).freqs= FinalData(Runs2Use(r)).Pos.Freqs(1).Gamma;%unique([FinalData(mdx).tbldata{:,6}]);
        freqstruct(r).pows= FinalData(Runs2Use(r)).Pos.Freqs(1).GammaP;%unique([FinalData(mdx).tbldata{:,6}]);
end
myplotvals(isnan(myplotvals))=0;
gy=plotmebaronly(myname, mycolors, myplotvals,mysize,xL,yL,myfreqs,freqstruct,mytitle,myapp);
if printflag
printeps(gy(1),[savepath sl myname])
printeps(gy(2),[savepath sl myname 'Freq'])
end

% Ephys same and PV Conv
Runs2Use=[];
for r=1:length(FinalData)
    if ~isempty(strfind(FinalData(r).Comments,'Conv PV')) || (~isempty(strfind(FinalData(r).Comments,'Control')) && isempty(strfind(FinalData(r).Comments,'GrrControl'))) %~isempty(find(myfigs(r).figs==2))
        Runs2Use=[Runs2Use r];
    end
end

thetafreq=[];
thetapow=[];
thetanorm=[];
for bb=1:length(Runs2Use)
    thetafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Theta;
    thetapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP;
    thetanorm(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP/FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
end

myname='Converge2PVB';
mytitle='Interneurons Converge to PV+ B. cells';
myapp='';

if usenorm==1
    myplotvals=thetanorm;
else
    myplotvals=thetapow;
end

mysize=[length(Runs2Use) myszht];
xL={};
myconvname={'E''phys.','+input wgt','+input #','All PV+B'};
chkr=1;
for r=1:length(Runs2Use)
    if ~isempty(strfind(FinalData(Runs2Use(r)).Comments,'Control')) && isempty(strfind(FinalData(Runs2Use(r)).Comments,'GrrControl'))
        xL{r}='Ctrl';
        chkr=chkr-1;
    else
        xL{r}=myconvname{chkr};%char(str2num(FinalData(Runs2Use(r)).Name(end))-1+'A');
    end
    chkr=chkr+1;
end

%xL={myfigs(Runs2Use).label};
%xL={myfigs(1).mydata.name};
if usenorm==1
yL={'Norm. Theta Power','Relative to Control'};
else
yL='Theta Power';
end

mycolors=[];
myfreqs=[];
freqstruct=[];
for r=1:length(Runs2Use)
    mycolors(r,:)=[0 0 0];
        myfreqs(r,:)=[FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Theta FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).ThetaP FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Overall FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).OverallP  FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Gamma FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).GammaP];
        freqstruct(r).freqs= FinalData(Runs2Use(r)).Pos.Freqs(1).Gamma;%unique([FinalData(mdx).tbldata{:,6}]);
        freqstruct(r).pows= FinalData(Runs2Use(r)).Pos.Freqs(1).GammaP;%unique([FinalData(mdx).tbldata{:,6}]);
end
myplotvals(isnan(myplotvals))=0;
gy=plotmebaronly(myname, mycolors, myplotvals,mysize,xL,yL,myfreqs,freqstruct,mytitle,myapp);
if printflag
printeps(gy(1),[savepath sl myname])
printeps(gy(2),[savepath sl myname 'Freq'])
end

% Disinh
Runs2Use=[];
for r=1:length(FinalData)
    if ~isempty(strfind(FinalData(r).Comments,'Disinh')) || (~isempty(strfind(FinalData(r).Comments,'Control')) && isempty(strfind(FinalData(r).Comments,'GrrControl'))) %~isempty(find(myfigs(r).figs==2))
        Runs2Use=[Runs2Use r];
    end
end

thetafreq=[];
thetapow=[];
thetanorm=[];
for bb=1:length(Runs2Use)
    thetafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Theta;
    thetapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP;
    thetanorm(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP/FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
end

myname='Disinh';
mytitle='Disinhibition';
myapp='';

if usenorm==1
    myplotvals=thetanorm;
else
    myplotvals=thetapow;
end

mysize=[length(Runs2Use) myszht];
xL={};
for r=1:length(Runs2Use)
    if ~isempty(strfind(FinalData(Runs2Use(r)).Comments,'Control')) && isempty(strfind(FinalData(Runs2Use(r)).Comments,'GrrControl'))
        xL{r}='Ctrl';
    else
        whichnum=regexp(FinalData(Runs2Use(r)).Name,'[0-9][0-9]','match');
        xL{r}= [FinalData(Runs2Use(r)).Name(14:end-4) '-' whichnum{1} '%'];%strrep(strrep(FinalData(Runs2Use(r)).Comments(7:end),'.',''),' ','');
    end
end
if usenorm==1
yL={'Norm. Theta Power','Relative to Control'};
else
yL='Theta Power';
end

mycolors=[];
myfreqs=[];
freqstruct=[];
for r=1:length(Runs2Use)
    mycelltype=lower(strrep(strrep(strrep(strrep(xL{r},'-',''),'.',''),'+',''),' ',''));
    if strcmp(xL{r}(1:2),'pv')
        mycolors(r,:)=[0 1 0];
    elseif strcmp(xL{r}(1:2),'so')
        mycolors(r,:)=[.8 .2 1];
    else
        mycolors(r,:)=[0 0 0];
    end
    myfreqs(r,:)=[FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Theta FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).ThetaP FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Overall FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).OverallP  FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Gamma FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).GammaP];
    freqstruct(r).freqs= FinalData(Runs2Use(r)).Pos.Freqs(1).Gamma;%unique([FinalData(mdx).tbldata{:,6}]);
    freqstruct(r).pows= FinalData(Runs2Use(r)).Pos.Freqs(1).GammaP;%unique([FinalData(mdx).tbldata{:,6}]);
end
myplotvals(isnan(myplotvals))=0;
gy=plotmebaronly(myname, mycolors, myplotvals,mysize,xL,yL,myfreqs,freqstruct,mytitle,myapp);
if printflag
printeps(gy(1),[savepath sl myname])
printeps(gy(2),[savepath sl myname 'Freq'])
end


% Recurrent Collaterals
Runs2Use=[];
for r=1:length(FinalData)
    if ~isempty(strfind(FinalData(r).Comments,'RecColl')) || (~isempty(strfind(FinalData(r).Comments,'Control')) && isempty(strfind(FinalData(r).Comments,'GrrControl'))) %~isempty(find(myfigs(r).figs==2))
        Runs2Use=[Runs2Use r];
    end
end

thetafreq=[];
thetapow=[];
thetanorm=[];
for bb=1:length(Runs2Use)
    thetafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Theta;
    thetapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP;
    thetanorm(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP/FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
end

myname='RecColl';
mytitle='Recurrent Collaterals';
myapp='';

if usenorm==1
    myplotvals=thetanorm;
else
    myplotvals=thetapow;
end

mysize=[length(Runs2Use) myszht];
xL={};
for r=1:length(Runs2Use)
    if ~isempty(strfind(FinalData(Runs2Use(r)).Comments,'Control')) && isempty(strfind(FinalData(Runs2Use(r)).Comments,'GrrControl'))
        xL{r}='Ctrl';
    else
        xL{r}= FinalData(Runs2Use(r)).Comments(length(RecColl)+1:end);%strrep(strrep(FinalData(Runs2Use(r)).Comments(7:end),'.',''),' ','');
    end
end
if usenorm==1
yL={'Norm. Theta Power','Relative to Control'};
else
yL='Theta Power';
end

mycolors=[];
myfreqs=[];
freqstruct=[];
for r=1:length(Runs2Use)
    mycelltype=lower(strrep(strrep(strrep(strrep(xL{r},'-',''),'.',''),'+',''),' ',''));
    mycolors(r,:)=[0 0 0];
    myfreqs(r,:)=[FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Theta FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).ThetaP FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Overall FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).OverallP  FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Gamma FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).GammaP];
    freqstruct(r).freqs= FinalData(Runs2Use(r)).Pos.Freqs(1).Gamma;%unique([FinalData(mdx).tbldata{:,6}]);
    freqstruct(r).pows= FinalData(Runs2Use(r)).Pos.Freqs(1).GammaP;%unique([FinalData(mdx).tbldata{:,6}]);
end
myplotvals(isnan(myplotvals))=0;
gy=plotmebaronly(myname, mycolors, myplotvals,mysize,xL,yL,myfreqs,freqstruct,mytitle,myapp);
if printflag
printeps(gy(1),[savepath sl myname])
printeps(gy(2),[savepath sl myname 'Freq'])
end

% Muted
Runs2Use=[];
for r=1:length(FinalData)
    if ~isempty(strfind(FinalData(r).Comments,'Muted')) || (~isempty(strfind(FinalData(r).Comments,'Control')) && isempty(strfind(FinalData(r).Comments,'GrrControl'))) %~isempty(find(myfigs(r).figs==2))
        Runs2Use=[Runs2Use r];
    end
end

thetafreq=[];
thetapow=[];
thetanorm=[];
for bb=1:length(Runs2Use)
    thetafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Theta;
    thetapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP;
    thetanorm(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP/FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
end

myname='2All';
mytitle='Outputs Muted From Each Cell Type';
myapp='';

if usenorm==1
    myplotvals=thetanorm;
else
    myplotvals=thetapow;
end

mysize=[length(Runs2Use) myszht];
xL={};
for r=1:length(Runs2Use)
    if ~isempty(strfind(FinalData(Runs2Use(r)).Comments,'Control')) && isempty(strfind(FinalData(Runs2Use(r)).Comments,'GrrControl'))
        xL{r}='Ctrl';
    else
        xL{r}=strrep(strrep(FinalData(Runs2Use(r)).Comments(7:end),'.',''),' ','');
    end
end
if usenorm==1
yL={'Norm. Theta Power','Relative to Control'};
else
yL='Theta Power';
end

mycolors=[];
myfreqs=[];
freqstruct=[];
for r=1:length(Runs2Use)
    mycelltype=lower(strrep(strrep(strrep(strrep(xL{r},'-',''),'.',''),'+',''),' ',''));
    tt=strmatch(mycelltype(1:3),mytypes);
    if isempty(tt)
        mycolors(r,:)=[0 0 0];
    else
        mycolors(r,:)=colorstruct.(mytypes{tt}).color;
    end
        myfreqs(r,:)=[FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Theta FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).ThetaP FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Overall FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).OverallP  FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Gamma FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).GammaP];
        freqstruct(r).freqs= FinalData(Runs2Use(r)).Pos.Freqs(1).Gamma;%unique([FinalData(mdx).tbldata{:,6}]);
        freqstruct(r).pows= FinalData(Runs2Use(r)).Pos.Freqs(1).GammaP;%unique([FinalData(mdx).tbldata{:,6}]);
end
myplotvals(isnan(myplotvals))=0;
gy=plotmebaronly(myname, mycolors, myplotvals,mysize,xL,yL,myfreqs,freqstruct,mytitle,myapp);
if printflag
printeps(gy(1),[savepath sl myname])
printeps(gy(2),[savepath sl myname 'Freq'])
end

% PV Exc
Runs2Use=[];
for r=1:length(FinalData)
    if ~isempty(strfind(FinalData(r).Comments,'Exc PV')) %|| ~isempty(strfind(FinalData(r).Comments,'Control')) %~isempty(find(myfigs(r).figs==2))
        Runs2Use=[Runs2Use r];
    end
end
exclevel=[FinalData(Runs2Use).DegreeStim];
[exclevel, sorti]=sort(exclevel);
Runs2Use=[Runs2Use(sorti) mdx];
mdx=1;

    thetafreq=[];
    thetapow=[];
    thetanorm=[];
    
    for bb=1:length(Runs2Use)
        thetafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Theta;
        thetapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP;
        thetanorm(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP/FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
    end

myname='DiffPVExcLevel';
mytitle='Pyr & PV+ B. Network: Tonic Excitation Level (Hz)';
myapp='';

if usenorm==1
    myplotvals=thetanorm;
else
    myplotvals=thetapow;
end
mysize=[length(Runs2Use) myszht];
xL={};
for r=1:length(Runs2Use)-1
    xL{r}=sprintf('%.2f',exclevel(r));
end
xL{length(Runs2Use)}='Ctrl 0.65';

if usenorm==1
yL={'Norm. Theta Power','Relative to Control'};
else
yL='Theta Power';
end

mycolors=[];
myfreqs=[];
freqstruct=[];
for r=1:length(Runs2Use)
    %tt=strmatch(myfigs(Runs2Use(r)).label,{myfigs(1).mydata.name});
    %if isempty(tt)
        mycolors(r,:)=[0 0 0];
    %else
    %    mycolors(r,:)=myfigs(1).mydata(tt).color;
    %end

        myfreqs(r,:)=[FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Theta FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).ThetaP FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Overall FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).OverallP  FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).Gamma FinalData(Runs2Use(r)).Pos.Freqs(FinalData(Runs2Use(r)).Pyr).GammaP];
        freqstruct(r).freqs= FinalData(Runs2Use(r)).Pos.Freqs(1).Gamma;%unique([FinalData(mdx).tbldata{:,6}]);
        freqstruct(r).pows= FinalData(Runs2Use(r)).Pos.Freqs(1).GammaP;%unique([FinalData(mdx).tbldata{:,6}]);
end
myplotvals(isnan(myplotvals))=0;
gy=plotmebaronly(myname, mycolors, myplotvals,mysize,xL,yL,myfreqs,freqstruct,mytitle,myapp);
if printflag
printeps(gy(1),[savepath sl myname])
printeps(gy(2),[savepath sl myname 'Freq'])
end


if printtable
    fprintf(fid,'\\end{tabular}\n');
    fprintf(fid,'\\caption[Spectral Analysis of Network Activity]{Peak, theta and gamma frequencies and powers of the pyramidal cell spike density function using Welch''s Periodogram.}\n');
    fprintf(fid,'\\label{tab:oscpower}\n');
    fprintf(fid,'\\end{table}\n');
    fclose(fid);
end


function formatter(ax,varargin)
global myFontSize myFontName myFontWeight
if isempty(varargin)
    set(ax,'LineWidth',1.0,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)  
elseif varargin{1}==0
    set(ax,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)  
else
    set(ax,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)  
end
box off

function GetAllSDFs(FinalData,handles)
global savepath sl
    load('FinalFig6Data.mat','cellnumbers','LayerHeights','LongitudinalLength','TransverseLength','historical')
load(['data' sl 'MyOrganizer.mat'],'general')
handles.general=general;
timewindow=2000;
timeflag=1;
myElectrode=historical.myElectrode ;%[2000 500];
maxdist=historical.maxdist ;%500;
fid=fopen([savepath sl 'Pyramidal_SDF_All_Conditions.txt'],'w');
fid2=fopen([savepath sl 'All_SDF_Ctrl_Condition.txt'],'w');
fprintf(fid,'Time');
fprintf(fid2,'Time');
for r=1:length(FinalData)-1
    if r==1 && isfield(FinalData(1),'sdftimevec') && ~isempty(FinalData(1).sdftimevec)
        mymat=zeros(length(FinalData(1).sdftimevec),length(FinalData)+1);
        onemat=zeros(length(FinalData(1).onetimevec),length(FinalData)+1);
        mymat(:,1)=FinalData(r).sdftimevec(:);
        onemat(:,1)=FinalData(r).onetimevec(:);
        onemat(:,r+1)=FinalData(r).tmpsdf(:);
        mymat(:,r+1)=FinalData(r).sdf(:);
    elseif r>1 &&  isfield(FinalData(r),'sdf') && ~isempty(FinalData(r).sdf) &&  isfield(FinalData(r),'tmpsdf') && ~isempty(FinalData(r).tmpsdf)
        onemat(:,r+1)=FinalData(r).tmpsdf(:);
        mymat(:,r+1)=FinalData(r).sdf(:);
    else

    if 1==1 %infig6
        eval(['global ' FinalData(r).Name])
        evalin('caller',['global ' FinalData(r).Name])
        evalin('base',['global ' FinalData(r).Name])
        if exist(FinalData(r).Name,'var')==0 || eval(['isempty(' FinalData(r).Name ')']) || eval(['isstruct(' FinalData(r).Name ')==0'])
            % export again
            handles.curses=[];
            handles.curses.ind=FinalData(r).Idx;
            guidata(handles.btn_generate, handles)
            h=load_results([],handles);
        end
        eval(['handles.curses=' FinalData(r).Name '.curses;']);
        eval(['handles.runinfo=' FinalData(r).Name '.runinfo;']);
        if exist('myhandles','var')==0 || isempty(myhandles) || isstruct(myhandles)==0
            myhandles = CustomNearElectrode(handles,cellnumbers,LayerHeights,LongitudinalLength,TransverseLength,myElectrode,maxdist);
        end
        handles.curses.cells=myhandles.curses.cells;
        
        fprintf(fid,'\t%s',FinalData(r).Name);
        
        pyrIdx=7;
        idx = find(ismember(handles.curses.spikerast(:,2),handles.curses.cells(pyrIdx).mygids)==1);
        if timeflag
            idt = find(handles.curses.spikerast(idx,1)>handles.general.crop & handles.curses.spikerast(idx,1)<timewindow);
        else
            idt = find(handles.curses.spikerast(idx,1)>handles.general.crop);
        end

        if timeflag
            binned= histc(handles.curses.spikerast(idx(idt),1),[handles.general.crop:1:timewindow]); % binned by 1 ms
            timevec=[handles.general.crop:1:timewindow];
        else
            binned= histc(handles.curses.spikerast(idx(idt),1),[handles.general.crop:1:handles.runinfo.SimDuration]); % binned by 1 ms
            timevec=[handles.general.crop:1:handles.runinfo.SimDuration];
        end
        window=3; % ms
        kernel=mynormpdf(-floor((window*6+1)/2):floor((window*6+1)/2),0,window);
        sdf = conv(binned,kernel,'same');
        if r==1
            mymat=zeros(length(timevec),length(FinalData)+1);
            mymat(:,1)=timevec(:);
            FinalData(r).sdftimevec=timevec(:);
            
            for q=1:length(handles.curses.cells)
                fprintf(fid2,'\t%s',handles.curses.cells(q).name);
                idx = find(ismember(handles.curses.spikerast(:,2),handles.curses.cells(q).mygids)==1);
                if timeflag
                    idt = find(handles.curses.spikerast(idx,1)>handles.general.crop & handles.curses.spikerast(idx,1)<timewindow);
                else
                    idt = find(handles.curses.spikerast(idx,1)>handles.general.crop);
                end

                if timeflag
                    binned= histc(handles.curses.spikerast(idx(idt),1),[handles.general.crop:1:timewindow]); % binned by 1 ms
                    timevec=[handles.general.crop:1:timewindow];
                else
                    binned= histc(handles.curses.spikerast(idx(idt),1),[handles.general.crop:1:handles.runinfo.SimDuration]); % binned by 1 ms
                    timevec=[handles.general.crop:1:handles.runinfo.SimDuration];
                end
                window=3; % ms
                kernel=mynormpdf(-floor((window*6+1)/2):floor((window*6+1)/2),0,window);
                tmpsdf = conv(binned,kernel,'same');                
                if q==1
                    onemat=zeros(length(timevec),length(FinalData)+1);
                    onemat(:,1)=timevec(:);
                    FinalData(r).onetimevec=timevec(:);
                end
                onemat(:,r+1)=tmpsdf(:);%onemat=[onemat tmpsdf(:)];
                FinalData(r).tmpsdf=tmpsdf(:);
                save('FinalFig6Data.mat','FinalData','-append')
            end
        end
        mymat(:,r+1)=sdf(:);%mymat=[mymat sdf(:)];
        FinalData(r).sdf=sdf(:);
        save('FinalFig6Data.mat','FinalData','-append')
    end
    end
end
fclose(fid);
fclose(fid2);
%fprintf(fid,'\n');
dlmwrite([savepath sl 'Pyramidal_SDF_All_Conditions.txt'],mymat,'-append',...
'delimiter','\t','roffset',1)
dlmwrite([savepath sl 'All_SDF_Ctrl_Condition.txt'],onemat,'-append',...
'delimiter','\t','roffset',1)


function h=printExpData(myresultsstruct)
global colorvec myFontSize myFontWeight myFontName sl
showmod = 0;
try
        ftmp = figure('GraphicsSmoothing','off', 'Renderer', 'painters','Color','w','Name','Theta ExpHist','Units','inches','Position',[.1 .1 1.75 6],'PaperUnits','inches','PaperPosition',[0 0 1.75 6],'PaperSize',[1.75 6]);
catch
        ftmp = figure('Renderer', 'painters','Color','w','Name','Theta ExpHist','Units','inches','Position',[.1 .1 1.75 6],'PaperUnits','inches','PaperPosition',[0 0 1.75 6],'PaperSize',[1.75 6]);
end
        nrninput = gettheta(-2);
        nrninputf = gettheta(0);
        
        a1=subplot('Position',[.05 1-1/10+.01 .9 1/10-.04]);
        nrninputmp=nrninput;
        for fg=1:length(nrninputmp)
            nrninputmp(fg).name=nrninputmp(fg).tech;
        end
        %h=Phases_figure('GraphicsSmoothing','off', 'Renderer', 'painters',nrninputmp,myRuns(myridx).Theta.freq,gca);
        period = 125;
        Hzval = 8;
        trace.data=0:.025:(period*2);
        trace.data=trace.data';

        excell = -sin((Hzval*(2*pi))*trace.data(:,1)/1000 + pi/2);  % -13.8/125);  %  - handles.phasepref +   -  (0.25)*Hzval*2*pi
        plot(trace.data,excell,'k','LineWidth',1.5)
        ylim([-1.1 1.1])
        xlim([0 max(trace.data)])
        set(gca,'XTickLabel',{})
        set(gca,'YTickLabel',{})
        set(gca,'XColor','w')
        set(gca,'YColor','w')
        title('Experimental','FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
        %celltypeNice={'Pyr.','PV+ B.','CCK+ B.','S.C.-A.','Axo.','Bis.','O-LM','Ivy','NGF.'};
        celltypeNice={'Pyr','PV+B','CCK+B','SC-A','Axo','Bis','O-LM','Ivy','NGF'};
        for r=1:length(myresultsstruct.curses.cells)
            celltype = myresultsstruct.curses.cells(r).name; %celltypevec{r};
            z = strmatch(celltype,{'pyramidalcell','pvbasketcell','cckcell','scacell','axoaxoniccell','bistratifiedcell','olmcell','ivycell','ngfcell'});
            if isempty(z)
                continue;
            end
                subplot('Position',[.05 1-(z+1)/10+.05+.01+.02 .72 .03]);
                % subplot('Position',[.05 1-marg-(z+1)/10+(1/10-.05) .72 1/10-.05]);
                axis off
                b(r) = text(0,.45,celltypeNice{z});
                if strcmp(celltype,'scacell')==1
                    set(b(r),'Color','w','FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
                else
                    set(b(r),'Color',colorvec(z,:),'FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
                end
%                 else
%                     set(b(r),'Color','w','FontWeight',myFontWeight,'FontName',myFontName)
%                 end
                subplot('Position',[.95 1-(z+1)/10+.05+.01+.02 .04 .03]);
                % subplot('Position',[.95 1-marg-(z+1)/10+(1/10-.05) .04 1/10-.05]);
                myz=strmatch(celltype,{nrninput.tech},'exact');
                if isempty(myz)
                    myz=strmatch('cckcell',{nrninput.name},'exact');
                end
                if isempty(myz)
                    b(r) = text(0,.5,[num2str(round(0)) '^o']);
                    set(b(r),'HorizontalAlignment','right')
                    axis off
                    set(b(r),'Color','w','FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
                elseif showmod==1
                    myx=strmatch(celltype,{nrninputf.tech});
                    mystr='';
                    if ~isempty(myx)
                    mystr=sprintf('%0.0f^o,',unique([nrninputf(myx).phase]));
                    mystr=['(awake='  mystr(1:end-1)  ') '];
                    end
                    b(r) = text(0,.5,[mystr  num2str(round(nrninput(myz).phase)) '^o']);
                    set(b(r),'HorizontalAlignment','right')
                    axis off
                    set(b(r),'Color',colorvec(z,:),'FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
                else
                    % myx=strmatch(celltype,{nrninput.tech});
                    if ~isempty(myz)
                    b(r) = text(0,.5,[ num2str(round(nrninput(myz).phase)) '^o']);
                    set(b(r),'HorizontalAlignment','right')
                    set(b(r),'Color',colorvec(z,:),'FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
                    else
                        'uh oh'
                    end
                    axis off
                end

            subplot('Position',[.05 1-(z+1)/10+.01 .9 .05+.015])
            %myim=imread(['C:\Users\maria_000\Documents\MATLAB\work\SimTracker\KS08\' celltype(1:end-4) '.bmp']);
            if exist(['.' sl 'KS08' sl celltype(1:end-4) '.bmp'],'file')
                myim=imread(['.' sl 'KS08' sl celltype(1:end-4) '.bmp']);
                image(myim)
                hold on
                xL=get(gca,'xLim');
                yL=get(gca,'yLim');
                fst = (xL(2) - xL(1))/4;
                snd= (xL(2) - xL(1))*3/4;
                if strcmp(celltype(1:end-4),'sca')
                    plot([fst fst], yL, 'Color',[1 1 1],'LineStyle','--')
                    plot([snd snd], yL, 'Color',[1 1 1],'LineStyle','--')
                    plot(xL, [yL(end) yL(end)],'Color',[1 1 1])
                else
                    plot([fst fst], yL, 'Color',[.5 .5 .5],'LineStyle','--')
                    plot([snd snd], yL, 'Color',[.5 .5 .5],'LineStyle','--')
                    plot(xL, [yL(end) yL(end)],'Color',[.75 .75 .75])
                end
                
%                 hold on
%                 xr=xlim();
%                 yr=ylim();
%                 plot([diff(xr)*.25+xr(1) diff(xr)*.25+xr(1)],yr)
%                 plot([xr(2)-diff(xr) xr(2)-diff(xr)],yr)
                
                axis off
            end
        end
        tt=MiniPhaseShift_Figure(nrninputmp,8,[]);
        ylim([-1.1-.8*(9-1) 1.1])
        bf = findall(tt,'Type','text');
        for b=1:length(bf)
            set(bf(b),'FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
        end
        bf = findall(tt,'Type','axis');
        for b=1:length(bf)
            set(bf(b),'FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
        end
        h=[ftmp tt];

function myh=tgcross(handles)
global myFontName myFontWeight myFontSize
%plot(handles.curses.lfp(:,1),handles.curses.lfp(:,2))
        mydata=handles.curses.lfp;
        
        
 dt=diff(mydata(1:2,1))/1000;
d=mydata(:,2)';  % -1 to convert to pyramidal layer LFP

t = (1:length(d))*dt;               %Define a time axis.
N = length(t);  Fs = 1/dt;          %Define other useful parameters.

if 1==0
deg=2;                              %Sets the filter order.
LowWn = [5*2/Fs,10*2/Fs]  ;             %Define the low frequency window of interest.
[B,A] = butter(deg,LowWn,'bandpass');  %Apply the filter to isolate the band.
dlo = filtfilt(B,A,d);

HighWn = [25*2/Fs,80*2/Fs] ;           %Define the high frequency window of interest.
[B,A] = butter(deg,HighWn,'bandpass');  %Apply the filter to isolate the band.
dhi = filtfilt(B,A,d);

%   How well does the filter work?  Let's compare the original signal with
%   the fitlered data.
else
    dlo=mikkofilter(handles.curses.lfp,1000/handles.runinfo.lfp_dt,[5 10]);
    dhi=mikkofilter(handles.curses.lfp,1000/handles.runinfo.lfp_dt,[25 80]);
    t=dlo(:,1);
    dlo=dlo(:,2);
    dhi=dhi(:,2);
end

%   We might notice that funny effects happen at the edges of the data.  So,
%   let's ignore these edge effects due to the filtering and concentrate on
%   the center of the data.

dlo = dlo(round(N/4):end-round(N/4)-1);
dhi = dhi(round(N/4):end-round(N/4)-1);
taxis = t(round(N/4):end-round(N/4)-1);

phi = angle(hilbert(dlo));    %Phase of low frequency signal.
amp = abs(hilbert(dhi));      %Amplitude envelope of high frequency signal.

p_bins = (-pi:0.2:pi);              %Define the phase bins.

a_mean = zeros(length(p_bins)-1,1); %Vector to hold avg amp env results.
theta_mean = zeros(length(p_bins)-1,1); %Vector to hold avg amp env results.
p_mean = zeros(length(p_bins)-1,1); %Vector to hold center of phase bins.
  
for k=1:length(p_bins)-1
    pL = p_bins(k);                         %Phase lower limit for this bin.
    pR = p_bins(k+1);                       %Phase upper limit for this bin.
    indices = find(phi >= pL & phi < pR);   %Find phase values in this range.
    a_mean(k) = mean(amp(indices));         %Compute mean amplitude at these phases.
    theta_mean(k) = -mean(amp(indices));         %Compute mean amplitude at these phases.
    p_mean(k) = mean([pL, pR]);             %Label the phase bin with the center phase.
end
h = max(a_mean)-min(a_mean);                %The difference between the max and min modulation.
[R,P] = corrcoef(phi,amp);
%   Plot the mean envelope versus phase.
myh=figure('Color','w','Name','ThetaGammaCross','Units','inches','PaperUnits','inches','Position',[1 1 3.75 3],'PaperPosition',[0 0 3.75 3],'PaperSize',[3.75 3]);
% subplot(2,1,1)
% plot([p_mean; p_mean+2*pi]*180/pi+180, [theta_mean; theta_mean], 'k', 'LineWidth', 2);
% set(gca,'XTick',[0:90:720],'xlim',[0 720])
% xlabel('Theta Phase (degrees)');ylabel('calculated LFP analog (mV)')
% subplot(2,1,2)
plot([p_mean; p_mean+2*pi]*180/pi+180, [a_mean; a_mean], 'k', 'LineWidth', 2);
xlabel('Theta Phase (^o)');ylabel('Gamma Amplitude (mV)')
bf = findall(myh,'Type','text');
for b=1:length(bf)
    set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
end
bf = findall(myh,'Type','axis');
for b=1:length(bf)
    set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
end
set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
box off
set(gca,'XTick',[0:90:720],'xlim',[0 720])

%xlabel('Low frequency phase');  ylabel('High frequency envelope height difference');
disp(['crosscorr Metric h=' num2str(h)])


function otherCtrlFigs(FinalData,mycontrolstruct,colorstruct,plotCtrlFigs)
global myFontSize sl colorvec myFontWeight myFontName myhandles repodir savepath printflag RunArray

mrr=ver;

if ismember(3,plotCtrlFigs)
% Figure 3a: Spectrogram
mycontrolstruct.optarg='lfp';
h(1)=plot_spectro(mycontrolstruct);
set(h(1),'Name','Spectrogram')
title('LFP')
colormap(jet)
cbar_axes=colorbar;
set(get(cbar_axes,'ylabel'),'string','Power')
bf = findall(h(1),'Type','text');
for b=1:length(bf)
    set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
end
bf = findall(h(1),'Type','axis');
for b=1:length(bf)
    set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
end
set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
ylim([0 15])
set(h(1),'Units','inches','PaperUnits','inches','PaperSize',[3.75 3],'PaperPosition',[0 0 3.75 3],'Position',[.5 .5 3.75 3])

% Figure 3b: SDF Pwelch Periodogram
    mycontrolstruct.optarg='sdf,type,heat,norm,pwelch,[0 16]';
    h(2)=plot_spectral(mycontrolstruct);
    set(h(2),'Name','ThetaSpect')
    title('Theta')
    bf = findall(h(2),'Type','text');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    bf = findall(h(2),'Type','axis');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    xlim([0 15])
    set(h(2),'Units','inches','PaperUnits','inches','PaperSize',[2.5 3],'PaperPosition',[1 1 2.5 3],'Position',[.5 .5 2.5 3])
%     Gram_color = findall(h(1),'Type','colorbar');
%     Gram_axis = findall(h(1),'Type','axes');
%     set(Gram_axis,'Units','inches');
%     set(Gram_color,'Units','inches');
%     cpos=get(Gram_color,'Position');
%     apos=get(Gram_axis,'Position');
%     Gram_fig=get(h(1),'Position');
%     FFT_fig=get(h(2),'Position');
%     FFT_color = findall(h(2),'Type','colorbar');
%     FFT_axis = findall(h(2),'Type','axes');
%     set(FFT_axis,'Units','inches');
%     a2pos=get(FFT_axis,'Position');
%     if ~isempty(FFT_color) && ~isempty(cpos)
%         set(FFT_color,'Units','inches' ,'Position',[cpos(1)-Gram_fig(3)+FFT_fig(3) cpos(2)-Gram_fig(4)+FFT_fig(4) cpos(3)/2 cpos(4)]);
%     end
%     marg=Gram_fig(3)-(apos(1)+apos(3));
%     set(FFT_axis,'Units','inches' ,'Position',[a2pos(1) a2pos(2) FFT_fig(3)-marg-a2pos(1) a2pos(4)]);
    

        mycontrolstruct.optarg='sdf,type,heat,norm,pwelch,[25 81]';
    h(3)=plot_spectral(mycontrolstruct);
    set(h(3),'Name','GammaSpect')
    title('Gamma')
    bf = findall(h(3),'Type','text');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    bf = findall(h(3),'Type','axis');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    xlim([25 80])
    set(h(3),'Units','inches','PaperUnits','inches','PaperSize',[2.5 3],'PaperPosition',[1 1 2.5 3],'Position',[.5 .5 2.5 3])
%     Gram_color = findall(h(1),'Type','colorbar');
%     Gram_axis = findall(h(1),'Type','axes');
%     set(Gram_axis,'Units','inches');
%     set(Gram_color,'Units','inches');
%     cpos=get(Gram_color,'Position');
%     apos=get(Gram_axis,'Position');
%     Gram_fig=get(h(1),'Position');
%     FFT_fig=get(h(3),'Position');
%     FFT_color = findall(h(3),'Type','colorbar');
%     FFT_axis = findall(h(3),'Type','axes');
%     set(FFT_axis,'Units','inches');
%     a2pos=get(FFT_axis,'Position');
%     set(FFT_color,'Units','inches' ,'Position',[cpos(1)-Gram_fig(3)+FFT_fig(3) cpos(2)-Gram_fig(4)+FFT_fig(4) cpos(3)/2 cpos(4)]);
%     marg=Gram_fig(3)-(apos(1)+apos(3));
%     set(FFT_axis,'Units','inches' ,'Position',[a2pos(1) a2pos(2) FFT_fig(3)-marg-a2pos(1) a2pos(4)]);
%     %close(h(2))
% %     
% %     figure(h(1))
% %     ylim([0 200])
% %     set(h(1),'Units','inches','PaperUnits','inches','PaperSize',[5 6],'PaperPosition',[0 0 5 6],'Position',[.5 .5 5 6])
    
        mycontrolstruct.optarg='sdf,type,heat,norm,pwelch,[0 200]';
    h(4)=plot_spectral(mycontrolstruct);
    set(h(4),'Name','FullSpect')
    title('Welch''s Periodogram of SDF')
    bf = findall(h(4),'Type','text');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    bf = findall(h(4),'Type','axis');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    xlim([0 200])
    set(h(4),'Units','inches','PaperUnits','inches','PaperSize',[3.25 3],'PaperPosition',[0 0 3.25 3],'Position',[.5 .5 3.25 3])
    Gram_color = findall(h(1),'Type','colorbar');
    Gram_axis = findall(h(1),'Type','axes');
    set(Gram_axis,'Units','inches');
    set(Gram_color,'Units','inches');
    cpos=get(Gram_color,'Position');
    apos=get(Gram_axis,'Position');
    Gram_fig=get(h(1),'Position');
    FFT_fig=get(h(4),'Position');
    FFT_color = findall(h(4),'Type','colorbar');
    FFT_axis = findall(h(4),'Type','axes');
    set(FFT_axis,'Units','inches');
    a2pos=get(FFT_axis,'Position');
    set(FFT_color,'Units','inches' ,'Position',[cpos(1)-Gram_fig(3)+FFT_fig(3) cpos(2)-Gram_fig(4)+FFT_fig(4) cpos(3) cpos(4)]);
    marg=Gram_fig(3)-(apos(1)+apos(3));
    set(FFT_axis,'Units','inches' ,'Position',[a2pos(1) a2pos(2) FFT_fig(3)-marg-a2pos(1) a2pos(4)]);

if printflag
    for tt=1:length(h)
        printeps(figure(h(tt)),[savepath sl get(h(tt),'Name')])
    end
end
    
myh=tgcross(mycontrolstruct);
if printflag
printeps(figure(myh),[savepath sl get(myh,'Name')])
end
end

    if ismember(20,plotCtrlFigs)
        mycontrolstruct.optarg='pyremph,down';
        
        tmpstruct = mycontrolstruct;
        if isfield(tmpstruct.curses.cells, 'mygids')
            tmpstruct.curses.cells = rmfield(tmpstruct.curses.cells, 'mygids');
        end
        
        if datenum(mrr(1).Date)<735749 || olderflag==1
            tmph=figure('Name','fullraster','Renderer', 'painters','Visible','on','Color','w','Units','inches','PaperUnits','inches','PaperSize',[6 3],'PaperPosition',[0 0 6 3],'Position',[.5 .5 6 3]);
        else
            tmph=figure('Name','fullraster','GraphicsSmoothing','off', 'Renderer', 'painters','Visible','on','Color','w','Units','inches','PaperUnits','inches','PaperSize',[6 3],'PaperPosition',[0 0 6 3],'Position',[.5 .5 6 3]);
        end
        tmph=plot_raster(tmpstruct,gca,colorstruct); % myresultsstruct
        bf = findall(tmph,'Type','text');
        for b=1:length(bf)
            set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
        end
        bf = findall(tmph,'Type','axes');
        for b=1:length(bf)
            set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
        end
        set(tmph,'Units','inches','PaperUnits','inches','PaperSize',[6 3],'PaperPosition',[0 0 6 3],'Position',[.5 .5 6 3])
        if printflag
            printeps(figure(tmph),[savepath sl get(tmph,'Name')])
        end
    end

    
if ismember(14,plotCtrlFigs) 
    olderflag=1;
    handles = mycontrolstruct;
    eval(['global ca1oVgH_ScaledUp_04'])
    evalin('caller',['global ca1oVgH_ScaledUp_04'])
    evalin('base',['global ca1oVgH_ScaledUp_04'])
    handles.curses=[];
    handles.curses.ind=strmatch('ca1oVgH_ScaledUp_04',{RunArray.RunName},'exact');
    if ~isempty(handles.curses.ind)
    guidata(handles.btn_generate, handles)
    h=load_results([],handles);

    myzoomstruct = ca1oVgH_ScaledUp_04;
    myzoomstruct.general=handles.general;
    myzoomstruct.btn_generate = handles.btn_generate;
    myzoomstruct.curses.cells=myhandles.curses.cells;
    if datenum(mrr(1).Date)<735749 || olderflag==1
        hfull=figure('Name','SpkRstrZoom','Renderer', 'painters','Visible','on','Color','w','Units','inches','PaperUnits','inches','PaperSize',[6 3],'PaperPosition',[0 0 6 3],'Position',[.5 .5 6 3]);
    else
        hfull=figure('Name','SpkRstrZoom','GraphicsSmoothing','off', 'Renderer', 'painters','Visible','on','Color','w','Units','inches','PaperUnits','inches','PaperSize',[6 3],'PaperPosition',[0 0 6 3],'Position',[.5 .5 6 3]);
    end
    myzoomstruct.optarg='pyremph,posfil,down';
    hfull=plot_raster(myzoomstruct,gca,colorstruct); % myresultsstruct
    xlim([100 600])

    %set(get(hh(newi),'Children'),'xLim',zoomrange)
    bf = findall(hfull,'Type','text');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    bf = findall(hfull,'Type','axes');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
if printflag
    printeps(figure(hfull),[savepath sl get(hfull,'Name')])
end
    end
end
    
olderflag=1;
% Figure 4: LFP, Raster, Traces
    zoomrange=[1000 2000];
if ismember(4,plotCtrlFigs)

    [hh, colorstruct]=basictrace(mycontrolstruct.runinfo.RunName,mycontrolstruct,zoomrange); % myfigs(1).name,myresultsstruct
    for hz=1:4
        bf = findall(hh(hz),'Type','text');
        for b=1:length(bf)
            set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
        end
        bf = findall(hh(hz),'Type','axes');
        for b=1:length(bf)
            set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
        end
        if printflag
        printeps(figure(hh(hz)),[savepath sl get(hh(hz),'Name')])
        end
    end

    % next
    
    newi=length(hh)+1;
    
    if datenum(mrr(1).Date)<735749 || olderflag==1
        hh(newi)=figure('Name','SpkRstr','Renderer', 'painters','Visible','on','Color','w','Units','inches','PaperUnits','inches','PaperSize',[6 3],'PaperPosition',[0 0 6 3],'Position',[.5 .5 6 3]);
    else
        hh(newi)=figure('Name','SpkRstr','GraphicsSmoothing','off', 'Renderer', 'painters','Visible','on','Color','w','Units','inches','PaperUnits','inches','PaperSize',[6 3],'PaperPosition',[0 0 6 3],'Position',[.5 .5 6 3]);
    end
    mycontrolstruct.optarg='pyremph,posfil,down';
    hh(newi)=plot_raster(mycontrolstruct,gca,colorstruct); % myresultsstruct
    xlim(zoomrange)

    %set(get(hh(newi),'Children'),'xLim',zoomrange)
    bf = findall(hh(newi),'Type','text');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    bf = findall(hh(newi),'Type','axes');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize,'xLim',zoomrange)
    end
if printflag
    printeps(figure(hh(newi)),[savepath sl get(hh(newi),'Name')])
end

    consolflag=1;
    if consolflag
        % find highest position on intratrace
        vv=get(hh(1),'Children');
        for v=1:length(vv)
            set(vv(v),'Units','inches')
            vp=get(vv(v),'Position');
            vppos(v)=vp(2);
            vpht(v)=vp(4);
        end

        [maxp, maxi]=max(vppos);
        maxheight = maxp+vpht(maxi);
        mymarg=.08;
        figp=get(hh(1),'Position');
        maxht=7.5; %figp(4)+6+mymarg;
        set(hh(1),'Position',[0 0 figp(3) maxht],'PaperPosition',[0 0 figp(3) maxht],'PaperSize',[figp(3) maxht]);

        mych=get(hh(5),'Children');
        figrelpos=get(hh(5),'Position');
        for mv=1:length(mych)
            set(mych(mv),'Units','inches')
            bp=get(mych(mv),'Position');
            set(mych(mv),'Parent',hh(1),'Position',[bp(1) bp(2)+maxheight+mymarg bp(3) bp(4)])
        end
        close(hh(5));
        bp=get(mych(1),'Position');
        newmaxheight = bp(2)+mymarg + bp(4);

        lfpg=get(hh(4),'Children');
        lfpt=get(hh(3),'Children');
    
        set(lfpg(2),'Units','inches')
        mp=get(lfpg(2),'Position');
        set(lfpg(2),'Parent',hh(1),'Position',[mp(1) mp(2)+newmaxheight mp(3) mp(4)])

        newmaxheight = mp(2)+newmaxheight + mp(4);
        set(lfpt(2),'Units','inches')
        mp=get(lfpt(2),'Position');
        set(lfpt(2),'Parent',hh(1),'Position',[mp(1) mp(2)+newmaxheight mp(3) mp(4)])
        close(hh(3));
        close(hh(4));
        
        for hzz=1
            bf = findall(hh(hzz),'Type','text');
            for b=1:length(bf)
                set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
            end
            bf = findall(hh(hzz),'Type','axes');
            for b=1:length(bf)
                set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
            end
            if printflag

            printeps(figure(hh(hzz)),[savepath sl get(hh(hzz),'Name') 'New'])
            end
        end
    end
end   

if ismember(5,plotCtrlFigs)
% Figure 5: Firing rates and phases

    tmporg=mycontrolstruct;
    tmporg.curses.cells=mycontrolstruct.curses.cells([7 8 3 9 1 2 6 4 5]);
    tmporg.optarg='subset';
    [hhh(1), firedata] = plot_firingrates(tmporg);
    NiceAbbrev = {'Pyr','PV+B','CCK+B','SC-A','Axo','Bis','O-LM','Ivy','NGF'};

    set(hhh(1),'color','w','Name','Active Cell Firing Rates','units','inches','PaperUnits','inches','PaperSize',[3 2],'Position',[.5 .5 3 2],'PaperPosition',[0 0 3 2]);
    formatter(gca)
    ty=ylabel('Firing frequency (Hz)');
    formatter(ty)
    set(gca,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
    set(ty,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
    legstr={'Model','Exp. Anesth.','Exp. Awake'};
    bfm=legend(legstr,'Position',[0.5705    0.6719    0.4387    0.3003]);
    legend boxoff
    formatter(bfm);
    set(bfm,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
    title('')
    bf = findall(hhh(1),'Type','text');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    bf = findall(hhh(1),'Type','axis');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    set(gca,'xticklabel','','xlim',[.5 9.25])
    for ar=1:length(NiceAbbrev)
        myr(ar)=text(ar,-3,NiceAbbrev{ar},'Rotation',30,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize,'HorizontalAlignment','right');
    end
    set(gca,'Position',[.13 .15 .85 .83])
if printflag

    printeps(figure(hhh(1)),[savepath sl get(hhh(1),'Name')])
end

    hexp=printExpData(mycontrolstruct);
    close(hexp(2));
    hexp(2)=[];

    mycontrolstruct.optarg='hist';
    mdx=find([FinalData.Idx]==mycontrolstruct.curses.ind);

    if isfield(FinalData(mdx),'mydata') && ~isempty(FinalData(mdx).mydata)
        [tmp, mydata, peaks] = plot_phasemod([],mycontrolstruct,1,FinalData(mdx).mydata);
    else
        [tmp, mydata, peaks] = plot_phasemod([],mycontrolstruct,1);
        FinalData(mdx).mydata = mydata;
    end

    myHz=8; % Hz
    thetaper=1000/myHz;
    for r=1:length(mydata)
        [angle, magnitude, modtimes]=getspikephase(mydata(r).activitytimes, thetaper, peaks);
        mycellSD(r).name=mydata(r).tech;
        mycellSD(r).sptimes=mydata(r).activitytimes;
        mycellSD(r).spphases=modtimes/thetaper*360;
        
        fid=fopen([savepath sl 'Spike_Phase_Time_' mydata(r).tech '.txt'],'w');
        for q=1:length(mycellSD(r).spphases)
            fprintf(fid,'%f\t%.1f\n',mycellSD(r).sptimes(q),mycellSD(r).spphases(q));
        end
        fclose(fid);
    end

    myfg=figure('Color','w','Name','PhasexMod','units','inches','PaperUnits','inches','Position',[0 0 3 2.5],'PaperPosition',[0 0 3 2.5],'PaperSize',[3 2.5]);
    tmpphase=[];
    for m=1:9
        tmpphase(m)=mod(mydata(m).phase-90,360)+90;
        plot(mydata(m).phase,mydata(m).mod,'Color',mydata(m).color,'Marker','o','MarkerSize',5,'LineWidth',2,'LineStyle','none')
        if mydata(m).mod>.5
            btm=text(mydata(m).phase,mydata(m).mod,['    ' mydata(m).name],'Color',mydata(m).color,'HorizontalAlignment','right');
        else
            btm=text(mydata(m).phase,mydata(m).mod,[' ' mydata(m).name],'Color',mydata(m).color,'HorizontalAlignment','left');
        end
        if mod(length(find(mydata(m).mod>[mydata(1:9).mod])),2)==0
            set(btm,'VerticalAlignment','bottom')
        else
            set(btm,'VerticalAlignment','top')
        end
        hold on
    end
    p = polyfit(tmpphase,[mydata(1:9).mod],1);
    yfit = polyval(p,tmpphase);
    yresid = [mydata(1:9).mod] - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length([mydata(1:9).mod])-1) * var([mydata(1:9).mod]);
    rsq = 1 - SSresid/SStotal;
    plot(tmpphase,yfit,'k','LineWidth',2)
    text(mean(tmpphase),mean(yfit),['    r^2: ' sprintf('%.2f',rsq)])
    ylim([-.1 1])
    xlim([0 360])
    xlabel('Theta Phase Preference (^o)')
    ylabel('Theta Modulation of Spiking')

    bf = findall(myfg,'Type','text');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    bf = findall(myfg,'Type','axis');
    for b=1:length(bf)
        set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    end
    set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
    set(gca,'xtick',[0 90 180 270 360])
    box off
    if printflag

    printeps(figure(myfg),[savepath sl get(myfg,'Name')])
    end
    
    for r=1:length(FinalData(mdx).mydata)
        mycontrolstruct.curses.cells(r).phase=FinalData(mdx).mydata(r).phase;
        mycontrolstruct.curses.cells(r).color=FinalData(mdx).mydata(r).color;
        mycontrolstruct.curses.cells(r).offset=[0 0];
    end
    
    swap1=tmp;
    tmp=hexp(1);
    hexp(1)=swap1;
    
    bc=get(hexp(1),'Children');
    for b=1:length(bc)
        set(bc(b),'Units','inches')
    end
    set(hexp(1),'Units','inches','PaperUnits','inches')
    pos=get(hexp(1),'Position');
    set(hexp(1),'Position',[pos(1) pos(2) pos(3)*2 pos(4)],'PaperPosition',[pos(1) pos(2) pos(3)*2 pos(4)],'PaperSize',[pos(3)*2 pos(4)])
    movew=pos(3);
    bd=get(tmp,'Children');
    for b=1:length(bd)
        set(bd(b),'Units','inches')
        pos=get(bd(b),'Position');
        set(bd(b),'Position',[pos(1)+movew pos(2) pos(3) pos(4)])
        set(bd(b),'Parent',hexp(1))
    end
    close(tmp)

    hexp(2)=MiniPhaseShift_Figure(mycontrolstruct.curses.cells(1:9),FinalData(mdx).Pos.Freqs(FinalData(mdx).Pyr).Theta,[]);
    set(hexp(2),'Units','inches','Position',[3 1 3 2],'Units','inches','PaperPosition',[0 0 3 2],'PaperSize',[3 2])

%     for g=1:2
%         bc=get(hexp(g),'Children');
%         for b=1:length(bc)
%             set(bc(b),'Units','normalized')
%             rr=get(bc(b),'Children');
%             for r=1:length(rr)
%                 if isprop(rr(r),'Units')
%                     set(rr(r),'Units','normalized')
%                 end
%             end
%         end
%         pos=get(hexp(g),'Position');
%         sc=4/pos(3);
%         set(hexp(g),'Position',[pos(1) pos(2) pos(3)*sc pos(4)*sc*.9])
%     end

        for h=1:2
            bf = findall(hexp(h),'Type','text');
            for b=1:length(bf)
                set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
            end
            bf = findall(hexp(h),'Type','axes');
            for b=1:length(bf)
                set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
            end
            if printflag

            printeps(figure(hexp(h)),[savepath sl get(hexp(h),'Name')])
            end
        end
end

function GABAFig(FinalData,handles)
global myFontSize myFontWeight myFontName savepath printflag printtable fid sl
% This needs additional datasets so I sequestered this code for now unless
% you want to run it too, just let me know if you do.

% GABA
myh=1.1*1.5;
hsc=1.5;
hg=figure('GraphicsSmoothing','off', 'Renderer', 'painters','color','w','Name','GABA','units','inches','PaperUnits','inches','PaperSize',[2.5 myh],'Position',[.5 .5 2.5 myh],'PaperPosition',[0 0 2.5 myh]);
conds={'ctrl','a','w'};
xL={'Control','No GABA_B','Equiv. CT'};
xLline={'Control','No GABA$_B$','No GABA$_B$, Equiv. CT'};

if printtable
    fprintf(fid,'\\multicolumn{7}{|l|}{\\textbf{%s}} \\\\ \n', 'GABA$_{B}$');
    fprintf(fid,'\\hline\n');
end

usenorm=0;

mdx=1;

allgabadata=[];
grouping=[];
totspks=[];
totvar=[];
for x=1:length(conds)
    xa(x)=x;
    Runs2Use=[];
    
    for r=1:length(FinalData)
        if ~isempty(strfind(FinalData(r).Comments,['gaba ' conds{x}]))
            Runs2Use=[Runs2Use r];
        end
    end
    if strcmp(conds{x},'ctrl') % the control run will not have 'ctrl' in its comments, will instead have 'control', so add that in manually.
        Runs2Use=[Runs2Use mdx];
    end
        
    thetafreq=[];
    thetapow=[];
    thetanorm=[];
    gammafreq=[];
    gammapow=[];
    overallfreq=[];
    tmpspks=[];
    spkvar=[];
    
    for bb=1:length(Runs2Use)
        thetafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Theta;
        thetapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP;
        thetanorm(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).ThetaP/FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
        gammafreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Gamma;
        gammapow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).GammaP;
        overallfreq(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).Overall;
        overallpow(bb)=FinalData(Runs2Use(bb)).Pos.Freqs(FinalData(Runs2Use(bb)).Pyr).OverallP;
                
        handles.curses=[];
        handles.curses.ind=FinalData(Runs2Use(bb)).Idx;
        if exist(FinalData(Runs2Use(bb)).Name,'var')
            eval(['handles=' FinalData(Runs2Use(bb)).Name ';']);
        else
            guidata(handles.btn_generate, handles)

            spikeraster(handles.btn_generate,guidata(handles.btn_generate))
            handles=guidata(handles.btn_generate);

            getcelltypes(handles.btn_generate,guidata(handles.btn_generate))
            handles=guidata(handles.btn_generate);
        end

        mdidx=find(handles.curses.spikerast(:,2)>=handles.curses.cells(7).range_st & handles.curses.spikerast(:,2)<=handles.curses.cells(7).range_en);
        tmpspks(bb)=length(mdidx);
        handles.curses.spikerast = handles.curses.spikerast(mdidx,:);
        
        handles.runinfo.SimDuration = FinalData(Runs2Use(bb)).SimDuration;
        
        tmphandles=getoscphaseamp(handles,'sdf',0);
        pyrphase=tmphandles.curses.oscspikerast(:,end-2);        % tmphandles.curses.oscspikerast(:,3)==6
        spkvar(bb) = std(pyrphase);
    end
    
    if usenorm==1
        ya(x)=mean(thetanorm);
        ea(x)=std(thetanorm);
    else
        ya(x)=mean(thetapow);
        ea(x)=std(thetapow);
    end
    fa(x)=mean(thetafreq);
    ga(x)=mean(gammafreq);
    gp(x)=mean(gammapow);
    fo(x)=mean(overallfreq);
    fp(x)=mean(overallpow);
    
    if printtable
        fprintf(fid,'%s & %.1f & %.1e & %.1f & %.1e & %.1f & %.1e \\\\ \n', xLline{x}, mean(thetafreq), mean(thetapow), mean(gammafreq), mean(gammapow), mean(overallfreq), mean(overallpow));
        fprintf(fid,'\\hline\n');  
    end
    
    
    
    % for x=1:length(conds)
    %eval(['newmag = max(synfigdata_' corsyndata{x} '.fig.axis.data(11).y)']);
    if usenorm==1
        allgabadata = [allgabadata; thetanorm(:)];
        %IPSCmag = [IPSCmag; newmag*ones(length(thetanorm),1)];
        grouping = [grouping; x*ones(length(thetanorm),1)];
    else
        allgabadata = [allgabadata; thetapow(:)];
        %IPSCmag = [IPSCmag; newmag*ones(length(thetapow),1)];
        grouping = [grouping; x*ones(length(thetapow),1)];
    end
    totspks = [totspks; tmpspks(:)];
    totvar = [totvar; spkvar(:)];
end
gabacolor={'k','r','g'};

%[h,atab,ctab,stats] = aoctool(IPSCmag,allgabadata,grouping);
% OR
try
P = anova1(allgabadata,grouping);

% OR
[hp,atab,ctab,stats] = aoctool(totvar,allgabadata,grouping);

end

for r=1:3:9
    p = polyfit(totvar(r:r+2),allgabadata(r:r+2),1);
    yfit = polyval(p,totvar(r:r+2));
    yresid = allgabadata(r:r+2) - yfit;
    SSresid = sum(yresid.^2);
    SStotal = (length(allgabadata(r:r+2))-1) * var(allgabadata(r:r+2));
    rsq = 1 - SSresid/SStotal;
%    disp(['r^2 = ' num2str(rsq) ' and rho = ' num2str(corr(totvar(r:r+2),allgabadata(r:r+2)))])
end


for i=1:length(ya)
    h(i) = bar(i, ya(i));
    if i == 1, hold on, end
    set(h(i), 'FaceColor', gabacolor{i},'EdgeColor','none','BarWidth',.4) % get(h(i),'BarWidth')*.7
    errorbar(i,ya(i),ea(i),'Color',gabacolor{i},'Marker','none','MarkerSize',15,'LineStyle','none','LineWidth',1.5)
end
%bar(xa,ya,'EdgeColor',[0 0 0],'FaceColor',[0 0 0],'BarWidth',.4)
set(gca,'Position',[.2 .18 .77 .77],'LineWidth',1.5)
hold on

try
sprintf('%.6f\n',P)
if P<.001
    H = sigstar({[1,2],[2,3],[3,1]});
end
end
%errorbar(xa,ya,ea,'k','Marker','none','MarkerSize',15,'LineStyle','none','LineWidth',1.5)
set(gca,'XTickLabel',{});
ypos = -max(ylim)/45;% -max(ylim)/10;
N=length(xL);
gw=text(1:N,repmat(ypos,N,1),xL','verticalalignment','top','horizontalalignment','center','Rotation',0,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName);
formatter(gw)
set(gw,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
ylabel('Theta Power','FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
    for h=1:1
        bf = findall(hg(h),'Type','text');
        for b=1:length(bf)
            set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
        end
        bf = findall(hg(h),'Type','axes');
        for b=1:length(bf)
            set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
        end
    end
    
xlim([0.5 3.5])
bb=get(gca,'Position');
set(gca,'Position',[bb(1) bb(2) bb(3) .72])
if printflag
printeps(figure(gcf),[savepath sl get(gcf,'Name')])
end

mycolors=[];
mycolors(1,:)=[0 0 0];
mycolors(2,:)=[1 0 0];
mycolors(3,:)=[0 1 0];
yL='Theta Power';

myfreqs=[];
freqstruct=[];
for r=1:length(fa)
    myfreqs(r,:)=[fa(r) ya(r) fo(r) fp(r) ga(r) gp(r)];
    freqstruct(r).freqs= ga(r);%unique([FinalData(mdx).tbldata{:,6}]);
    freqstruct(r).pows= gp(r);%unique([FinalData(mdx).tbldata{:,6}]);
end
mysize=[3 1.1];
mytitle='';
myapp='';
myplotvals=ya;
myplotvals(isnan(myplotvals))=0;

gy=plotmebaronly('GABAFreq', mycolors, myplotvals,mysize,xL,yL,myfreqs,freqstruct,mytitle,myapp);
if printflag
printeps(gy(2),[savepath sl 'GABAFreq'])
end

function gy=plotmebaronly(myname, mycolors, myplotvals, mysize,xL,yL,varargin)
global myFontSize myFontName myFontWeight printtable fid

alg=.5;
sc=.37;
hsc=1.5;
clipin=.5;


if isempty(varargin)
    myfreqs=[];
    freqstruct=[];
    mytitle='';
    myapp='';
    gy(1)=figure('Color','w','Name',myname,'Units','inches','PaperUnits','inches','PaperSize',[alg+mysize(1)*sc mysize(2)],'Position',[.5 .5 alg+mysize(1)*sc mysize(2)],'PaperPosition',[0 0 alg+mysize(1)*sc mysize(2)]);
    gy(2)=figure('Color','w','Name',[myname 'Freq'],'Units','inches','PaperUnits','inches','PaperSize',[alg+mysize(1)*sc mysize(2)],'Position',[.5 .5 alg+mysize(1)*sc mysize(2)],'PaperPosition',[0 0 alg+mysize(1)*sc mysize(2)]);
else
    myfreqs=varargin{1};
    freqstruct=varargin{2};
    mytitle=varargin{3};
    myapp=varargin{4};
    gy(1)=figure('Color','w','Name',myname,'Units','inches','PaperUnits','inches','PaperSize',[alg+mysize(1)*sc mysize(2)*hsc],'Position',[.5 .5 alg+mysize(1)*sc mysize(2)*hsc],'PaperPosition',[0 0 alg+mysize(1)*sc mysize(2)*hsc]);
    gy(2)=figure('Color','w','Name',[myname 'Freq'],'Units','inches','PaperUnits','inches','PaperSize',[alg+mysize(1)*sc mysize(2)*hsc],'Position',[.5 .5 alg+mysize(1)*sc mysize(2)*hsc],'PaperPosition',[0 0 alg+mysize(1)*sc mysize(2)*hsc]);
    if 1==0 % don't sort again
        [~, sorti]=sort(myfreqs(:,3));
        sorti=[1; sorti(sorti~=1)];
        mycolors=mycolors(sorti,:);
        myplotvals=myplotvals(sorti);
        myfreqs=myfreqs(sorti,:);
        xL=xL(sorti);
        freqstruct=freqstruct(sorti);
    end
end
if ~isempty(myapp)
    for x=1:length(xL)
        xL{x}={xL{x},myapp};
    end
end

if printtable
fprintf(fid,'\\multicolumn{7}{|l|}{\\textbf{%s}} \\\\ \n', strrep(strrep(strrep(mytitle,'&','\&'),'#','\#'),'GABA_B','GABA$_B$'));
fprintf(fid,'\\hline\n');
end

N = numel(myplotvals);
% if isempty(myfreqs)
%     bargraph=subplot('Position',[alg/(alg+sc*mysize(1)) .14 1-1.1*alg/(alg+sc*mysize(1)) .8]);
% else
figure(gy(1))
    bargraph=subplot('Position',[alg/(alg+sc*mysize(1)) 0.1120+.1 1-1.1*alg/(alg+sc*mysize(1)) .7]);
% end
for i=1:N
    if printtable
        fprintf(fid,'%s & %.1f & %.1e & %.1f & %.1e & %.1f & %.1e \\\\ \n', strrep(strrep(strrep(xL{i},'&','\&'),'#','\#'),'GABA_B','GABA$_B$'), myfreqs(i,1), myfreqs(i,2), myfreqs(i,5), myfreqs(i,6), myfreqs(i,3), myfreqs(i,4));
        fprintf(fid,'\\hline\n'); 
    end
    h(i) = bar(i, myplotvals(i));
    if i == 1, hold on, end
    set(h(i), 'FaceColor', mycolors(i,:),'EdgeColor','none','BarWidth',.8) % get(h(i),'BarWidth')*.7
end   
set(gca, 'XTickLabel', '') 
yy=ylim;
set(gca,'ylim',[-yy(2)*.02 yy(2)])
for i=1:N
    h(i).BaseLine.BaseValue=-yy(2)*.02;
end
% mrr=get(gca,'YTickLabel');
% %mystr=length(strfind(num2str(mrr{2}),'0'));
% for m=2:length(mrr)
%     mystr=length(strfind(num2str(mrr{m}),'0'));
%     mrr{m}=[mrr{m}(1:end-mystr) 'e' num2str(mystr)];
%     if strcmp(mrr{m},'0')
%         m
%     end
% end
%set(gca,'YTickLabel',mrr)
set(gca,'YTickLabelMode','auto')
btmp2=xlabel(mytitle);
%set(btmp2,'Position',get(btmp2,'Position')- [0 diff(get(gca,'ylim'))/20 0])
set(btmp2,'Position',get(btmp2,'Position')- [0 diff(get(gca,'ylim'))/10 0])
mm=get(gca,'xlabel');
 mp=get(mm,'Position');
% set(mm,'Position',[mp(1) mp(2)-.2*mp(2) mp(3)])
ypos = -diff(ylim)/13;% -max(ylim)/10;
rotflag=0;
        hA='center';
for x=1:length(xL)
    if length(xL{x})>4
        ypos = -diff(ylim)/16;% -max(ylim)/10;
        rotflag=25;
        hA='right';
        set(mm,'Position',[mp(1) mp(2)+.6*mp(2) mp(3)])
        set(bargraph,'Position',[alg/(alg+sc*mysize(1)) 0.1120+.04+.1 1-1.1*alg/(alg+sc*mysize(1)) .7-.04]);
    end
end
gw=text(1:N,repmat(ypos,N,1),xL','horizontalalignment',hA,'Rotation',rotflag,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName);
formatter(gw)
set(gw,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
btmp=ylabel(yL);
%btmp2 = get(gca,'XLabel');

formatter(btmp)
formatter(btmp2)
set(btmp,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
set(btmp2,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
formatter(gca)
try
    set(get(gy,'Children'),'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
end

bf = findall(gy,'Type','text');
for b=1:length(bf)
    set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
end
bf = findall(gy,'Type','axes');
for b=1:length(bf)
    set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
end

figure(gy(2))    
thetagraph=subplot('Position',[alg/(alg+sc*mysize(1)) 0.1120+.2 1-1.1*alg/(alg+sc*mysize(1)) .6]);
%thetagraph=subplot('Position',[alg/(alg+sc*mysize(1)) .43 1-1.1*alg/(alg+sc*mysize(1)) .44]);
thetapos=get(thetagraph,'Position');
%set(thetagraph,'XLim',get(bargraph,'XLim'),'YLim',[0 80])
set(thetagraph,'XLim',get(bargraph,'XLim'),'YLim',[0 25])
th=.02*diff(get(thetagraph,'YLim'))*(.2/thetapos(4));
formatter(ylabel('Modul. Freq. of SDF (Hz)'),'VerticalAlignment','top')
%formatter(ylabel({'Frequency of','Modulation of SDF (Hz)'}),'VerticalAlignment','top')
set(get(gca,'ylabel'),'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
%formatter(ylabel('         Theta (Hz)'),'VerticalAlignment','top')
xrange=get(bargraph,'XLim');
patch([xrange(1) xrange(2) xrange(2) xrange(1) xrange(1)],[5 5 10 10 5],[.94 .94 .94],'EdgeColor','none')
hold on
%     patch([xrange(1) xrange(2) xrange(2) xrange(1) xrange(1)],[25 25 80 80 25],[.94 .94 .94],'EdgeColor','none')
for i=1:N
    ht=myfreqs(i,3);
    %ht=myfreqs(i,1);
    st=get(h(i),'XData')-get(h(i),'BarWidth')/2;
    en=get(h(i),'XData')+get(h(i),'BarWidth')/2;
    if ht>25
        gty=text(st,25,num2str(round(ht*10)/10));
        set(gty,'Color',mycolors(i,:),'HorizontalAlignment','left','VerticalAlignment','top','Interpreter','none')
    else
        plot([st en],[ht ht],'Color',mycolors(i,:))
    end
    hold on
%         for r=1:length(freqstruct(i).freqs)
%             hold on
%             plot([get(h(i),'XData')-get(h(i),'BarWidth')/2 get(h(i),'XData')+get(h(i),'BarWidth')/2],[freqstruct(i).freqs(r) freqstruct(i).freqs(r)],'Color',mycolors(i,:),'LineStyle','-.') %mycolors(i,:)
%         end
end

set(thetagraph, 'XTickLabel', '','XTick',[]) 

formatter(thetagraph)
% try
% set(thetagraph,'XColor','none');
% end
set(thetagraph,'Layer','top');
set(thetagraph,'Clipping','Off','FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)


    btmp2=xlabel(mytitle);
%set(btmp2,'Position',get(btmp2,'Position')- [0 diff(get(gca,'ylim'))/20 0])
set(btmp2,'Position',get(btmp2,'Position')- [0 diff(get(gca,'ylim'))/10 0])
mm=get(gca,'xlabel');
 mp=get(mm,'Position');
% set(mm,'Position',[mp(1) mp(2)-.2*mp(2) mp(3)])
ypos = -diff(ylim)/13;% -max(ylim)/10;
rotflag=0;
        hA='center';
for x=1:length(xL)
    if length(xL{x})>4
        ypos = -diff(ylim)/16;% -max(ylim)/10;
        rotflag=25;
        hA='right';
        set(mm,'Position',[mp(1) mp(2)+.6*mp(2) mp(3)])
        set(thetagraph,'Position',[alg/(alg+sc*mysize(1)) 0.1120+.04+.1 1-1.1*alg/(alg+sc*mysize(1)) .7-.04]);
    end
end
gw=text(1:N,repmat(ypos,N,1),xL','horizontalalignment',hA,'Rotation',rotflag,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName);
formatter(gw)
set(gw,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
btmp=ylabel('Frequency (Hz)');
%btmp2 = get(gca,'XLabel');
%xlabel('')
formatter(btmp)
formatter(btmp2)
set(btmp,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
set(btmp2,'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
formatter(gca)
try
    set(get(gy,'Children'),'FontSize',myFontSize,'FontWeight',myFontWeight,'FontName',myFontName)
end

bf = findall(gy,'Type','text');
for b=1:length(bf)
    set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
end
bf = findall(gy,'Type','axes');
for b=1:length(bf)
    set(bf(b),'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
end


function handles = CustomNearElectrode(handles,cellnumbers,LayerHeights,LongitudinalLength,TransverseLength,ElectrodeLoc,maxdist)

layers=regexp(LayerHeights,';','split');
LayerVec=str2double(layers(2:end-1));
layind=cellnumbers.data(:,2)+1;
ZHeight=zeros(1,length(handles.curses.cells));
for r=1:length(handles.curses.cells)
    BinInfo(r) = setBins(handles.curses.cells(r).numcells,LongitudinalLength,TransverseLength,LayerVec(layind(r)));
    ZHeight(r) = sum(LayerVec(1:layind(r)-1));
end

for postype=1:length(handles.curses.cells)
    handles.curses.cells(postype).mygids=[];
    for gid=handles.curses.cells(postype).range_st:handles.curses.cells(postype).range_en
        pos = getpos(gid, handles.curses.cells(postype).range_st, BinInfo(postype), ZHeight(postype));
        mydist=sqrt((pos.x - ElectrodeLoc(1)).^2 + (pos.y - ElectrodeLoc(2)).^2);
        if mydist<maxdist
            handles.curses.cells(postype).mygids=[handles.curses.cells(postype).mygids gid];
        end
    end
    disp(['Just finished cell type r=' num2str(postype)])
end

function [h, mytrace]=basictrace(resultname,handles,zoomrange)
global myFontSize sl repodir myFontWeight myFontName savepath

totrange=0;
myw=6;
olderflag=1;

posflag=0;
if isfield(handles.curses.cells(1),'mygids')
    posflag=1;
end

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'pyramidalcell');
end
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_pyramidalcell99180.dat']);
    disp('using a somatic recording that may be outside of the electrode range')
end
mytrace.pyramidal.trace=zztrace.data;
mytrace.pyramidal.color=[.0 .0 .6];
mytrace.pyramidal.pos=1;
mytrace.pyramidal.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.pyramidal.range;

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'pvbasketcell');
end
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_pvbasketcell333500.dat']);
end
mytrace.pvbasket.trace=zztrace.data;
mytrace.pvbasket.color=[.0 .75 .65];
mytrace.pvbasket.pos=2;
mytrace.pvbasket.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.pvbasket.range;

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'cckcell');
end
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_cckcell5840.dat']);
end
mytrace.cck.trace=zztrace.data;
mytrace.cck.color=[1 .75 .0];
mytrace.cck.pos=3;
mytrace.cck.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.cck.range;

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'scacell');
end
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_scacell338400.dat']);
end
mytrace.sca.trace=zztrace.data;
mytrace.sca.color=[1 .5 .3];
mytrace.sca.pos=4;
mytrace.sca.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.sca.range;

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'axoaxoniccell');
end
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_axoaxoniccell360.dat']);
end
mytrace.axoaxonic.trace=zztrace.data;
mytrace.axoaxonic.color=[1 .0 .0];
mytrace.axoaxonic.pos=5;
mytrace.axoaxonic.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.axoaxonic.range;

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'bistratifiedcell');
end
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_bistratifiedcell1800.dat']);
end
mytrace.bistratified.trace=zztrace.data;
mytrace.bistratified.color=[.6 .4 .1];
mytrace.bistratified.pos=6;
mytrace.bistratified.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.bistratified.range;

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'olmcell');
end
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_olmcell20900.dat']); % 20818
end
mytrace.olm.trace=zztrace.data;
mytrace.olm.color=[.5 .0 .6];
mytrace.olm.pos=7;
mytrace.olm.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.olm.range;

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'ivycell');
end
zztrace=[];
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_ivycell10360.dat']); % 13000
end
mytrace.ivy.trace=zztrace.data;
mytrace.ivy.color=[.6 .6 .6];
mytrace.ivy.pos=8;
mytrace.ivy.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.ivy.range;

zztrace=[];
if posflag==1
    zztrace=loadmyfile(handles,resultname,'ngfcell');
end
if isempty(zztrace)
    zztrace=importdata([repodir sl 'results' sl resultname sl 'trace_ngfcell17870.dat']);
end
mytrace.ngf.trace=zztrace.data;
mytrace.ngf.color=[1 .1 1];
mytrace.ngf.pos=9;
mytrace.ngf.range=max(zztrace.data(:,2))-min(zztrace.data(:,2))+2;
totrange = totrange + mytrace.ngf.range;

matchstr={'pyramidalcell','olmcell','bistratifiedcell','axoaxoniccell','pvbasketcell','cckcell','scacell','ivycell','ngfcell'};
NiceAbbrev = {'Pyr.','O-LM','Bis.','Axo.','PV+ B.','CCK+ B.','S.C.-A.','Ivy','NGF.'};
handles.formatP.left = .065;
handles.formatP.bottom=.065;

intracehgt=4;
if olderflag==1
    h=figure('Renderer', 'painters','Visible','on','Color','w','Units','inches','PaperUnits','inches','PaperSize',[myw intracehgt],'PaperPosition',[0 0  myw intracehgt],'Position',[.5 .5  myw intracehgt],'Name','IntracellularTrace');
else
    h=figure('GraphicsSmoothing','off', 'Renderer', 'painters','Visible','on','Color','w','Units','inches','PaperUnits','inches','PaperSize',[myw intracehgt],'PaperPosition',[0 0  myw intracehgt],'Position',[.5 .5  myw intracehgt],'Name','Trace');
end
%pos=get(gcf,'Units');
%set(gcf,'Units','normalized','Position',[0.1 0.1 .9 .9]);
%set(gcf,'Units',pos);

prevend=.96;%+handles.formatP.bottom*2;
mycells=fieldnames(mytrace);
for m=1:length(mycells)
    %subplot(length(mycells),1,mytrace.(mycells{m}).pos)
    tidx=find(mytrace.(mycells{m}).trace(:,1)>=zoomrange(1) & mytrace.(mycells{m}).trace(:,1)<=zoomrange(2));
    g(m)=subplot('Position',[handles.formatP.left*2 prevend-mytrace.(mycells{m}).range/totrange*(.96-handles.formatP.bottom*2) .9-handles.formatP.left .95*mytrace.(mycells{m}).range/totrange*(.96-handles.formatP.bottom*2)]);
    prevend=prevend-mytrace.(mycells{m}).range/totrange*(.96-handles.formatP.bottom*2);
    plot(g(m),mytrace.(mycells{m}).trace(tidx,1),mytrace.(mycells{m}).trace(tidx,2),'Color',mytrace.(mycells{m}).color,'LineWidth',1.0)
    if mytrace.(mycells{m}).pos==9
        ylim([min(mytrace.(mycells{m}).trace(tidx,2))-5 max(mytrace.(mycells{m}).trace(tidx,2))+1])
    else
        ylim([min(mytrace.(mycells{m}).trace(tidx,2))-1*3 max(mytrace.(mycells{m}).trace(tidx,2))+1])
    end
    z=strmatch([mycells{m} 'cell'],matchstr,'exact');
    bd=ylabel(NiceAbbrev{z},'Rotation',0,'HorizontalAlignment','Right','FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize);
    set(bd, 'Units', 'Normalized', 'Position', [-0.03, 0.3, 0]);
    box off
    set(g(m),'Clipping','off','ycolor',mytrace.(mycells{m}).color,'xcolor','w','YTick',[],'YTickLabel',{},'LineWidth',0.5)
%     if mytrace.(mycells{m}).pos==9
%         xlabel('Time (ms)','FontName','ArialMT','FontWeight','bold','FontSize',myFontSize)
%         set(g(m),'FontName','ArialMT','FontWeight','bold','FontSize',myFontSize,'xcolor','w')
%     else
        set(g(m),'XTick',[],'XTickLabel',{})
%     end
    ax2(m) = axes('Position', get(g(m), 'Position'),'Color','none');
    set(ax2(m),'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top','LineWidth',1.1)
end
linkaxes(g,'x');

axes(g(end))
hold on
yL=get(gca,'YLim');
xL=get(gca,'XLim');
plot([xL(2)-3 xL(2)-3],[yL(1)+2 yL(1)+52],'k','LineWidth',2)
plot([xL(2)-103 xL(2)-3],[yL(1)+2 yL(1)+2],'k','LineWidth',2)
%text(xL(2)-30,yL(2),{'100 ms','20 mV'},'HorizontalAlignment','right','FontName','ArialMT','FontSize',myFontSize)
pospos=get(gca,'Position');
ngfdiff=(50/diff(ylim))*pospos(4);

if olderflag==1
    h(2)=figure('Renderer', 'painters','Visible','on','Color','w','PaperUnits','inches','PaperSize',[myw 1.2],'PaperPosition',[.3 .3 myw 1.2],'Name','LFP','Units','inches','Position',[.5 .5 myw 1.2]);
else
    h(2)=figure('GraphicsSmoothing','off', 'Renderer', 'painters','Visible','on','Color','w','PaperUnits','inches','PaperSize',[myw 1.2],'PaperPosition',[.3 .3 myw 1.2],'Name','LFP','Units','inches','Position',[.5 .5 myw 1.2]);
end
%pos=get(gcf,'Units');
%set(gcf,'Units','normalized','Position',[0.1 0.1 .9 .3]);
%set(gcf,'Units',pos);
subplot('Position',[handles.formatP.left*2 handles.formatP.bottom*2 .9-handles.formatP.left .95*(1-handles.formatP.bottom*2)]);
tidx=find(handles.curses.lfp(:,1)>=zoomrange(1) & handles.curses.lfp(:,1)<=zoomrange(2));


% if posflag==1
%     plot(handles.curses.lfp(tidx,1),handles.curses.epos.lfp(tidx),'Color','k','LineWidth',2)
% else
    plot(handles.curses.lfp(tidx,1),handles.curses.lfp(tidx,2),'Color','k','LineWidth',2)
% end
bd=ylabel('LFP','Rotation',0,'HorizontalAlignment','Right','FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize);
set(bd, 'Units', 'Normalized', 'Position', [-0.03, 0.3, 0]);
box off
xlabel('Time (ms)','FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
set(gca,'Clipping','off','ycolor','k','xcolor','w','YTick',[],'YTickLabel',{})
ax2 = axes('Position', get(gca, 'Position'),'Color','none');
set(ax2,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','on','layer','top')

lfpht=.7; %.55;

filteredlfp=mikkofilter(handles.curses.lfp,1000/handles.runinfo.lfp_dt);
thetalfp=filteredlfp;
if olderflag==1
    h(3)=figure('Renderer', 'painters','Visible','on','Color','w','PaperUnits','inches','PaperSize',[myw lfpht],'PaperPosition',[0 0 myw lfpht],'Name','Filtered LFP','Units','inches','Position',[.5 .5 myw 1]);
else
    h(3)=figure('GraphicsSmoothing','off', 'Renderer', 'painters','Visible','on','Color','w','PaperUnits','inches','PaperSize',[myw lfpht],'PaperPosition',[0 0 myw lfpht],'Name','Filtered LFP','Units','inches','Position',[.5 .5 myw lfpht]);
end
subplot('Position',[handles.formatP.left*2 handles.formatP.bottom .9-handles.formatP.left .95*(1-handles.formatP.bottom*2)]);
tidx=find(filteredlfp(:,1)>=zoomrange(1) & filteredlfp(:,1)<=zoomrange(2));
plot(filteredlfp(tidx,1),filteredlfp(tidx,2),'Color','k','LineWidth',2)
bd=ylabel({'Theta','Filtered','LFP'},'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize,'Rotation',0,'HorizontalAlignment','Right');
set(bd, 'Units', 'Normalized', 'Position', [-0.03, 0.3, 0]);
set(bd,'Units','character')
box off
set(gca,'XTick',[],'XTickLabel',{},'XColor','w')
ylim([min(filteredlfp(tidx,2))-(max(filteredlfp(tidx,2))-min(filteredlfp(tidx,2)))*.1 max(filteredlfp(tidx,2))+(max(filteredlfp(tidx,2))-min(filteredlfp(tidx,2)))*.1])
hold on

disp(['Scale difference: 1 mV on ngf scale = ' num2str(ngfdiff/(.95*(1-handles.formatP.bottom*2))*diff(ylim)) ' mV on Theta LFP trace'])

set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
set(gca,'Clipping','off','ycolor','k','xcolor','k','YTick',[],'YTickLabel',{})
myax=gca;
ax2 = axes('Position', get(gca, 'Position'),'Color','none');
set(ax2,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','off','layer','top','LineWidth',2)
yL=get(myax,'YLim');
xL=get(myax,'XLim');
hold(myax,'on')
plot(myax,[xL(2)-3 xL(2)-3],[yL(1)+.02 yL(1)+.2+.02],'k','LineWidth',2)
plot(myax,[xL(2)-103 xL(2)-3],[yL(1)+.02 yL(1)+.02],'k','LineWidth',2)
disp('LFP scale = .2 mV on theta LFP')
%text(xL(2)-30,yL(2),{'100 ms','5 mV'},'HorizontalAlignment','right','FontName','ArialMT','FontSize',myFontSize)

filteredlfp=mikkofilter(handles.curses.lfp,1000/handles.runinfo.lfp_dt,[25 40]);
if olderflag==1
    h(4)=figure('Renderer', 'painters','Visible','on','Color','w','PaperUnits','inches','PaperSize',[myw lfpht],'PaperPosition',[0 0 myw lfpht],'Name','Filtered LFP Gamma','Units','inches','Position',[.5 .5 myw lfpht]);
else
    h(4)=figure('GraphicsSmoothing','off', 'Renderer', 'painters','Visible','on','Color','w','PaperUnits','inches','PaperSize',[myw lfpht],'PaperPosition',[0 0 myw lfpht],'Name','Filtered LFP','Units','inches','Position',[.5 .5 myw lfpht]);
end
subplot('Position',[handles.formatP.left*2 handles.formatP.bottom .9-handles.formatP.left .95*(1-handles.formatP.bottom*2)]);
tidx=find(filteredlfp(:,1)>=zoomrange(1) & filteredlfp(:,1)<=zoomrange(2));
plot(filteredlfp(tidx,1),filteredlfp(tidx,2),'Color','k','LineWidth',2)
bd=ylabel({'Gamma','Filtered','LFP'},'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize,'Rotation',0,'HorizontalAlignment','Right');
set(bd, 'Units', 'Normalized', 'Position', [-0.03, 0.3, 0]);
set(bd,'Units','character')
box off
set(gca,'XTick',[],'XTickLabel',{},'XColor','w')
ylim([min(filteredlfp(tidx,2))-(max(filteredlfp(tidx,2))-min(filteredlfp(tidx,2)))*.1 max(filteredlfp(tidx,2))+(max(filteredlfp(tidx,2))-min(filteredlfp(tidx,2)))*.1])
hold on

disp(['Scale difference: 1 mV on ngf scale = ' num2str(ngfdiff/(.95*(1-handles.formatP.bottom*2))*diff(ylim)) ' mV on Gamma LFP trace'])

set(gca,'FontName',myFontName,'FontWeight',myFontWeight,'FontSize',myFontSize)
set(gca,'Clipping','off','ycolor','k','xcolor','k','YTick',[],'YTickLabel',{})
myax=gca;
ax2 = axes('Position', get(gca, 'Position'),'Color','none');
set(ax2,'XTick',[],'YTick',[],'XColor','w','YColor','w','box','off','layer','top','LineWidth',2)
disp(['LFP scale = ' num2str(diff(get(myax,'ylim'))/diff(yL)*.2) ' mV on gamma LFP'])

% get(myax,'YLim') yL);
% set(myax,'XLim',xL);
% hold(myax,'on')
% plot(myax,[xL(2)-3 xL(2)-3],[yL(1)+.02 yL(1)+.2+.02],'k','LineWidth',2)
% plot(myax,[xL(2)-103 xL(2)-3],[yL(1)+.02 yL(1)+.02],'k','LineWidth',2)
% text(xL(2)-30,yL(2),{'100 ms','.2 mV'},'HorizontalAlignment','right','FontName','ArialMT','FontSize',myFontSize)

fid=fopen([savepath sl 'LFP.txt'],'w');
ploffset=(length(handles.curses.lfp)-length(thetalfp))/2;
fprintf(fid,'Time\tRaw\n');
for f=1:length(thetalfp)
    fprintf(fid,'%f\t%f\n', handles.curses.lfp(f+ploffset,1), handles.curses.lfp(f+ploffset,2));
end
fclose(fid);

fid=fopen([savepath sl 'FilteredLFP.txt'],'w');
ploffset=(length(handles.curses.lfp)-length(thetalfp))/2;
fprintf(fid,'Time\tTheta\tGamma\n');
for f=1:length(thetalfp)
    fprintf(fid,'%f\t%f\t%f\n', handles.curses.lfp(f+ploffset,1), thetalfp(f,2), filteredlfp(f,2));
end
fclose(fid);

fid=fopen([savepath sl 'Membrane_Potentials.txt'],'w');
fprintf(fid,'Time');
mystr='%f\t';
myvars='mytrace.(mycells{1}).trace(q,1),';
for m=1:length(mycells)
    mystr=[mystr '%f\t'];
    myvars=[myvars 'mytrace.' mycells{m} '.trace(q,2),'];
    fprintf(fid,'\t%s',mycells{m});
end
fprintf(fid,'\n');

mystr=[mystr(1:end-1) 'n'];
myvars=myvars(1:end-1);

for q=1:length(mytrace.(mycells{1}).trace(:,1))
    eval(['fprintf(fid,''' mystr ''',' myvars ');']);
end
fclose(fid);
    %subplot(length(mycells),1,mytrace.(mycells{m}).pos)
    



function trace=loadmyfile(handles,resultname,celltype)
global sl repodir
    trace=[];
    fd=strmatch(celltype,{handles.curses.cells.name});
    myfile='';
    r=1;
    gg=0;
    while ~isempty(fd) && r<=length(handles.curses.cells(fd).mygids)
        if exist([repodir sl 'results' sl resultname sl 'trace_' celltype num2str(handles.curses.cells(fd).mygids(r)) '.dat'],'file')==2
            gg=gg+1;
            if isempty(myfile)
                myfile=[repodir sl 'results' sl resultname sl 'trace_' celltype num2str(handles.curses.cells(fd).mygids(r)) '.dat'];
            end
        end
        r=r+1;
    end
    if isempty(myfile)
        disp(['no ' celltype ' cells within area were recorded, using a recording from outside the Area of Interest!!'])
    else
        trace=importdata(myfile);
        disp(['Out of ' num2str(gg) ' recordings for ' num2str(r-1) ' '  celltype 's in area, Importing: ' myfile])
    end
    
    
function h=performancefig(FinalData)
global disspath savepath repodir sl

disspath = savepath;
    dd=strmatch('ca1_olmVar_gH_30_03',{FinalData.Name},'exact');
    if isempty(dd)
        disp('ca1_olmVar_gH_30_03 needs to be in FinalData')
        return
    end

    %ncrunmem=[FinalData(dd).memvec.memstruct([1 3 66 69]).VIRT]/10^6;
if exist([repodir sl 'networkclamp_results' sl 'ca1_olmVar_gH_30_03' sl '00002' sl 'topoutput.dat'],'file')
    filename = [repodir sl 'networkclamp_results' sl 'ca1_olmVar_gH_30_03' sl '00002' sl 'topoutput.dat'];
    fid = fopen(filename);
    c=textscan(fid,'%s','Delimiter','\t');
    fclose(fid);
    DivBy=1024*1024;
    ofint=[1 2 13 16];
    tempmem=[];
    for p=1:length(ofint)
        medrr=regexp(c{:}{p*2},'\s*','split');
        switch medrr{5}(end)
            case 'm'
                tempmem(p) = str2num(medrr{5}(1:end-1))*1024; % get into kb (from mb)
            case 'g'
                tempmem(p) = str2num(medrr{5}(1:end-1))*1048576; % get into kb (from gb)
            otherwise
                tempmem(p) = str2num(medrr{5}(1:end-1));  % assume kb
        end
    end


    ncrunmem=tempmem/DivBy; % VIRT is column 5



    yL='Memory (GB)';
    mf=perfgraphs([FinalData(dd).memvec.memstruct([1 3 66 69]).VIRT]*FinalData(dd).NumProcessors/(1024^2),ncrunmem,yL,'grouped');
    set(mf,'Position',[.5 .5 4 3],'PaperPosition',[0 0 4 3],'PaperSize',[4 3])
    setFonts(mf) 
    printeps(mf,[disspath 'ModelTypeMemory'])
else
    disp('Network clamp results not available to grab memory from, using hard coded value')
    yL='Memory (GB)';
    mf=perfgraphs([392.11 1917.89 3013.83 3000.70],ncrunmem,yL,'grouped');
    set(mf,'Position',[.5 .5 4 3],'PaperPosition',[0 0 4 3],'PaperSize',[4 3])
    setFonts(mf) 
    printeps(mf,[disspath 'ModelTypeMemory'])
end
ncrun=[27.22 1.75 370.17 192.53 0.23];

yL='Wall-Clock Time (s)';
mf=perfgraphs(FinalData(dd).timevec,ncrun,yL,'stacked');
    setFonts(mf) 
printeps(mf,[disspath 'ModelTypeWallTimes'])


yL='CPU Time (s)';
mf=perfgraphs(FinalData(dd).timevec*FinalData(dd).NumProcessors,ncrun,yL,'stacked');
    setFonts(mf) 
printeps(mf,[disspath 'ModelTypeCPUTimes'])

h=getperf(FinalData);


function setFonts(tt) 
global myFontWeight myFontName myFontSize

bf = findall(tt,'Type','text');
for b=1:length(bf)
    set(bf(b),'FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
end
bf = findall(tt,'Type','axis');
for b=1:length(bf)
    set(bf(b),'FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
end
bf = findall(tt,'Type','axes');
for b=1:length(bf)
    set(bf(b),'FontWeight',myFontWeight,'FontName',myFontName,'FontSize',myFontSize)
end
