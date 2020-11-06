
%% Setting preprocessing parameters, paths, etc, compatible with QPPv0420
%% Dataset name, sbj IDs, etc
clear; clc; nmD='HCPR3'; % dataset name
p2prep=... % path-name to (pth2/p2) SAVE all params, paths, etc, set here
    ['Prep_Params_' nmD '.mat']; %% >> KEEP THIS FORMAT OR CHANGE P1_P~.m
ID=dlmread(['ID_' nmD '.txt']); nsbj=size(ID,1); % IDs & consequent #sbj

nD=2; nPD=2; % #days, #scans/day (>>HCP-BASED)
iD=1:nD; iPD=nPD; % indices of day&scan to process (all scans iPD=1:nPD;)
nscn=length(iD)*length(iPD); % conseq. #scans/sbj to process
nscng=nsbj*nscn; % conseq. #scans/grp

p2u='../Utils/'; addpath(genpath(p2u)); % pth2 cifti-/nii-related utilities

%% Pth2 read scans & WM&CSF signals (nuisance regressors) & get params/scan
p1='/keilholz-lab/SharedFiles/HCP_S900_rsfMRI';
p2=cell(nsbj,1);
for i=1:nsbj   
    if ID(i,2)==100, p2{i}='/HCP100';
    elseif ID(i,2)==500, p2{i}='/HCP500minus100';
    elseif ID(i,2)==900, p2{i}='/HCP900minus500';
    end
end
p3='_RSfMRI_1and2_FIXDenoised_Extended1/';
p4='/MNINonLinear/Results/rfMRI_REST';
RLLR={'RL','LR','LR','RL'}; RLLR=repmat(RLLR,nsbj,1); % order of scans/day
for i=1:nsbj, if ID(i,3)==1, RLLR{i,3}='RL'; RLLR{i,4}='LR'; end; end
ID=single(ID(:,1)); % no more needing columns 2-3
p5='/rfMRI_REST';
nmb='_Atlas_hp2000_clean.dtseries.nii'; % p1-p5 & nmb for BOLD fMRI scans
p6='/RestingStateStats'; 
p7='_Atlas_Cleaned';
nmrg={'WMtc.txt','CSFtc.txt'}; % p1-p5 & p6-7 & frg for nuisance regressors

% is=1; id=1; ipd=1; % >> sbj1scan1 to CHECK reading path & get params/scan
% iscn=(id-1)*nD+ipd;
% pc=[p1 p2{is} p3 num2str(ID(is)) p4 num2str(id) '_' RLLR{is,iscn} ];
% pb=[pc p5 num2str(id) '_' RLLR{is,iscn}]; 
% b=ft_read_cifti([pb nmb]); 
% b0=b.dtseries;
% 
% [nVX,nt]=size(b0); % >> spatial & temporal dim, >> HCP scans are 2D
% ivx=find(~isnan(b0(:,1))); % >> indices of brain areas (ft_read_cifti)
% nvx=length(ivx); 
% tres=b.time(2)-b.time(1); clear b b0 % >> TR(s)

nt=1200; tres=0.72; % another way to set params/scan 
load('myHCPcft.mat','nVX','ivx','nvx'); % see [p2u 'U1myHCPcft.m']

nT=nt*nscn; nTg=nt*nscng; % conseq. # all-timepoints/sbj & /grp

%% Directory to save preprocessed scans in grayordinates
d2prep={['./WCR_' nmD '/'],... % WM&CSF Regressed scans
['./GWCR_' nmD '/']};  % Gray-Matter(GM)&WM&CSF Regressed scans
for i=1:2, if ~exist(d2prep{i},'dir'), mkdir(d2prep{i}); end; end

%% Putting path-names to read&save grayordinate scans in cell arrays
c=cell(nsbj,nscn); p2b=c; p2rg=c; p2bwcr=c; p2bgwcr=c; IDF=c; clear c
for is=1:nsbj
    iscn=1; IDS=num2str(ID(is));
for id=iD
for ipd=iPD
    a=[num2str(id) '_' RLLR{is,(id-1)*nD+ipd}];
    pc=[p1 p2{is} p3 IDS p4 a];
    pb=[pc p5 a]; p2b{is,iscn}=[pb nmb]; % pth2 read scans
    prg=[pc p6 p5 a p7];
    for i=1:2, p2rg{is,iscn}{i}=[prg nmrg{i}]; end % pth2 read regressors
    a=[IDS '_' num2str(id) RLLR{is,(id-1)*nD+ipd}];
    IDF{is,iscn}=a; a=[a '.mat'];
    p2bwcr{is,iscn}=[d2prep{1} a]; p2bgwcr{is,iscn}=[d2prep{2} a]; % p2 save
    iscn=iscn+1;
end
end
end; clear p1 p2 p3 p4 p5 p6 p7 pc pb prg nmb nmrg a ...
    RLLR ID IDS d2prep i is id ipd iscn nD nPD
INFO.IDF=IDF; INFO.iD=iD; INFO.iPD=iPD; clear IDF iD iPD

%% Bandpass filter & zeropad size
h=fdesign.bandpass('N,F3dB1,F3dB2',4,0.008,0.125,1/tres); % F1dB: 0.01&0.1
H=design(h,'butter','SystemObject', true);
Fltr{1}=H.SOSMatrix; Fltr{2}=H.ScaleValues; clear h H
nzFltr=(2^nextpow2(nt)-nt)/2; % for nt=1200, nzFltr=~400, which is good

%% Indices of Glasser's parcellation
load('myGlssr.mat','ixG'); nX=length(ixG); % see [p2u 'U2myGssr.m']

%% Pth2 save parcellated scans
p2BWCR=['B_WCR_' nmD '.mat']; p2BGWCR=['B_GWCR_' nmD '.mat'];

%% Clearing extras & SAVING all 
clear p2u nmD
save(p2prep);
