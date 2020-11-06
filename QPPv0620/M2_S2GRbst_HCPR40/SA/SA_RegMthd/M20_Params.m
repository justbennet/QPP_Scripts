
%% Setting paths, parameters & indices
%% 
clear; clc;
nmD='HCPR40'; % dataset name 
p2prep=... % path to (pth2/p2) preproc params, see ../../All_Pre~/P0_Pre~.m 
    ['../../../../All_Preproc_Data/Prep_Params_' nmD '.mat'];

%% For saving & reading path-names
nmMD=['S2GQPP_Rbst_Sgm_' nmD]; nmMDt=[nmMD '.mat']; % method & dataset name
p2p=['Params_' nmMDt]; % p2 SAVE SETTINGS HERE >> SHOULD BEGIN WITH Params_
p2O=nmMDt; % pth2 save primary outputs (>> mkdir if needed)
d2O='./Fls/'; % dir2 save secondary outputs or temporary variables
if ~exist(d2O,'dir'), mkdir(d2O); end

load(p2prep,'p2BGWCR','p2bgwcr','nsbj','nscn'); a=strfind(p2prep,'/'); 
p2B=[p2prep(1:a(end)) p2BGWCR]; % pth2 parcellated scans
p2bpr=cell(nsbj,nscn); % pth2 preproc-ed grayordinate scans
for is=1:nsbj, for i=1:nscn
        p2bpr{is,i}=[p2prep(1:a(end)) p2bgwcr{is,i}(3:end)]; end; end

p2qppf='../../../QPPfv0620/'; addpath(p2qppf); % pth2 QPP functions
p2u='../../../../Utils/'; addpath(genpath(p2u)); % pth2 cifti utilities

%% For detecting QPPs of subjects (which also involves regressing QPPs)
PL=30; % QPP pattern length, ~20s in humans, 30 timepoints wiht TR=0.72s
PLh=round(PL/2)+[0 -rem(PL,2)]; % pad length for temporally extending a QPP
PLc=round(PL/2)+(1:PL); % range of interest in an extended QPP, matches PLh
PLe=PL+sum(PLh); % length of an extended QPP (derivable but saved/read)

nP=3; % #QPPs to detect
cth=cell(nP,1); % correlation threshold in the main algorithm
cth(1:min(nP,3))={[0.1 0.2]}; cth(4:end)={[0.1 0.2]};
ncth1=3*ones(nP,1); if nP>=4, ncth1(4:end)=3; end % #iters with lower cth
nitr=15; % max #iterations in the main algorithm

load(p2prep,'nt');
esg=nt-PL+1; % end segment/scan when correlating a QPP template
ssg=ones(nP,1); ssg(2:end)=PL; % starting segment/scan ~, SEE QPPf4regscn
ITP=cell(nP,1); % all possible initial-segments/sbj for robust QPPdetection
for ip=1:nP
for i=1:nscn
    a=(i-1)*nt+(ssg(ip):esg); 
	ITP{ip}=[ITP{ip}; single(a')]; 
end
end
ITPstp=50; % step to show progress of algorithm when running ITP times

tsh=floor(PL/4); % max timeshift when comparing QPPs or any 2 templates 
cthph=0.88; % similarity threshold when phase-adjusting (phadj) a QPP
load('myGlssr.mat','ilv2','ilpcc5','ilpf','ilm1'); % set of seeds for phadj
sdph=cell(nP,1); sdph{1}={ilv2;ilpcc5}; sdph{2}={ilv2;ilpf}; 
sdph{3}={ilv2;ilm1}; sdph(4:end)={{ilv2}}; clear ilv2 ilpcc5 ilpf ilm1
nsd=zeros(nP,1); for ip=1:nP, nsd(ip)=length(sdph{ip}); end; nSD=max(nsd);

a=[d2O 'SbjQPP/']; % dir2 save SbjQPPs
if ~exist(a,'dir'), mkdir(a); end
p2S=cell(nsbj,nP); for is=1:nsbj, for ip=1:nP % pth2 save SbjQPPs
	p2S{is,ip}=[a 'S' num2str(is) '_' num2str(ip)]; end; end

nnll=15; % number of null patterns for significance testing of QPP activity

%% For functional connectivity (FC)
load('myGlssr.mat',... % indices to reorder Glasser's parcels based on ...
   'nY','G2Y','ibY','iG2Y','YLB'); %  Yeo's 7 RSNs, see [p2u 'U2myGlssr.m']
fcbn=-1:0.01:1; % histogram bins for FC
fcth=0.1; % FC-values with magnitude > fcth are found after regressing QPPs

%% For obtaining group QPPs
nitrg=1; % #iters in group averaging to tune the prior template
p2G=cell(nP,nSD); % pth2 save secondary outputs for group QPPs
for ip=1:nP, for i=1:nsd(ip), p2G{ip,i}=[d2O  ...
            'GrpQPP' num2str(ip) '_SD' num2str(i) '.mat']; end; end

%% For saving output plots
d2plt='./Plts/'; % dir2 save plots
if ~exist(d2plt,'dir'), mkdir(d2plt); end

%% Clearing extras, loading needed from prep.mat & saving all set here
clear a i ip is nmD nmMD nmMDt p2BGWCR p2bgwcr
load(p2prep,'nT','ivx','nvx','nVX','nX','tres','INFO'); % nsbj,nscn,nt
save(p2p);
