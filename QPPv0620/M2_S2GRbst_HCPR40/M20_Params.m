
%% Setting paths, parameters & indices
%% 
clear; clc;
nmD='HCPR40'; % dataset name 
p2prep=... % path to (pth2/p2) preproc params, see ../../All_Pre~/P0_Pre~.m 
    ['../../All_Preproc_Data/Prep_Params_' nmD '.mat'];

%% For saving & reading path-names
nmMD=['S2GQPP_Rbst_' nmD]; nmMDt=[nmMD '.mat']; % method & dataset name
p2p=['Params_' nmMDt]; % p2 SAVE SETTINGS HERE >> SHOULD BEGIN WITH Params_
p2O=nmMDt; % pth2 save primary outputs (>> mkdir if needed)
d2O='./Fls/'; % dir2 save secondary outputs or temporary variables
if ~exist(d2O,'dir'), mkdir(d2O); end

load(p2prep,'p2BGWCR','p2bgwcr','nsbj','nscn'); a=strfind(p2prep,'/'); 
p2B=[p2prep(1:a(end)) p2BGWCR]; % pth2 parcellated scans
p2bpr=cell(nsbj,nscn); % pth2 preproc-ed grayordinate scans
for is=1:nsbj, for i=1:nscn
        p2bpr{is,i}=[p2prep(1:a(end)) p2bgwcr{is,i}(3:end)]; end; end

p2qppf='../QPPfv0620/'; addpath(p2qppf); % pth2 QPP functions
p2u='../../Utils/'; addpath(genpath(p2u)); % pth2 cifti utilities

%% For detecting QPPs of subjects (which also involves regressing QPPs)
PL=30; % QPP pattern length, ~20s in humans, 30 timepoints wiht TR=0.72s
PLh=round(PL/2)+[0 -rem(PL,2)]; % pad length for temporally extending a QPP
PLc=round(PL/2)+(1:PL); % range of interest in an extended QPP, matches PLh
PLe=PL+sum(PLh); % length of an extended QPP (derivable but saved/read)

nP=5; % #QPPs to detect
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

nt2PL=round(nt/PL); % a way for choosing limited #initial-segments/scan
nITPps=1*nt2PL; % #initial-segments/scan 
ITPl=cell(nsbj,nP); % random selection of nITPps/scan
for is=1:nsbj
for ip=1:nP
for i=1:nscn
    a=(i-1)*nt+(ssg(ip):esg); a=a(randperm(length(a))); a=a(1:nITPps); 
    ITPl{is,ip}=[ITPl{is,ip}; single(a')]; 
end
end
end
ITPlClk=datestr(now); ITPlClk(end-2:end)=[]; % clocking random ITP to log

tsh=floor(PL/4); % max timeshift when comparing QPPs or any 2 templates 
cthph=0.88; % similarity threshold when phase-adjusting (phadj) a QPP
load('myGlssr.mat','ilv2','ilpcc5','ilpf','ilm1'); % set of seeds for phadj
sdph=cell(nP,1); sdph{1}={ilv2;ilpcc5}; sdph{2}={ilv2;ilpf}; 
sdph{3}={ilv2;ilm1}; sdph(4:end)={{ilv2}}; clear ilv2 ilpcc5 ilpf ilm1
nsd=zeros(nP,1); for ip=1:nP, nsd(ip)=length(sdph{ip}); end; nSD=max(nsd);

nCNT=3; % #QPPs to included when counting transitions

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

%% For obtaining QPPs in grayordinate & summary maps
p2V=['Vx' nmMDt]; % pth2 save GrpQPPs & summaries in grayordinates
% p2VG2S=cell(nP,nSD); % p2 save GrpQPPs project back to Sbjs for statistics
% for ip=1:nP, for i=1:nsd(ip), p2VG2S{ip,i}=[d2O  ...
%             'VxGrp2SbjQPP' num2str(ip) '_SD' num2str(i) '.mat']; end; end
p2NV=[d2O 'NllVx' nmMDt]; % pth2 save null patterns in grayordinates

load('myHCPcft.mat','irgn','nrgn','nmrgn1');% indices for brain regions
load('myGlssr.mat',... % grayordinate indices of areas used as ref.s to ...
    'ivxlpcc5','ivxlv2','ivxlm1','ivxlpf'); %  "describe" btwn regions
irf={ivxlpcc5; ivxlv2; ivxlm1; ivxlpf}; nrf=length(irf);

tshclst=2; % allowed timeshift for clustering GrpQPP's timecourses
nclst=10; % to simplify, 1st few largest clusters kept after h. clustering
clstth=0.1; % cutt-off for h. clustering
clstclr=[.6 0 .5; .25 .25 .5; % cluster colors as workbench Videen wointerp
    .37 .37 .75; .5 .5 1; 0 1 0; 0 .75 0; 1 1 0; 1 .75 0; 1 .5 0; 1 0 0]; 

load('myPrcls','PR','nmPR','nPR'); % existing parcellations see [p2u 'U3~']
ntst=50; % number of repetition when matching group sizes in paired t-tests

%% For saving output plots
d2plt='./Plts/'; % dir2 save plots
if ~exist(d2plt,'dir'), mkdir(d2plt); end 
d2cft=[d2plt 'cfts/']; % dir2 save cifti files
if ~exist(d2cft,'dir'), mkdir(d2cft); end

d2SA='./Fls/';  % d2 save SA files used by codes in SA fldr
a=['./SA' d2SA(2:end)]; if ~exist(a,'dir'), mkdir(a); end
d2SAplt='./Plts/'; % ~ SA plots
a=['./SA' d2SAplt(2:end)]; if ~exist(a,'dir'), mkdir(a); end
d2SAcft=[d2SAplt 'cfts/']; % ~ SA cifti files
a=['./SA' d2SAcft(2:end)]; if ~exist(a,'dir'), mkdir(a); end

%% Clearing extras, loading needed from prep.mat & saving all set here
clear a i ip is nmD nmMD nmMDt p2BGWCR p2bgwcr
load(p2prep,'nT','ivx','nvx','nVX','nX','tres','INFO'); % nsbj,nscn,nt
save(p2p);
