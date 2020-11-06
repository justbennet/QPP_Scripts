
%% Parcellations: reading nii files & saving into cifti format
%%
clear; clc; close all; th1=0.1; th2=10; 
load('../myHCPcft.mat','sc','nVX','nrgn'); 
sc=[{[]};sc]; PR=cell(nrgn,1); nmPR=cell(nrgn,1);
addpath '../NIfTI_20140122/' '../cifti-matlab-master/';

e=ft_read_cifti('../empty.dtseries.nii'); 
ps=e.pos-diag(e.transform(1:3,1:3))';
st=e.brainstructure;

%% Cerebral cortex (Yeo et al 2011)
%% http://surfer.nmr.mgh.harvard.edu/fswiki/
%%
irgn=1; 
load('../myGlssr.mat','Y','YLB'); nmPR{irgn}=YLB'; 
f2s='A1_Cerebrum_Yeo';
PR{irgn}=Y;

e.dtseries=Y; 
e.dtseries(Y==3)=4; e.dtseries(Y==4)=3; % for DA&VA colors in workbench
e.time=1; e.hdr.dim(6)=1; 
ft_write_cifti(f2s,e,'parameter','dtseries');

%% Cerebellum (Buckner et al 2011)
%% http://surfer.nmr.mgh.harvard.edu/fswiki/CerebellumParcellation_Buckner2011
%%
irgn=2; 
nmPR{irgn}=YLB'; 
nPR=length(nmPR{irgn});
f2r='Buckner2011_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz';
f2s='A2_Cerebellum_Buckner';

A2=myfn1(f2r,nVX,nPR,ps,st,sc{irgn}); 
PR{irgn}=A2;

e.dtseries=A2; 
e.dtseries(A2==3)=4; e.dtseries(A2==4)=3; % for DA&VA colors in workbench
e.time=1; e.hdr.dim(6)=1; 
ft_write_cifti(f2s,e,'parameter','dtseries');
 
%% Thalamus (Behrens et al 2003)
%% https://www.lead-dbs.org/ - After downloading, based on v2018, go to:
% /lead/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/Oxford Thalamic
% Connectivity Atlas (Behrens 2003)/mixed/
%%
irgn=3;
nmPR{irgn}={'PFC';'PreM';'M1';'S1';'PPC';'OCL';'TL'}; 
nPR=length(nmPR{irgn});
p0='./Oxford Thalamic Connectivity Atlas (Behrens 2003)/'; 
p1={'prefrontal';'premotor';'primarymotor';'sensory';'postparietal';
    'occipital';'temporal'};
f2r=cell(nPR,1); for i=1:nPR, f2r{i}=[p0 p1{i} '.nii.gz']; end
f2s='A3_Thalamus_Behrens';

AT=myfn2(f2r,nVX,nPR,ps,st,sc{irgn},th2); 
[~,imx]=max(AT,[],2);
AT(isnan(AT))=0; io=find(sum(AT,2)>0);
A2=nan(nVX,1,'single'); A2(io)=imx(io); 
PR{irgn}=A2;

e.dtseries=A2; e.time=1; e.hdr.dim(6)=1; 
ft_write_cifti(f2s,e,'parameter','dtseries');

%% Thalamus (Brainnetome, Fan et al 2016)
%% https://atlas.brainnetome.org/bnatlas.html
%%
% % irgn=3; 
% % nmPR{irgn}={'LmPF';'RmPF';'LPreM';'RPreM';'LS1';'RS1';'LaTL';'RaTL';
% %     'LPP';'RPP';'LOC';'ROC';'LpTL';'RpTL';'LlPF';'RlPF'}; shft=230; 
% % nPR=length(nmPR{irgn});
% % f2r='BN_Atlas_246_2mm.nii.gz';
% % f2s='A3_Thalamus_BNT';
% % 
% % A2=myfn1p5(f2r,nVX,nPR,shft,ps,st,sc{irgn}); 
% % PR{irgn}=A2;
% % 
% % %%%% reordering just for plotting: 1)PF 2)PreM 3)S1 4)PP 5)OC 6)TL
% % a=round(A2/2); b=a; b(a==8)=1; b(a==4|a==7)=6; b(a==5)=4; b(a==6)=5; 
% % e.dtseries=b; e.time=1; e.hdr.dim(6)=1; 
% % ft_write_cifti(f2s,e,'parameter','dtseries');

%% Hippocampus (Robinson et al HBM 2016)
%% http://aucanlab.com/research/data/ 
%%
irgn=4;
% nmPR{irgn}={'left_anterior'; 'left_mid'; 'left_posterior'; 'right_anterior';
% 'right_anteromedial'; 'right_anterolateral'; 'right_mid'; 'right_posterior'}; 
nmPR{irgn}={'LA'; 'LM'; 'LP'; 'RA';'RAM'; 'RAL'; 'RM'; 'RP'}; 
nPR=length(nmPR{irgn});
p0='./HippocampusSegmentations/';
f2r=cell(nPR,1); for i=1:3, f2r{i}=[p0 'left' num2str(i) '.nii.gz']; end
for i=4:nPR, f2r{i}=[p0 'right' num2str(i-3) '.nii.gz']; end
f2s='A4_Hipppcampus_Robinson'; 

AT=myfn2(f2r,nVX,nPR,ps,st,sc{irgn},th1); 
[~,imx]=max(AT,[],2);
AT(isnan(AT))=0; io=find(sum(AT,2)>0);
A2=nan(nVX,1,'single'); A2(io)=imx(io); 
PR{irgn}=A2;

e.dtseries=A2; e.time=1; e.hdr.dim(6)=1; 
ft_write_cifti(f2s,e,'parameter','dtseries');

%% Amygdala (Tyszka & Pauli 2016)
%% https://osf.io/dxsz5/
%%
irgn=5;
nmPR{irgn}={'La'; 'BL'; 'BM'; 'BLV'; 'CEN-CMN'; 'ASTR-ATA-AAA-AMY'};
nPR=length(nmPR{irgn});
iPR={1; 2; 3; 6; 4:5; 7:10}; % combined parcels named in nmPR
for i=1:nPR, nmPR{irgn}{i}=['P' num2str(i)]; end % renaming to simplify
f2s='A5_Amygdala_Pauli';
load('../myHCPcft.mat','ivx','nvxc');

% addpath '../cifti-gifti/'; % >> needs workbench installed in C:/workbench
% f2r='CIT168_AmyNuc.dlabel.nii';
% A=ciftiopen(f2r,'C:/workbench/bin_windows64/wb_command');
% save('CIT168_AmyNuc','A');
load('CIT168_AmyNuc.mat','A');

A1=nan(nVX,1); ivxsc=ivx(nvxc+1:end); A1(ivxsc)=A.cdata;
A2=nan(nVX,1,'single');
for i=1:nPR, for j=1:length(iPR{i}), A2(A1==iPR{i}(j))=i; end; end
ia=find(~isnan(A2));
st1=st(ia); u=unique(st1); [~,iu,~]=intersect(u,sc{irgn}); u(iu)=[];
ix=[]; for i=1:length(u), ix=[ix; find(st1==u(i))]; end
A2(ia(ix))=nan; clear A A1 ivxsc ia st1 u iu ix
PR{irgn}=A2;

e.dtseries=A2; e.time=1; e.hdr.dim(6)=1; 
ft_write_cifti(f2s,e,'parameter','dtseries');

%% Brainstem & deep brain nuclei
%%
irgn=6;
nmPR{irgn}={'LC';'VTA';'SNc';'SNr';'BNM';'PPN';'DR';'MR';'PAG';'STH';'GPe';
    'GPi';'HTH'}; 
nPR=length(nmPR{irgn});
sci=sc{irgn};
AT=nan(nVX,nPR,'single');
f2s='A6_Brainstem_Deepbrain_Nuclei';

%% LC (Keren et al 2009)
%% http://eckertlab.org/lc/ (needs registration by emailing)
%%
f2r={'LC_2SD_BINARY_TEMPLATE.nii'};
AT(:,1)=myfn2(f2r,nVX,1,ps,st,sci,th1);

%% VTA,SNc,SNr,STH,GPe,GPi&HTH (Pauli et al 2018)
%% https://www.lead-dbs.org/
% /lead/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/CIT168_Reinf_Learn
% (Pauli 2017) [imported from CIT168]/mixed/
%%
p0='./CIT168_Reinf_Learn (Pauli 2017) [imported from CIT168]/'; 
p1={'VTA';'SNc';'SNr';'STH';'GPe';'GPi';'HTH';'RN'}; 
n=length(p1);
f2r=cell(n,1); for i=1:n, f2r{i}=[p0 p1{i} '.nii.gz']; end
AT(:,[2:4 10:14])=myfn2(f2r,nVX,n,ps,st,sci,th1);

%% BNM (Li et al 2014, Zaborszky group)
%% Shared by first author after emailing (similar to Xiao Liu et al 2018)
%%
f2r={'BNM.nii'};
AT(:,5)=myfn2(f2r,nVX,1,ps,st,sci,th1);

%% PPN,DR,MR&PAG (Eldow et al 2012)
%% https://www.lead-dbs.org/
% /lead/templates/space/MNI_ICBM_2009b_NLIN_ASYM/atlases/Harvard AAN atlas
% (Edlow 2012)/mixed/
%%
p0='./Harvard AAN atlas (Edlow 2012)/';
p1={'PPN';'DR';'MR';'PAG'}; 
n=length(p1);
f2r=cell(n,1); for i=1:n, f2r{i}=[p0 p1{i} '.nii.gz']; end
AT(:,6:9)=myfn2(f2r,nVX,n,ps,st,sci,th1);

%% Combining all brainstem & deep brain nuclei
%%
[~,imx]=max(AT,[],2);
AT(isnan(AT))=0;
io=find(sum(AT,2)>0);
A2=nan(nVX,1,'single'); A2(io)=imx(io); 
PR{irgn}=A2;

e.dtseries=A2; 
e.time=1; e.hdr.dim(6)=1; 
ft_write_cifti(f2s,e,'parameter','dtseries');

%% Striatum (Choi et al 2012)
%% http://surfer.nmr.mgh.harvard.edu/fswiki/StriatumParcellation_Choi2012
%%
irgn=7;
nmPR{irgn}=YLB'; 
nPR=length(nmPR{irgn});
f2r='Choi2012_7Networks_MNI152_FreeSurferConformed1mm_LooseMask.nii.gz';
f2s='A7_Striatum_Choi';

A2=myfn1(f2r,nVX,nPR,ps,st,sc{irgn}); 
PR{irgn}=A2;

e.dtseries=A2; 
e.dtseries(A2==3)=4; e.dtseries(A2==4)=3; % for DA&VA colors in workbench
e.time=1; e.hdr.dim(6)=1; 
ft_write_cifti(f2s,e,'parameter','dtseries');

%%
%%
nPR=zeros(nrgn,1,'single'); for i=1:nrgn, nPR(i)=length(nmPR{i}); end
save('../myPrcls.mat','PR','nmPR','nPR');

%% 
%%
function A2=myfn1(f2r,nVX,nA,ps,st,sci)
A=load_nii(f2r); A1=A.img; 
AO=A.hdr.hist.originator(1:3); AD=A.hdr.dime.pixdim(2);
A2=nan(nVX,1,'single');
for i=1:nA
    [I,J,K]=ind2sub(size(A1),find(A1==i));
    clear a; a(:,1)=I; a(:,2)=J; a(:,3)=K; 
    a=(a-AO)*AD;
    [~,~,ips]=intersect(a,ps,'rows');
    st1=st(ips); u=unique(st1); [~,iu,~]=intersect(u,sci); u(iu)=[];
    ix=[]; for j=1:length(u), ix=[ix; find(st1==u(j))]; end; ips(ix)=[];
    A2(ips)=i;
end
end

function A2=myfn2(f2r,nVX,n,ps,st,sci,th)
    A2=nan(nVX,n,'single');
for i=1:n
    A=load_nii(f2r{i}); A1=A.img;
    AO=A.hdr.hist.originator(1:3); AD=A.hdr.dime.pixdim(2);
    iA1=find(A1>=th);
    [I,J,K]=ind2sub(size(A1),iA1);
    clear a; a(:,1)=I; a(:,2)=J; a(:,3)=K;
    a=round((a-AO)*AD);
    [~,ia,ips]=intersect(a,ps,'rows');
    st1=st(ips); u=unique(st1); [~,iu,~]=intersect(u,sci); u(iu)=[];
    ix=[]; for j=1:length(u), ix=[ix; find(st1==u(j))]; end
    ips(ix)=[]; ia(ix)=[];
    A2(ips,i)=A1(iA1(ia));
end
end

% % function A2=myfn1p5(f2r,nVX,nA,shft,ps,st,sci)
% % A=load_nii(f2r); A1=A.img; 
% % AO=A.hdr.hist.originator(1:3); AD=A.hdr.dime.pixdim(2);
% % A2=nan(nVX,1,'single');
% % for i=(1:nA)+shft
% %     [I,J,K]=ind2sub(size(A1),find(A1==i));
% %     clear a; a(:,1)=I; a(:,2)=J; a(:,3)=K; 
% %     a=(a-AO)*AD;
% %     [~,~,ips]=intersect(a,ps,'rows');
% %     st1=st(ips); u=unique(st1); [~,iu,~]=intersect(u,sci); u(iu)=[];
% %     ix=[]; for j=1:length(u), ix=[ix; find(st1==u(j))]; end; ips(ix)=[];
% %     A2(ips)=i-shft;
% % end
% % end