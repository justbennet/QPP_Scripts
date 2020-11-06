
%% Correlating QPP fine summary maps with the corresponding FCG at cortex
%%
clear; clc; close all; p1='../'; isd=1; 
irf=[1 1 2]; % indices of "rf" as seed for correlation map within QPP
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2V','nP','PLc','nvxc','p2u','nVX','ivx','d2SAcft');
p2V=[p1 p2V]; p2u=[p1 p2u]; addpath(genpath(p2u));
load(p2V,'QPPv_te','QPPv','QPPv_io','QPPv_rf');

a=ft_read_cifti('hcp.embed.dscalar.nii');
b=a.x100307_tfmri_motor_level2_cue_hp200_s2; g(:,1)=b(ivx(1:nvxc));
b=a.x100307_tfmri_motor_level2_lf_hp200_s2; g(:,3)=b(ivx(1:nvxc));
b=a.x100307_tfmri_motor_level2_lh_hp200_s2; g(:,2)=b(ivx(1:nvxc));
clear a b

%%
cTG=zeros(nP,2); pTG=cTG; CWT=nan(nvxc,nP,'single');
for ip=1:nP
    tp1=QPPv_te(1:nvxc,1,ip,isd); % tp within QPP
    ix=find(isnan(tp1)); tp1(ix)=[];
    g1=g(:,ip); g1(ix)=[]; 
    [cTG(ip,1),pTG(ip,1)]=corr(tp1,g1); % corr btwn QPP's-tp & FCG
    
    io=QPPv_io{ip,isd}; io(io>nvxc)=[];
    T=QPPv{ip,isd}(io,PLc)'; rf=QPPv_rf(PLc,irf(ip),ip,isd);
    c1=corr(T,rf); CWT(io,ip)=c1; % correlatiom within QPP with a ref seed
    g1=g(io,ip);
    [cTG(ip,2),pTG(ip,2)]=corr(c1,g1); % corr btwn QPP's-corr & FCG
end

%% writing QPP's-corr into cifti
M=nan(nVX,nP);
for ip=1:nP, M(ivx(1:nvxc),ip)=CWT(:,ip); end
e=ft_read_cifti('empty.dtseries.nii');
e.time=1:nP; e.hdr.dim(6)=nP; e.dtseries=M;
ft_write_cifti([d2SAcft 'Corr_within_QPP_Ref'],e,'parameter','dtseries');
