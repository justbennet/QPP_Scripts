
%% Preprocessing: averaging over non-preprocessed scans
%%
clear; clc; nIS=2; p1='../'; p2s='SA02.mat';
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'d2SA','nsbj'); 
IS=cell(nIS,1); p=IS; r=floor(nsbj/nIS);
for i=1:nIS, IS{i}=(1:r)+(i-1)*r; p{i}=[d2SA 'f' num2str(i) '_' p2s]; end
IS{nIS}=[IS{nIS} IS{nIS}(end)+1:nsbj];  

%% running fXPreprocon on servers (~nsbj/2 per node)
% i=1; fXPreproc(IS{i},p{i}); 
% i=2; fXPreproc(IS{i},p{i}); 

%% 
load(p2p,'p2O','p2u','d2SAcft','ivx','nVX','nP','nsd','PLe'); 
load([p1 p2O],'NMX'); addpath(genpath([p1 p2u]));

%% Combining & writing to cifti
load(p{1},'QPPvxp');
for i=2:nIS
    a=load(p{i},'QPPvxp'); 
    for ip=1:nP, for isd=1:nsd(ip)
            QPPvxp{ip,isd}=QPPvxp{ip,isd}+a.QPPvxp{ip,isd}; end; end
end

for ip=1:nP, for isd=1:nsd(ip)
    QPPvxp{ip,isd}=QPPvxp{ip,isd}/sum(NMX(ip,isd,:)); end; end
save([d2SA p2s],'QPPvxp');

e=ft_read_cifti('empty.dtseries.nii'); e.time=1:PLe; e.hdr.dim(6)=PLe; 
for ip=1:nP
for isd=1:nsd(ip)
    T=nan(nVX,PLe,'single'); T(ivx,:)=QPPvxp{ip,isd}; e.dtseries=T;
    ft_write_cifti([d2SAcft 'QPPxp' num2str(ip) 'SD' num2str(isd)],e,...
        'parameter','dtseries'); 
end
end
