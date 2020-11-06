 
%%  Plots in grayordinate
%%
clear; clc;
p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nVX','ivx','nvx','p2u','p2qppf','p2V','d2cft','d2plt','nP',...
    'PLe','PL','PLc','nrgn','irgn','nmrgn1'); 
addpath(genpath(p2u)); addpath(p2qppf);
e=ft_read_cifti('empty.dtseries.nii');
load(p2V,'QPPv','QPPv_io','QPPv_tp');

%% Writing QPPs into cifti
e.time=1:PLe; e.hdr.dim(6)=PLe; 
for i=1:nP
    T=nan(nVX,PLe,'single'); T(ivx,:)=QPPv{i}; e.dtseries=T;
    ft_write_cifti([d2cft 'QPP' num2str(i)],e,'parameter','dtseries');
end; clear T

%% Number & percentage of active vertices/voxels overall & per brain region
noa=zeros(nP,1,'single'); roa=noa;
no=zeros(nrgn,nP,'single'); ro=no; lg=zeros(nrgn,1,'single');
for ip=1:nP
    io=QPPv_io{ip};
    noa(ip)=length(io); 
    roa(ip)=noa(ip)/nvx*100;
    I=zeros(nvx,1,'single'); I(io)=1;
    for ig=1:nrgn
        J=zeros(nvx,1,'single'); J(irgn{ig})=I(irgn{ig});
        no(ig,ip)=length(find(J));
        lg(ig)=length(irgn{ig});
        ro(ig,ip)=no(ig,ip)/lg(ig)*100;
    end
end; clear I J io

a1=[no lg]; a1=[a1; sum(a1)]; % for copy-pasting into tables
a2=round(ro); a2=[a2; round(roa')];

%% Writing time-of-peak maps into ciftis
e.time=1:nP; e.hdr.dim(6)=nP; 
TP=nan(nVX,nP,'single'); TP(ivx,:)=QPPv_tp; e.dtseries=TP; clear TP
ft_write_cifti([d2cft 'tp_QPPs'],e,'parameter','dtseries');

