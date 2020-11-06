
%% 
%%
clear; clc; close all; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2V','p2u','ivx','nvx','nP','nclst','d2SAplt');
load([p1 p2V],'QPPv_io','QPPv_clst','QPPv_te');
addpath(genpath([p1 p2u])); isd=1;

%%
e=ft_read_cifti('empty.dtseries.nii');
es1=e.brainstructure(ivx); en=e.brainstructurelabel';
ex1=e.pos(ivx,1);

exB=ex1(es1==7|es1==12|es1==13);
xc=median(exB);
emn=min(exB); emx=max(exB); ebn=emn:emx;
h=hist(exB,ebn); hmx=ceil(1.25*max(h)/50)*50;
% th=[quantile(exB(exB<xc),0.5) quantile(exB(exB>xc),0.5)];
th=[xc-4 xc+4];

% figure; hist(exB,ebn); axis([emn-1 emx+1 0 hmx]); hold on; 
% plot([xc xc],[0 hmx],'k:');
% for i=[-2 2], plot([xc+i xc+i],[0 hmx],'r'); end
% for i=th, plot([i i],[0 hmx],'g'); end

%% Left vs Right
% I=find(es1==7); I(ex1(I)>=xc-2 )=[]; iB{1}=[I; find(es1==12)];
% I=find(es1==7); I(ex1(I)<=xc+2)=[]; iB{2}=[I; find(es1==13)];

%% Lateral vs Medial
I=find(es1==7|es1==12|es1==13); 
iB={I(ex1(I)<=th(1)); I(ex1(I)>th(1)&ex1(I)<th(2)); I(ex1(I)>=th(2))};

%%
n=length(iB);
nB=zeros(nP,n); nBI=zeros(nclst,n,nP); tpB=cell(nP,n); mtpB=zeros(nP,n);
for ip=1:nP
    io=QPPv_io{ip,isd}; 
    I=QPPv_clst(:,ip,isd); 
    tp=QPPv_te(:,1,ip,isd);
    for i=1:n
        iio=intersect(io,iB{i}); nB(ip,i)=length(iio);
        nBI(:,i,ip)=hist(I(iio),1:nclst); 
        t=tp(iio); t(isnan(t))=[]; tpB{i,ip}=t; 
        mtpB(ip,i)=median(t);
    end
end; clear io I tp t iio
