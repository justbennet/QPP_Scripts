
%% Grp_Fst (original method by Majeed et el 2011) vs Sbj2Grp_Rbst
%%
clear; clc; close all; p1='../';
p2fst='../../M1_GrpFst_HCPR40/GrpQPP_Fst_HCPR40.mat';

p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load([p1 p2p],'p2O','p2qppf','p2u','nP','nSD','nsd','PLc','tsh',...
    'G2Y','nY','nX','PL','d2SAplt');
load([p1 p2O],'QPP'); QPP_rbst=QPP; clear QPP
load(p2fst,'QPP','QPPa'); QPP_fst=QPP; QPP_fst_phadj=QPPa; clear QPP QPPa
addpath([p1 p2qppf]); addpath([p1 p2u]);

%%
c=nan(nP,nSD); ccw=c; iu=find(triu(ones(nY),1));
for ip=1:nP
for isd=1:nsd(ip)
    T{1}=QPP_rbst{ip,isd}; T{2}=QPP_fst_phadj{ip,isd};
    c(ip,isd)=Tcomp(T{1},T{2},PLc,tsh);
    cw=zeros(2,length(iu));
    for j=1:2
        P=zeros(nY,PL);
        for i=1:nY, P(i,:)=mean(T{j}(G2Y==i,PLc),1); end; c1=corr(P'); 
        cw(j,:)=c1(iu);
    end
    ccw(ip,isd)=corr(cw(1,:)',cw(2,:)');
end
end

c2=nan(nP,1); ccw2=c2; isd=1;
for ip=1:nP
    T{1}=QPP_rbst{ip,isd}; T{2}=QPP_fst{ip};
    c2(ip)=Tcomp(T{1},T{2},PLc,tsh);
    cw=zeros(2,length(iu));
    for j=1:2
        P=zeros(nY,PL);
        for i=1:nY, P(i,:)=mean(T{j}(G2Y==i,PLc),1); end; c1=corr(P'); 
        cw(j,:)=c1(iu);
    end
    ccw2(ip)=corr(cw(1,:)',cw(2,:)');
end

ccw=[ccw2 ccw]; [~,n,x]=myst(ccw(:),2);
figure; PLTC(ccw,[n x],1,{'','phadjSD1','~2'},1:nP); 
title('Grp QPP obtn by Fst vs Rbst mthd','fontweight','normal'); 
ylabel('QPP'); xlabel('QPP Fst');
saveas(gcf,[d2SAplt 'SA_Rbst_vs_Fst.png']); close
