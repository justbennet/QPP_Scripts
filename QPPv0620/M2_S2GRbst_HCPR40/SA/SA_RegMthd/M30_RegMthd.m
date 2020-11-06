
%%
clear; clc; close all; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'p2O','p2qppf','p2u','nP','nSD','nsd','nsbj','PLc','tsh','d2plt');
addpath(p2qppf); addpath(p2u);

load(p2O,'QPPsa','QPPs','QPP','pFC'); 
a=load('../../S2GQPP_Rbst_HCPR40.mat','QPPsa','QPPs','QPP','pFC'); 
Psa=a.QPPsa; Ps=a.QPPs; P=a.QPP; pfc=a.pFC; pfc(nP+2:end,:)=[];

%%
c=nan(nP,nSD+1,nsbj);
for is=1:nsbj
for ip=1:nP
for isd=1:nsd(ip)
    T1=QPPsa{is,ip,isd}; T2=Psa{is,ip,isd};
    if ~isempty(T1) && ~isempty(T2), c(ip,isd,is)=Tcomp(T1,T2,PLc,tsh); end
end
    c(ip,isd+1,is)=Tcomp(QPPs{is,ip},Ps{is,ip},PLc,tsh);
end
end

cg=nan(nP,nSD);
for ip=1:nP
for isd=1:nsd(ip)
    T1=QPP{ip,isd}; T2=P{ip,isd};
    if ~isempty(T1) && ~isempty(T2)
        cg(ip,isd)=Tcomp(T1,T2,PLc,tsh); 
    end
end
end

%%
figure; cbn=0:0.1:1;
for ip=1:nP
for isd=1:nsd(ip)+1
    c1=abs(c(ip,isd,:));
    subplot(nSD+1,nP,ip+(isd-1)*nP), hist(c1(:),cbn); xlim([0 1.05])
    title(['QPP' num2str(ip) 'SD' num2str(isd)]);
    if isd==nsd(ip)+1, title(['QPP' num2str(ip)]); end
end
end
saveas(gcf,[d2plt 'RegMthds_Sbjlvl.png']); close

figure; imagesc(abs(cg),[0 1]); axis image; colorbar; 
xticks(1:nSD); yticks(1:nP); xlabel('SD'); ylabel('QPP'); 
saveas(gcf,[d2plt 'RegMthds_Grplvl.png']); close

figure; clr='brk'; ymx=55;
subplot(1,2,1), hold on;for i=1:3, plot(pfc(:,i),[clr(i) 'o-']); end
axis([0 nP+2 0 ymx]); xticks(1:nP+1); xticklabels(0:nP); grid on; box on
title('dFC scan-wise regression');
subplot(1,2,2), hold on;for i=1:3, plot(pFC(:,i),[clr(i) 'o-']); end
axis([0 nP+2 0 ymx]); xticks(1:nP+1); xticklabels(0:nP); grid on; box on
title('dFC segment-wise regression');
saveas(gcf,[d2plt 'RegMthds_dFC.png']); close

%%
