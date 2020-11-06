
%% Variance reduction in timeseries after regressing QPPs
%%
clear; clc; p1='../';
p2p=dir([p1 'Params_*.mat']); p2p=p2p.name; 
load([p1 p2p],'p2B','p2S','p2u','d2SA','d2SAplt','d2SAcft',...
    'nsbj','nscn','nt','nT','nX','nP','PLc','PL','nVX');
p2s=[d2SA 'SA_dv.mat'];
load([p1 p2B],'B');
load([p1 p2u 'myGlssr.mat'],'iXG'); 
addpath(genpath([p1 p2u]));

%%
% dv=zeros(nX,nP,nsbj,nscn,'single');
% for is=1:nsbj
%     is
%     TT=cell(nP,1); CT=zeros(nP,nT,'single');
% for ip=1:nP
%     load(['.' p2S{is,ip}],'QPP','C'); TT{ip}=QPP(:,PLc); CT(ip,:)=C;
% for iscn=1:nscn
%     b=B{is,iscn}; c=CT(:,(1:nt)+(iscn-1)*nt); 
%     t=PL:nt; 
%     for ix=1:nX
%         r=zeros(nt-PL+1,ip,'single');
%         for ir=1:ip, r(:,ir)=conv(c(ir,:),TT{ir}(ix,:),'valid')'; end
%         beta=(r'*r)\r'*b(ix,t)';
%         dv(ix,ip,is,iscn)=std(b(ix,t)-(r*beta)',1)^2;
%     end
% end
% end
% end
% save(p2s,'dv');

%%
load(p2s,'dv');
mdv=mean(mean(dv,4),3);

figure; plot(0:nP,[ones(nX,1) mdv],'b.-'); 
hold on; plot(0:nP,[1 mean(mdv)],'ro-','linewidth',2); 
title('Variance of parcels \color{red}(average)'); xlabel('regressed QPPs');
axis([0 nP 0 1]); xticks(1:nP); axis square; grid on; box on;
saveas(gcf,[d2SAplt 'SA_dv.png']); close

DV=nan(nVX,nP,'single'); 
for ip=1:nP, for ix=1:nX, DV(iXG{ix},ip)=mdv(ix,ip); end; end
e=ft_read_cifti('empty.dtseries.nii'); e.time=1:nP; e.hdr.dim(6)=nP; 
e.dtseries=DV; ft_write_cifti([d2SAcft 'dv_QPP1-' num2str(nP) 'R'],e,...
    'parameter','dtseries');

