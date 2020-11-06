
%% Significant activation/deactivation within QPP
%% Part2Q) QPPs vs magnitude threshold
%%
clear; clc; close all; p2='../../'; p1='../';
p2s='SA2.mat'; pth=0.1; ath=0.025; % simplified conclusion from prev. parts
p2p=dir([p2 'Params_*.mat']); p2p=[p2 p2p.name];
load(p2p,'p2V','d2SA','nvx','nP','nSD','nsd','PLc','d2SAplt',...
    'PL','irgn','nrgn','nmrgn1');
p2V=[p2 p2V]; d2SA=[p1 d2SA(3:end)]; d2SAplt=[p1 d2SAplt(3:end)];
load(p2V,'QPPv');
ss=get(0,'Screensize'); lw=1.5; fs=20; 
set(0,'DefaultAxesTitleFontWeight','normal');

%%
AG=cell(nrgn,nP,nSD); amx=zeros(nP,nSD); 
for ip=1:nP
for isd=1:nsd(ip)
    T=QPPv{ip,isd}(:,PLc); p=sqrt(sum(T.^2,2)); T(p<pth,:)=nan; 
    a=T(:); a(isnan(a))=[]; amx(ip,isd)=max(abs(a));
    for ig=1:nrgn
        a=T(irgn{ig},:); a=a(:); a(isnan(a))=[]; AG{ig,ip,isd}=a;
    end
end
end; clear a p T

amx=ceil(max(amx(:))*100)/100; abn=-amx:0.01:amx; lbn=length(abn);
h=zeros(lbn,nrgn,nP,nSD); 
for ip=1:nP
for isd=1:nsd(ip)
for ig=1:nrgn
    h(:,ig,ip,isd)=hist(AG{ig,ip,isd},abn);
end
end
end

%%
isd=1; for i=1:nP, l{i}=['QPP' num2str(i)]; end
figure; s=get(gcf,'Position'); 
s2=s; s2([2 4])=ss([2 4]); set(gcf,'Position',s2);
for ig=1:nrgn
    hf=h(:,ig,:,isd); hmx=round(max(hf(:))/100)*100;
    subplot(nrgn,1,ig), hold on; 
    for ip=1:nP, plot(abn,h(:,ig,ip,isd),'linewidth',lw); end 
    for j=[-1 1], z=j*ath; plot([z z],[0 hmx],'k--','linewidth',lw); end
    axis([-amx amx 0 hmx]); box on; set(gca,'fontsize',fs);
    if ig~=nrgn, xticklabels([]); end; yticklabels([]); ylabel(nmrgn1{ig});
end; xlabel('amplitude of ~20s timecourse');
subplot(nrgn,1,1), legend(l,'fontsize',fs-5)
saveas(gcf,[d2SAplt 'SA2_2q_QPP-Amp-vs-Thsh.png']); close
