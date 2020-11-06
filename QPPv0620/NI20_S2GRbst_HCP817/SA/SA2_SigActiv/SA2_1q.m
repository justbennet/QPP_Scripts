
%% Significant activation/deactivation within QPP
%% Part1Q) Checking QPPs vs power threshold
%%
clear; clc; close all; p2='../../'; p1='../';
p2s='SA2.mat'; pth=0.1; % simplified conclusion from part1n
p2p=dir([p2 'Params_*.mat']); p2p=[p2 p2p.name];
load(p2p,'p2V','nvx','nP','nSD','nsd','PL','PLc','d2SA',...
    'd2SAplt','irgn','nrgn','nmrgn1');
p2V=[p2 p2V]; d2SA=[p1 d2SA(3:end)]; d2SAplt=[p1 d2SAplt(3:end)];
load(p2V,'QPPv'); 
ss=get(0,'Screensize'); ssw=ss; ssw(2)=ss(4)/4; ssw(4)=ss(4)/2;
set(0,'DefaultAxesTitleFontWeight','normal'); lw=1.5; fs=20;

%%
p=zeros(nvx,nP,nSD,'single');
r=zeros(nrgn,nP,nSD,'single');
for ip=1:nP
for isd=1:nsd(ip)
    T=QPPv{ip,isd}(:,PLc); 
    p(:,ip,isd)=sqrt(sum(T.^2,2));
    for ig=1:nrgn
        r(ig,ip,isd)=length(find(p(irgn{ig},ip,isd)>=pth))/length(irgn{ig});
    end
end
end

pmx=ceil(max(p(:))*10)/10; pbn=0:0.01:pmx; lbn=length(pbn);
r=round(r*100*10)/10;

h=zeros(lbn,nrgn,nP,nSD);
for ip=1:nP
for isd=1:nsd(ip)
    for ig=1:nrgn
        h(:,ig,ip,isd)=hist(p(irgn{ig},ip,isd),pbn);
    end
end
end

%%
isd=1; for i=1:nP, l{i}=['QPP' num2str(i)]; end
figure; s=get(gcf,'Position'); 
s2=s; s2([2 4])=ss([2 4]); set(gcf,'Position',s2);
for ig=1:nrgn
    hf=h(:,ig,:,isd); hmx=max(hf(:));
    subplot(nrgn,1,ig), hold on 
    for ip=1:nP, plot(pbn,h(:,ig,ip,isd),'Linewidth',lw); end
    plot([pth pth],[0 hmx],'k--','Linewidth',lw);
    axis([0 pbn(end) 0 hmx]); yticklabels([]); ylabel(nmrgn1{ig}); 
    if ig~=nrgn, xticklabels([]); end; set(gca,'fontsize',fs); box on
end; xlabel('root sum square of ~20s timecourse');
subplot(nrgn,1,1), legend(l,'fontsize',fs-5)
saveas(gcf,[d2SAplt 'SA2_1q_1_QPP-Pwr-vs-Thsh.png']); close

figure; imagesc(r(:,:,isd)); colormap('jet'); colorbar; axis image
xticks(1:nP); xlabel('QPP'); yticks(1:nrgn); yticklabels(nmrgn1);
set(gca,'fontsize',fs); 
saveas(gcf,[d2SAplt 'SA2_1q_2_Percent-Vox-above-PwrThsh.png']); close
