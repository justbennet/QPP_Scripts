
%% Random division to 2 subgroups & compare QPPs
%%
clear; clc; p1='../'; nrep=50; p2s='SA1.mat'; isd=1; 
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name]; 
load(p2p,'p2B','p2O','p2qppf','p2u','d2SA','d2SAplt',...
    'nsbj','nP','tsh','PLc','nitrg','cth','PLh');
addpath([p1 p2qppf]); addpath([p1 p2u]);
set(0,'DefaultAxesTitleFontWeight','normal');

%%
fprintf('Loading Data'); tic; 
load([p1 p2B],'B'); fprintf(' %ds\n',round(toc)); 
load([p1 p2O],'QPPsa');

n=nsbj; nh=round(n/2); 
c=zeros(nrep,nP,'single');
for irep=1:nrep
    s=randperm(n)'; IS={s(1:nh);s(nh+1:n)};
for ip=1:nP
    fprintf('QPP%d-rep%d',ip,irep); tic;
    QPPh=cell(2,1);
for ih=1:2
    QPPh{ih}=QPPsfave(QPPsa(IS{ih},ip,isd),B(IS{ih},:),...
        tsh,PLc,nitrg,cth{ip}(2),PLh);
end
    c(irep,ip)=Tcomp(QPPh{1},QPPh{2},PLc,tsh); fprintf(' %ds\n',round(toc));
end
end
save([d2SA p2s],'c');

%%
% load([d2SA p2s],'c');
[ix,ipx]=ind2sub(size(c),find(c<0)); iph=hist(ipx,1:nP);
c1=c; c1(ix,:)=[];
mc=myfshr(c1,1); mc=round(mc*100)/100; n=floor(min(c1(:))*100)/100;
figure; for i=1:nP, subplot(1,nP,i), hist(c1(:,i),n:0.01:1); 
    xlim([n-0.01 1.01]); title(mc(i)); axis square; end
saveas(gcf,[d2SAplt 'SA1_2SbGrps.png']); close
