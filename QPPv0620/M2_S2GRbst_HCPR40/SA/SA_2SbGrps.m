
%% Random division to 2 subgroups & compare (replicablity)
%%
clear; clc; p1='../';
p2p=dir([p1 'Params_*.mat']); p2p=p2p.name; 
load([p1 p2p],'p2B','p2O','p2qppf','p2u','d2SA','d2SAplt',...
    'nsbj','nP','tsh','PLc','nitrg','cth','PLh','G2Y','nY');
nrep=25; p2s=[d2SA 'SA_2SbGrp.mat'];
load([p1 p2B],'B');
load([p1 p2O],'QPPsa'); isd=1;
addpath([p1 p2qppf]); addpath([p2 p2u]);
set(0,'DefaultAxesTitleFontWeight','normal');

%%
n=nsbj; nh=n/2; c=zeros(nrep,nP,'single'); cw=c;
iu=find(triu(ones(nY),1)); lu=length(iu);
for irep=1:nrep
    s=randperm(n)'; IS={s(1:nh);s(nh+1:n)};
for ip=1:nP
    fprintf('QPP%d-rep%d\n',ip,irep); 
    QPPh=cell(2,1); CWQPPh=zeros(2,lu);
for ih=1:2
    [QPPh{ih},~,a]=QPPsfave(QPPsa(IS{ih},ip,isd),B(IS{ih},:),...
        tsh,PLc,nitrg,cth{ip}(2),PLh,G2Y);
    CWQPPh(ih,:)=a(iu);
end
    c(irep,ip)=Tcomp(QPPh{1},QPPh{2},PLc,tsh);
    cw(irep,ip)=corr(CWQPPh(1,:)',CWQPPh(2,:)');
end
end
save(p2s,'c');

%%
c=abs(c); 
mc=myfshr(c,1); mc=round(mc*100)/100;
mcw=myfshr(cw,1); mcw=round(mcw*100)/100;
figure; 
for i=1:nP
    subplot(2,nP,i), hist(c(:,i),0:0.1:1); title(mc(i)); 
    subplot(2,nP,i+nP), hist(cw(:,i),0:0.1:1); title(mcw(i)); 
end
subplot(2,nP,1), ylabel('corr btwn QPPs')
subplot(2,nP,1+nP), ylabel('corr btwn FCwQPPs')
saveas(gcf,[d2SAplt 'SA_2SbGrps.png']); close
