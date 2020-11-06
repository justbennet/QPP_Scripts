
%% Preprocessing: averaging over GM&WM&CSF-regressed scans (QPPs) vs 
%% WM&CSF-regressed scans (i.e., effect of GSR)
%%
clear; clc; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2prep','p2O','p2qppf','p2u','d2SAplt',...
    'nsbj','nscn','nX','nt','nT','nP','nSD','nsd','PL','PLh','PLc');
load([p1 p2prep],'p2BWCR'); a=strfind(p2prep,'/'); 
p2B=[p1 p2prep(1:a(end)) p2BWCR]; 
load([p1 p2O],'QPPs','TMXs','QPPsa','TMXsa'); 
addpath(genpath([p1 p2u])); addpath([p1 p2qppf]);
set(0,'DefaultAxesTitleFontWeight','normal'); fs=20;

%%
fprintf('Loading Data'); tic; 
load(p2B,'B'); fprintf(' %ds\n',round(toc));

c=nan(nP,nSD+1,nsbj);
for is=1:nsbj
    fprintf('Sbj%d\n',is);
    D=zeros(nX,nT,'single');
    for i=1:nscn, D(:,(i-1)*nt+(1:nt))=B{is,i}; end
for ip=1:nP
    P=QPPs{is,ip}(:,PLc);
    T=Tbld(D,TMXs{is,ip},PL,PLh,1); T=T(:,PLc);   
    c(ip,1,is)=corr(P(:),T(:));
for isd=1:nsd(ip)
    tmx=TMXsa{is,ip,isd}; n=length(tmx);
    if n
        P=QPPsa{is,ip,isd}(:,PLc);
        T=Tbld(D,tmx,PL,PLh,1); T=T(:,PLc); 
        c(ip,isd+1,is)=corr(P(:),T(:));
    end
end
end
end
[mc,~,mcn,mcx]=myfshr(c,3); [~,n]=myst(c(:),2);

%%
c1=squeeze(c(:,1,:))'; n=floor(min(c1(:))*100)/100; hmx=350; 
figure; s=get(gcf,'position'); s(3)=1.75*s(3); set(gcf,'position',s)
for ip=1:nP
    subplot(1,nP,ip), hist(c1(:,ip),n:0.01:1); axis square
    m=median(c1(:,ip)); m=round(m,3); 
    title(['QPP' num2str(ip)]); % xlabel('correlation'); 
    set(gca,'fontsize',fs); xlim([n-0.01 1.01]); ylim([0 hmx])
    text(n+0.02,hmx*0.9,['median:' num2str(m)],'fontsize',fs);
end