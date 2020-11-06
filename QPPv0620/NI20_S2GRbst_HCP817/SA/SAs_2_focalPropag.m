
%% Significance of difference of time of peak/dip for focal Propagations 
%%
clear; clc; close all; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2V','p2u','irgn','nvx','PLe','PLc','clstclr','p2qppf','d2SAplt');
nth=15; nrep=50;  
load('myGlssr.mat','ixG'); 
load('myHCPcft.mat','irgnlr','nmrgnlr1');

load([p1 p2V],'QPPv_te','QPPv_clst','QPPv','QPPv_rf');
addpath(genpath([p1 p2u])); addpath([p1 p2qppf]);

set(0,'DefaultAxesTitleFontWeight','normal'); flg=1; fs=20; 
ss=get(0,'screensize'); figure; ss2=get(gcf,'Position'); close; 

%%
ip=1; isd=1; 

%% Temporal lobe, hippocampus & amygdala 
%%
icR=[1 7:10]; nc=length(icR); clr=clstclr(icR,:); 
N11R=[175,125,129,128,130,176,123,107]; N11L=N11R+180;
N14R=[137,133,177,132,136,134,172,131,135]; N14L=N14R+180;
iTLR=sort([N11R N14R]); iTLL=sort([N11L N14L]); clear N11R N11L N14R N14L
I=[]; for i=1:length(iTLR), I=cat(1,I,ixG{iTLR(i)}); end; ivxTLR=sort(I);
I=[]; for i=1:length(iTLL), I=cat(1,I,ixG{iTLL(i)}); end; ivxTLL=sort(I);
iHAL=sort([irgnlr{7};irgnlr{9}]); iHAR=sort([irgnlr{8};irgnlr{10}]);
iA={ivxTLL; ivxTLR; iHAL; iHAR}; nA=length(iA); ncA=nc*nA;  
alph=0.01/(ncA*(ncA-1)/2);
nmA={'TL-L','TL-R','H&A-L','H&A-R'}; f2s='TL-HA'; 
s2=[ss2(1) 1 1.5*ss2(3) ss(4)]; ie2p=2;

%% V1
%%
% icR=1:6; nc=length(icR); clr=clstclr(icR,:);
% iV1=[1 181];
% I=[]; for i=1:2, I=cat(1,I,ixG{iV1(i)}); end; ivxV1=sort(I);
% iA={ivxV1}; nA=1; ncA=nc*nA; 
% alph=1e-6; % 0.01/(ncA*(ncA-1)/2);
% nmA={'V1'}; f2s='V1'; s2=ss2; s2(3)=2*ss2(3); s2(4)=ss2(4)/2; ie2p=1;

%% Common part for TL-H&A and V1 (uncomment whichever is to test & plot)
%%
c=cell(nA,nc); QPPv2clst_te=c; QPPv2clst=c; clear c
nvx2clst=zeros(nA,nc,'single'); mte=zeros(nA,nc,2,'single');
for ia=1:nA
    I=nan(nvx,1,'single'); I(iA{ia})=QPPv_clst(iA{ia},ip,isd);
    for ic=1:nc
        io=find(I==icR(ic)); te=QPPv_te(io,:,ip,isd);
        QPPv2clst_te{ia,ic}=te;
        QPPv2clst{ia,ic}=QPPv{ip,isd}(io,:);
        nvx2clst(ia,ic)=length(io);
        mte(ia,ic,:)=median(te,'omitnan');
    end
end

pttstf_te=nan(ncA,ncA,2,nrep,'single');
for irep=1:nrep
for ie=1:2
for ia1=1:nA
for ic1=1:nc
    t=QPPv2clst_te{ia1,ic1}(:,ie);
    t(isnan(t))=[]; l1=length(t);
    if l1>=nth
        te1=t; 
        for ia2=1:nA
        for ic2=1:nc
            t=QPPv2clst_te{ia2,ic2}(:,ie);
            t(isnan(t))=[]; l2=length(t);
            if l2>=nth
                te2=t;
                l=min(l1,l2);
                if l==l1, t1=te1; t2=te2(randperm(l2)); t2=t2(1:l);
                else, t2=te2; t1=te1(randperm(l1)); t1=t1(1:l); 
                end 
                [h,p]=ttest(t1,t2,'Alpha',alph);
                i=(ia2-1)*nc+ic2; j=(ia1-1)*nc+ic1;
                if ~h, pttstf_te(i,j,ie,irep)=1;
                else, pttstf_te(i,j,ie,irep)=p;
                end
            end
        end    
        end
    end
end
end
end
end

%% 
il=round(nc/2):nc:ncA; l=cell(ncA,1); 
for i=1:length(il), l{il(i)}=nmA{i}; end
ttl={'Time of peak','Time of dip'};
figure; ie=ie2p; p=mode(pttstf_te(:,:,ie,:),4);
im=imagesc(p); set(im,'AlphaData',~isnan(p)); axis square;
xticks(1:ncA); yticks(1:ncA); 
if ~strcmp(f2s(1),'V'), xticklabels(l); yticklabels(l); end
title(ttl{ie}); hold on; set(gca,'fontsize',fs);
for i=nc+0.5:nc:ncA-nc+0.5, plot([i i],[0 ncA+1],'k'); end
for i=nc+0.5:nc:ncA-nc+0.5, plot([0 ncA+1],[i i],'k'); end
saveas(gcf,[d2SAplt 'SAs_focalPropag_' f2s '_pval.png']); close
a1=mte(:,:,ie2p);

r=QPPv_rf(:,1,ip,isd); r=r/max(abs(r));
figure; set(gcf,'Position',s2)
for ia=1:nA
     a=0; for ic=1:nc, A=max(abs(QPPv2clst{ia,ic}(:)));
        if any(A), a=max(a,A); end; end; a=ceil(a*20)/20;
for ic=1:nc
    T=QPPv2clst{ia,ic}; n=nvx2clst(ia,ic);
    if n<nth, n=nan; end
    if ~isnan(n)
        if ~flg, iT=randperm(n); iT=iT(1:min(n,20)); T=T(iT,:); end
        subplot(nA,nc,ic+(ia-1)*nc), hold on;
        if size(T,1)==PLe, T(end,:)=[]; end
        plot(1:PLe,T,'color',clr(ic,:)); plot(r*a,'k');
        PLTUSD2(PLc,[-a a],0,~(ic-1)); 
        if ic==1, ylabel(nmA{ia}), end; if ia==1, title(icR(ic)); end
    elseif ic==1, subplot(nA,nc,ic+(ia-1)*nc), 
        PLTUSD2(PLc,[-a a],0,1); ylabel(nmA{ia}); 
    end
    set(gca,'Fontsize',fs);
end
end
saveas(gcf,[d2SAplt 'SAs_focalPropag_' f2s '.png']); close
