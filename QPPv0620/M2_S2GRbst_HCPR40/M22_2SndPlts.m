
%% Secondary plots
%%
clear; clc; close all; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nsbj','p2O','p2S','p2G','p2qppf','p2u','d2plt',...
    'nP','nSD','nsd','PLc','PL','tsh','cthph','nitr','cth',...
    'nt','nT','nX'); 
load(p2O,'QPPs','QPPsa','QPP','TMX');
addpath(p2qppf); addpath(p2u);

nsp1=5; nsp2=8; % #subplots set for nsbj=40
alm=repmat([-1 1],nP,1);
for i=4:nP, alm(i,:)=0.5*alm(i,:); end
set(0,'DefaultAxesTitleFontWeight','normal');
ss=get(0,'Screensize'); ssw=ss; ssw(2)=ss(4)/4; ssw(4)=ss(4)/2;

%% 5 SbjQPPs
%% 5.1 QPPi
for ip=1:nP
    figure; set(gcf,'Position',ss);
for is=1:nsbj
    subplot(nsp1,nsp2,is), PLTT(QPPs{is,ip},PLc,alm(ip,:));
end; saveas(gcf,[d2plt '51_SbjQPP' num2str(ip) '.png']); close;
end

%% 5.2 QPPi vs QPPj
c=nan(nP,nP,nsbj);
for is=1:nsbj
for ip1=1:nP-1
for ip2=ip1+1:nP
	c(ip1,ip2,is)=abs(Tcomp(QPPs{is,ip1},QPPs{is,ip2},PLc,tsh));
end
end
end
[~,~,x]=myst(c(:),2); 
[mc,~,mcmn,mcmx]=myfshr(c,3); 

figure; 
subplot(1,2,1), hist(c(:),0:0.01:x); 
title('QPPi vs QPPj, all pairs all sbj'); xlabel('corr'); 
xlim([0 x]); axis square
subplot(1,2,2), PLTC(mc,[mcmn mcmx],1);
title('med corr btwn pairs'); xlabel('QPP'); ylabel('QPP');
saveas(gcf,[d2plt '52_QPPi_vs_QPPj.png']); close

%% 5.3 Average over original vs residual scans for QPP2+ (QPP vs QPP0)
c=nan(nP-1,nSD+1,nsbj); 
for is=1:nsbj
for ip=2:nP
    load(p2S{is,ip},'QPP0','QPPa0')
	c(ip-1,1,is)=Tcomp(QPPs{is,ip},QPP0,PLc,tsh);
for isd=1:nsd(ip)
    T=QPPa0{isd}; 
    if any(T), c(ip-1,isd+1,is)=Tcomp(QPPsa{is,ip,isd},T,PLc,tsh); end
end
end
end; [~,n]=myst(c(:),2); [mc,~,mcmn,mcmx]=myfshr(c,3); 

figure;
subplot(1,2,1), hist(c(:),n:0.01:1); 
title({'QPP built by ave orig scans vs'; '~ resid scans, all QPPs all sbj'});
xlabel('corr'); xlim([n 1]); axis square
subplot(1,2,2), PLTC(mc,[mcmn mcmx],1,{'','phadjSD1','~2'},2:nP); 
title('med corr btwn pairs'); ylabel('QPP');
saveas(gcf,[d2plt '53_orig_vs_resid-scans.png']); close

%% 5.4 Phase-adjusting (phadj): QPPi vs other templates (cT1Tj & nsim)
c=cell(nP,1); nsm=zeros(nP,nsbj);
for ip=1:nP
for is=1:nsbj
    load(p2S{is,ip},'cT1Tj','nsim'); 
    a=cT1Tj; a(a==0)=[]; c{ip}=[c{ip}; a]; nsm(ip,is)=nsim;
end
end

cbn=0:0.01:1; figure; set(gcf,'Position',ssw);
for ip=1:nP
    subplot(2,nP,ip), hist(c{ip},cbn); 
    title(['QPP' num2str(ip)]); xlim([0 1]); xticks([0 cthph 1]);
    subplot(2,nP,ip+nP), hist(nsm(ip,:));
end
subplot(2,nP,1), title('QPP1 vs other templates'); xlabel('corr'); 
subplot(2,nP,nP+1), title('# sim templates'); 
saveas(gcf,[d2plt '54_QPPi_vs_other-templates.png']); close

%% 5.5 Seeds of QPP vs phase-adjusted QPPs
ip=1; 
for isd=1%:nsd(ip)
    figure; set(gcf,'Position',ss);
for is=1:nsbj
    load(p2S{is,ip},'SD','SDa');
    subplot(nsp1,nsp2,is), hold on
    plot(SD(isd,:),'b','linewidth',1.5);
    plot(SDa(isd,:),'g','linewidth',1.5); PLTUSD(PLc,[-1 1]);
end; saveas(gcf,[d2plt '55_QPP' num2str(ip) 'vs_phadj_at_SD' ...
    num2str(isd) '.png']); close;
end

for ip=2:nP 
for isd=1%:nsd(ip)
    figure; set(gcf,'Position',ss);
for is=1:nsbj
    load(p2S{is,ip},'SD0','SDa0'); 
    subplot(nsp1,nsp2,is), hold on
    plot(SD0(isd,:),'b','linewidth',1.5); 
    plot(SDa0(isd,:),'g','linewidth',1.5); PLTUSD(PLc,[-1 1]); 
end; saveas(gcf,[d2plt '55_QPP' num2str(ip) 'vs_phadj_at_SD' ...
    num2str(isd) '.png']); close;
end
end

%% 5.6 Phase-adjusted QPPs
for ip=1:nP
for isd=1%:nsd(ip)
    figure; set(gcf,'Position',ss);
for is=1:nsbj
    T=QPPsa{is,ip,isd}; 
    if any(T), subplot(nsp1,nsp2,is), PLTT(T,PLc,alm(ip,:)); end
end; saveas(gcf,[d2plt '56_SbjQPPphadj' num2str(ip) ...
    'SD' num2str(isd) '.png']); close
end
end

%% Iterations of the QPP algorithm
ITR=zeros(nP,nsbj);
for ip=1:nP
for is=1:nsbj
    load(p2S{is,ip},'ITER'); ITR(ip,is)=ITER;
end
end; m=round(median(ITR,2));
figure; set(gcf,'Position',ssw);
for ip=1:nP
    subplot(1,nP,ip),hist(ITR(ip,:),0:nitr+1); 
    title(m(ip)); xlim([0 nitr+2]); 
end; subplot(1,nP,1), title(['QPP1, med:' num2str(m(1))]); 
saveas(gcf,[d2plt '59_iteration-QPP-algorithm.png']); close

%% 6 Goodness of regression
cth1=0.2; 
cbn=cth1-0.05:0.01:0.8; l=length(cbn); h=zeros(nP,l,2);
for is=1:nsbj
    load(p2S{is,nP},'Cr');
for ip=1:nP
    load(p2S{is,ip},'C');
    c=C; c=c(c>cth1); h(ip,:,1)=h(ip,:,1)+hist(c,cbn);
    c=Cr(ip,:); c=c(c>cth1); h(ip,:,2)=h(ip,:,2)+hist(c,cbn);
end
end

a1=['# corr vals > ' num2str(cth1) ', btwn QPP_i & orig scans']; 
a2='~ QPP_i & resid.'; 
a=cell(1,nP); for i=1:nP, a{i}=['QPP' num2str(i)]; end
figure; set(gcf,'Position',ssw);
for i=1:2
    subplot(1,2,i), plot(cbn,h(:,:,i),'.-');
    title(eval(['a' num2str(i)])); xlabel('corr'); legend(a);
    axis([cbn(1) cbn(end) 0 max(h(:))/i]);
end; saveas(gcf,[d2plt '61_Goodness-of-Regression.png']); close

figure; set(gcf,'Position',ss); is=1; load(p2S{is,nP},'Cr'); t=1:nT;
for i=1:nP 
    load(p2S{is,i},'C'); 
    subplot(nP,1,i), plot(C(t),'b'); hold on; plot(Cr(i,t),'r'); 
    xlim([t(1) t(end)]); PLTUMET(length(t),length(t)/nt,[-0.4 0.6],cth1);
end
saveas(gcf,[d2plt '62_Goodness-of-Regressiong_2.png']); close

%% 7.1 Sbj2Grp: QPP of Sbj i vs Sbj j (cTiTj)
c=cell(nP,nSD); mc=nan(nP,nSD);
for ip=1:nP
for isd=1:nsd(ip)
    load(p2G{ip,isd},'cij'); 
    cij(eye(nsbj)>0)=nan; cij=cij(:); cij(cij==0)=nan; 
    c{ip,isd}=cij; 
    [~,mc(ip,isd)]=myfshr(cij,1);
end
end

cbn=0:0.01:1; figure; set(gcf,'Position',ssw);
for ip=1:nP
for isd=1:nsd(ip)
    subplot(nSD,nP,(isd-1)*nP+ip), hist(c{ip,isd},cbn); 
    title(mc(ip,isd)); xlim([0 1]);
    if ip==1, title(['QPP1SD' num2str(isd) ', med:' num2str(mc(1,1))]); end
end
end; saveas(gcf,[d2plt '71_Sbj-i_vs_Sbj-j.png']); close

%% 7.2 Sbj2Grp: SbjQPPs phase-matched to reference sbj
isd=1;
for ip=1:nP
    load(p2G{ip,isd},'tshij','sgnij','isref'); 
    figure; set(gcf,'Position',ss);
for is=1:nsbj
    T=QPPsa{is,ip,isd};
    if any(T)
        s=tshij(is,isref);
        if s<tsh+1, s=tsh+1-s; T=[T(:,s+1:end) zeros(nX,s,'single')]; 
        elseif s>tsh+1, s=s-(tsh+1); T=[zeros(nX,s,'single') T(:,1:end-s)];
        end
        subplot(nsp1,nsp2,is), PLTT(T,PLc,alm(ip,:));
    end
    if is==isref, set(gca,'linewidth',2); end
end; saveas(gcf,[d2plt '72_QPP' num2str(ip) ...
    '_phase-matched-to-QPPref_SD' num2str(isd) '.png']); close
end

%% 7.3 Prior GrpQPP vs GrpQPP
figure; set(gcf,'Position',ss);
for ip=1:nP
for isd=1:nsd(ip)
    load(p2G{ip,isd},'Tg0');
    subplot(2*nSD,nP,ip+(isd-1)*2*nP), PLTT(Tg0,PLc,alm(ip,:),1);
    if ip==1, title(['Prior GrpTMPL1 SD' num2str(isd)]); end
    subplot(2*nSD,nP,ip+nP+(isd-1)*2*nP), PLTT(QPP{ip,isd},PLc,alm(ip,:),1)
    if ip==1, title(['GrpQPP1 SD' num2str(isd)]); end
end
end; saveas(gcf,[d2plt '73_Prior-Grp-TMPL_vs_GrpQPP.png']); close;

c=nan(nP,nSD);
for ip=1:nP
for isd=1:nsd(ip)
    load(p2G{ip,isd},'Tg0');
    c(ip,isd)=Tcomp(Tg0,QPP{ip,isd},PLc,tsh);
end
end; [~,n,x]=myst(c(:),2);
figure; PLTC(c,[n x],1);
title('Prior GrpTMPL vs GrpQPP'); ylabel('QPP'); xlabel('SD');
saveas(gcf,[d2plt '74_Prior-Grp-TMPL_vs_GrpQPP.png']); close

%% 8 GrpQPP project back to sbj
for ip=1:nP
for isd=1%:nsd(ip)
    load(p2G{ip,isd},'QPPg2s','isref'); 
    figure; set(gcf,'Position',ss); 
for is=1:nsbj
    n=length(TMX{ip,isd}{is});
    if n, subplot(nsp1,nsp2,is), PLTT(QPPg2s{is},PLc,[-1 1]); title(n); end
    if is==isref, set(gca,'linewidth',2); end
end; saveas(gcf,[d2plt '81_Grp-project-back-to-Sbj_QPP' ...
    num2str(ip) 'SD' num2str(isd) '.png']); close
end
end
