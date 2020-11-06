 
%%  Plots in grayordinate
%%
clear; clc; close all; 
flg=0; % plotting timecourses/cluster, flg=0 (50/clst & fast), flg=1 (all)

%%
p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nVX','ivx','nvx','p2u','p2qppf','p2V','d2cft','d2plt','nP',...
    'nsd','nSD','PLe','PL','PLc','nrgn','nmrgn1','irgn','nrf','clstth',...
    'nclst','clstclr','YLB','nY','nmPR','nPR'); nth=15;
addpath(genpath(p2u)); addpath(p2qppf);

clh=colormap('jet'); close; % bin colors for histogram of tp/region (htp)
clh=clh(round(linspace(1,size(clh,1),PL)),:); lw=4; % width of bins in htp
irp={[1 2 3];[1 3 4];[2 3];[2 3];[2 3]}; % indices of "rf" as refs in htp
mrp={'-';'--';'-.';':'}; % line-type when plotting aboveset indices
alm=repmat(0.75*[-1 1],nP,1); % amp when plotting timecourses/cluster
cly=clstclr([1 2 6 4 7 9 10],:); % colors for Yeo's 7 RSNs
ss=get(0,'Screensize'); figure; ss2=get(gcf,'Position'); close;
set(0,'DefaultAxesTitleFontWeight','normal');
e=ft_read_cifti('empty.dtseries.nii');

load(p2V,'QPPv','QPPv_io','QPPv_tp','QPPv_tprf','QPPv_htp','QPPv_tm',...
    'pttst_tp2','mstp2','QPPv_Lnkg','QPPv_clst','QPPv_rf','QPPv2clst',...
    'nvx2clst','pttst_tp','mstp','ntp','QPPv2y','msytp','pttsty_tp',...
    'nytp','nvx2clst2prcl');

%% Writing QPPs into cifti
e.time=1:PLe; e.hdr.dim(6)=PLe; 
for i=1:nP
    T=nan(nVX,PLe,'single'); T(ivx,:)=QPPv{i}; e.dtseries=T;
    ft_write_cifti([d2cft 'QPP' num2str(i)],e,'parameter','dtseries');
end; clear T

%% Number & percentage of active vertices/voxels overall & per brain region
noa=zeros(nP,1,'single'); roa=noa;
no=zeros(nrgn,nP,'single'); ro=no; lg=zeros(nrgn,1,'single');
for ip=1:nP
    io=QPPv_io{ip};
    noa(ip)=length(io); 
    roa(ip)=noa(ip)/nvx*100;
    I=zeros(nvx,1,'single'); I(io)=1;
    for ig=1:nrgn
        J=zeros(nvx,1,'single'); J(irgn{ig})=I(irgn{ig});
        no(ig,ip)=length(find(J));
        lg(ig)=length(irgn{ig});
        ro(ig,ip)=no(ig,ip)/lg(ig)*100;
    end
end; clear I J io

% a1=[no lg]; a1=[a1; sum(a1)]; % for copy-pasting into tables
% a2=round(ro); a2=[a2; round(roa')];

%% Writing time-of-peak maps into ciftis
e.time=1:nP; e.hdr.dim(6)=nP; 
TP=nan(nVX,nP,'single'); TP(ivx,:)=QPPv_tp; e.dtseries=TP; clear TP
ft_write_cifti([d2cft 'tp_QPPs'],e,'parameter','dtseries');

%% Histogram of time-of-peak/-dip per brain regions
figure; set(gcf,'Position',ss);
for ip=1:nP
for ig=1:nrgn
    h=QPPv_htp(:,ig,ip); hmx=max(h);
    subplot(nrgn,nP,(ig-1)*nP+ip),hold on; 
    for k=1:PL, plot([k k],[0 h(k)],'linewidth',lw,'color',clh(k,:)); end 
    for k=1:length(irp{ip}), ia=irp{ip}(k); a=QPPv_tprf(ia); 
        plot([a a],[0 hmx],['k' mrp{ia}],'linewidth',1); end
    axis([1 PL 0 hmx]); xticks(5:5:PL); yticks([]); box on; 
    if ig~=nrgn, xticklabels([]); end; if ig==nrgn, xlabel('timepoint'); end
    if ig==1, title(['QPP' num2str(ip)]); end
    if ip==1, ylabel(nmrgn1{ig}); end
end
end; saveas(gcf,[d2plt '32_tp_hist.png']); close; clear a h hmx

%% Significant differences of time of peak of 2nd mode between regions
s=ss2; s(1)=1; s(3)=ss(3); figure; set(gcf,'Position',s);
for ip=1:nP
    p=max(pttst_tp2(:,:,ip,:),[],4);
    subplot(1,nP,ip),im=imagesc(p); axis square; colorbar;
    set(im,'AlphaData',~isnan(p)); title(['QPP' num2str(ip)]);
    xticks(1:nrgn); yticks(1:nrgn); xticklabels(nmrgn1); yticklabels(nmrgn1);
end; saveas(gcf,[d2plt '32_tp_hist_pval.png']); close

% a1=squeeze(mstp2(:,1,:)); % for table

%% Size of clusters of timecourses
figure; set(gcf,'Position',ss);
for ip=1:nP
    I=cluster(QPPv_Lnkg{ip},'Cutoff',clstth(ip),'Criterion','Distance');
    [h,ih]=sort(hist(I,1:max(I)),'descend');
    subplot(2,nP,ip), stem(h); 
    title(['QPP' num2str(ip)]); if ip==1, ylabel('cluster size'); end
    subplot(2,nP,nP+ip), stem(log10(h)); 
    xlim([0 20]); if ip==1, ylabel('log_{10} (cluster size)'); end    
end; saveas(gcf,[d2plt '33_clst_size.png']); close

%% Plotting timecourses per cluster & writing cluster maps into ciftis
MM=nan(nVX,nP); s=[1 ss2(2) ss(3) ss2(4)/2];
for ip=1:nP
    M=QPPv_clst(:,ip);
    figure; set(gcf,'Position',s);
    for i=1:nclst
        T=QPPv{ip}(M==i,:); 
        if ~flg, iT=randperm(size(T,1)); iT=iT(1:50); T=T(iT,:); end
        if size(T,1)==PLe, T(end,:)=[]; end
        subplot(1,nclst,i), plot(1:PLe,T,'color',clstclr(i,:));
        hold on; plot(QPPv_rf(:,1,ip),'k');
        PLTUSD2(PLc,alm(ip,:),~(i-1),~(i-1)); title(i)
    end
    f2s=[d2plt '34_clst_QPP' num2str(ip) '_1.png'];
    saveas(gcf,f2s); close;
    MM(ivx,ip)=M;
end
e.dtseries=MM; e.time=1:nP; e.hdr.dim(6)=nP; 
ft_write_cifti([d2cft 'clst_QPPs' ],e,'parameter','dtseries');

%% Size of clusters of timecourses per region
n=nvx2clst; n(n<nth)=nan; n(n==0)=nan; 
% ip=1; a1=n(:,:,ip)'; % for tables

%% Plotting timecourses per cluster per region
for ip=1:nP
    r=QPPv_rf(:,1,ip); r=r/max(abs(r));
    figure; set(gcf,'Position',ss);
for ig=1:nrgn
     a=0; for ic=1:nclst, A=max(abs(QPPv2clst{ip}{ig,ic}(:)));
        if any(A), a=max(a,A); end; end; a=ceil(a*20)/20;
for ic=1:nclst
    T=QPPv2clst{ip}{ig,ic}; n=nvx2clst(ig,ic,ip);
    if n<nth, n=nan; end
    if ~isnan(n)
        if ~flg, iT=randperm(n); iT=iT(1:min(n,200)); T=T(iT,:); end
        subplot(nrgn,nclst,ic+(ig-1)*nclst), hold on; 
        if size(T,1)==PLe, T(end,:)=[]; end
        plot(1:PLe,T,'color',clstclr(ic,:)); plot(r*a,'k');
        PLTUSD2(PLc,[-a a],0,~(ic-1)); 
        if ic==1, ylabel(nmrgn1{ig}), end
        if ig==1, title(ic); end
    elseif ic==1, subplot(nrgn,nclst,ic+(ig-1)*nclst), 
        PLTUSD2(PLc,[-a a],0,1); ylabel(nmrgn1{ig}); 
    end
end
end; saveas(gcf,[d2plt '34_clst_QPP' num2str(ip) '_3.png']); close
end

%% Significant difference of time of peak/dip between clusters & regions 
nrc=nrgn*nclst; il=nclst/2:nclst:nrc; l=cell(nrc,1); 
for i=1:length(il), l{il(i)}=nmrgn1{i}; end
for ip=1:nP
    figure; set(gcf,'Position',ss); 
    p=max(pttst_tp(:,:,ip,:),[],4); 
    im=imagesc(p); set(im,'AlphaData',~isnan(p)); colorbar; axis square; 
    xticks(1:nrc); yticks(1:nrc); xticklabels(l); yticklabels(l); 
    hold on; grid on;
    for i=nclst+0.5:nclst:nrc-nclst+0.5, plot([i i],[0 nrc+1],'k'); end
    for i=nclst+0.5:nclst:nrc-nclst+0.5, plot([0 nrc+1],[i i],'k'); end
    saveas(gcf,[d2plt '34_clst_QPP' num2str(ip) '_tp_pval.png']); close
end

for i=1:2, m=mstp(:,:,i,:); m(isnan(ntp))=nan; mstp(:,:,i,:)=m; end
% ip=1; a1=mstp(:,:,1,ip)'; % for tables

%% Plotting timecourses per cortical RSNs
s=ss2; s(3)=s(3)*1.25; xtk=[PLc(1) PL PLc(end)]; xtkl=xtk-PLc(1)+1; 
for ip=1:nP
    m=QPPv2y(:,:,1,ip); e=QPPv2y(:,:,2,ip);
    figure; set(gcf,'Position',s); hold on; 
    for iy=1:nY
        errorbar(m(iy,:),e(iy,:),'color',cly(iy,:)); 
    end
    xticks(xtk); xticklabels(xtkl); box on; grid on; legend(YLB);
    a=ceil(max(abs(m(:)))*100)/100; yticks([-a 0 a]); axis([1 PLe -a a]); 
    saveas(gcf,[d2plt '35_QPPv2RSN_' num2str(ip) '.png']); close
end

%% Times of peak for cortical RSNs
ip=1; a1=msytp(:,1,ip); % for tables

s=[1 ss2(2) ss(3) ss2(4)]; figure; set(gcf,'Position',s);
for ip=1:nP
    p=max(pttsty_tp(:,:,ip,:),[],4);
    subplot(1,nP,ip), im=imagesc(p); axis square; colorbar
    set(im,'AlphaData',~isnan(p)); title(['QPP' num2str(ip)]);
    xticks(1:nY); yticks(1:nY); xticklabels(YLB); yticklabels(YLB); 
end; saveas(gcf,[d2plt '35_QPPv2RSN_tp_pval.png']); close

%% Correlation between cortical RSNs
% figure; s=ss2; s(1)=1; s(3)=ss(3); set(gcf,'Position',s); id=find(eye(nY));
% for ip=1:nP
%     c=corr(QPPv2y(:,PLc,1,ip)'); c(id)=nan;
%     subplot(1,nP,ip),PLTFCY7(c,[-1 1],0,1,1); title(['QPP' num2str(ip)])
% end
% saveas(gcf,[d2plt '35_QPPv2RSN_zcorr.png']); close

% ip=1; a1=nytp(:,ip); % for tables

