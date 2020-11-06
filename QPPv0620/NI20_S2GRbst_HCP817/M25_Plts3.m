 
%%  Plots in grayordinate
%%
clear; clc; close all; 
flg=0; % plotting timecourses/cluster, flg=0 (50/clst & fast), flg=1 (all)

%%
p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nVX','ivx','nvx','p2u','p2qppf','p2V','d2cft','d2plt','nP',...
    'nsd','nSD','PLe','PL','PLc','nrgn','nmrgn1','irgn','nrf','clstth',...
    'nclst','clstclr','YLB','nY','ayth','nmPR','nPR'); nth=15;
addpath(genpath(p2u)); addpath(p2qppf);

clh=colormap('jet'); close; % bin colors for histogram of tp/region (htp)
clh=clh(round(linspace(1,size(clh,1),PL)),:); lw=7; % width of bins in htp
irp={[1 2 3];[1 3 4];[2 3]}; % indices of cell array "rf" as refs in htp
mrp={'-';'--';'-.';':'}; % line-type when plotting aboveset indices
nmte='pdz'; % for naming time-of-peak/-dip/-zero-crossing
alm=repmat(0.75*[-1 1],nP,1); % amp when plotting timecourses/cluster
cly=clstclr([1 2 6 4 7 9 10],:); % colors for Yeo's 7 RSNs
ss=get(0,'Screensize'); figure; ss2=get(gcf,'Position'); close;
set(0,'DefaultAxesTitleFontWeight','normal'); fs=20; % fontsize
e=ft_read_cifti('empty.dtseries.nii');

load(p2V,'QPPv','QPPv_io','QPPv_te','QPPv_terf','QPPv_hte','QPPv_tm',...
    'pttst_tp2','mstp2','QPPv_Lnkg','QPPv_clst','QPPv_rf','QPPv2clst',...
    'nvx2clst','pttst_te','mste','nte','QPPv2y','msyte','pttsty_te',...
    'nyte','nvx2clst2prcl');

%% Writing QPPs into cifti
e.time=1:PLe; e.hdr.dim(6)=PLe; 
for i=1:nP
for j=1:nsd(i)
    T=nan(nVX,PLe,'single'); T(ivx,:)=QPPv{i,j}; e.dtseries=T;
    ft_write_cifti([d2cft 'QPP' num2str(i) 'SD' num2str(j)],e,...
        'parameter','dtseries');
end
end; clear T

%% Number & percentage of active vertices/voxels overall & per brain region
noa=zeros(nP,nSD,'single'); roa=noa;
no=zeros(nrgn,nP,nSD,'single'); ro=no; lg=zeros(nrgn,1,'single');
for ip=1:nP
for isd=1:nsd(ip)
    io=QPPv_io{ip,isd};
    noa(ip,isd)=length(io); 
    roa(ip,isd)=noa(ip,isd)/nvx*100;
    I=zeros(nvx,1,'single'); I(io)=1;
    for ig=1:nrgn
        J=zeros(nvx,1,'single'); J(irgn{ig})=I(irgn{ig});
        no(ig,ip,isd)=length(find(J));
        lg(ig)=length(irgn{ig});
        ro(ig,ip,isd)=no(ig,ip,isd)/lg(ig)*100;
    end
end
end; clear I J io
rlm=[floor(min(ro(:))/10)*10 ceil(max(ro(:))/10)*10];

figure; ro=round(ro);
for i=1:nSD
    subplot(1,nSD,i), imagesc(ro(:,:,i),rlm); colormap('jet'); colorbar; 
    axis image; xticks(1:nP); yticks(1:nrgn); yticklabels(nmrgn1); 
    xlabel('QPP'); title(['SD' num2str(i)]);
end; subplot(1,nSD,1), ylabel('region'); 
saveas(gcf,[d2plt '30_3_Percentage-of-active-vox.png']); close;

% a1=[no(:,:,1) lg]; a1=[a1; sum(a1)]; % for copy-pasting into tables
% a2=round(ro(:,:,1)); a2=[a2; round(roa(:,1)')];

%% Writing time-of-peak/-dip maps into ciftis
n=sum(nsd)*2; TE=nan(nVX,n,'single'); cnt=1;
for i=1:nP
for j=1:nsd(i)
    for k=1:2, TE(ivx,cnt)=QPPv_te(:,k,i,j); cnt=cnt+1; end
end
end
e.time=1:n; e.hdr.dim(6)=n; e.dtseries=TE; clear TE
ft_write_cifti([d2cft 'te_QPPs'],e,'parameter','dtseries');

%% Histogram of time-of-peak/-dip per brain regions
for isd=1:nSD
for ie=1:2
    figure; set(gcf,'Position',ss);
for ip=1:nP
    a=QPPv_hte(:,:,ie,ip,isd);
if sum(a(:))
for ig=1:nrgn
    h=QPPv_hte(:,ig,ie,ip,isd); hmx=max(h);
    subplot(nrgn,nP,(ig-1)*nP+ip),hold on; 
    for k=1:PL, plot([k k],[0 h(k)],'linewidth',lw,'color',clh(k,:)); end 
    for k=1:length(irp{ip}), ia=irp{ip}(k); a=QPPv_terf(ia,ie,ip,isd); 
        plot([a a],[0 hmx],['k' mrp{ia}],'linewidth',1.5); end
    if ie==1, a=QPPv_tm(ip,isd); 
        plot([a a],[0 hmx],'color',clh(a,:),'linewidth',1.5); end
    axis([1 PL 0 hmx]); xticks(5:5:PL); yticks([]); box on; 
    if ig~=nrgn, xticklabels([]); end; set(gca,'Fontsize',fs);
    if ip==1, ylabel(nmrgn1{ig}); end
end; clear a h mx
end
end
    f2s=[d2plt '32_' num2str(ie) '_t' nmte(ie) '_hist_SD' num2str(isd) '.png'];
    print(gcf,f2s,'-dpng','-r1000'); close
end
end

%% Significant differences of time of peak of 2nd mode between regions
s=ss2; s(1)=1; s(3)=ss(3); 
for isd=1:2
    figure; set(gcf,'Position',s);
for ip=1:nP
    p=max(pttst_tp2(:,:,ip,isd,:),[],5);
    subplot(1,nP,ip),im=imagesc(p); axis square;
    set(im,'AlphaData',~isnan(p));
    xticks(1:nrgn); yticks(1:nrgn); xticklabels(nmrgn1); yticklabels(nmrgn1);
    set(gca,'fontsize',fs)
end
    saveas(gcf,[d2plt '32_3_tp2_pval_SD' num2str(isd) '_4ppr.png']); close
end

% isd=1;
% a1=squeeze(mstp2(:,1,:,isd)); % for tables
% a2=round(squeeze(mstp2(:,2,:,isd))*100)/100;

%% Size of clusters of timecourses
for l=1:2  
    figure; set(gcf,'Position',ss);
for ip=1:nP
for isd=1:nsd(ip)
    I=cluster(QPPv_Lnkg{ip,isd},'Cutoff',clstth,'Criterion','Distance');
    [h,ih]=sort(hist(I,1:max(I)),'descend'); 
    subplot(2*nSD,nP,(isd-1)*nP+ip), stem(h); 
    title(['QPP' num2str(ip) 'SD' num2str(isd)]);
    if ip==1 && isd==1, ylabel('cluster size'); end
    if l==2, xlim([0 20]); end
    subplot(2*nSD,nP,(isd+1)*nP+ip), stem(log10(h)); 
    title(['QPP' num2str(ip) 'SD' num2str(isd)]);
    if ip==1 && isd==1, ylabel('log_{10} (cluster size)'); end
    if l==2, xlim([0 20]); end 
end
end; saveas(gcf,[d2plt '33_clst_size_' num2str(l) '.png']); close
end

%% Plotting timecourses per cluster & writing cluster maps into ciftis
n=sum(nsd); MM=nan(nVX,n); cnt=1;
s=[1 ss2(2) ss(3) ss2(4)/2];
for ip=1:nP
for isd=1:nsd(ip)
    fprintf('QPP%dSD%d\n',ip,isd);
    M=QPPv_clst(:,ip,isd);
    figure; set(gcf,'Position',s);
    for i=1:nclst
        T=QPPv{ip,isd}(M==i,:); 
        if ~flg, iT=randperm(size(T,1)); iT=iT(1:50); T=T(iT,:); end
        if size(T,1)==PLe, T(end,:)=[]; end
        subplot(1,nclst,i), plot(1:PLe,T,'color',clstclr(i,:));
        hold on; plot(QPPv_rf(:,1,ip,isd),'k');
        PLTUSD2(PLc,alm(ip,:),0,0); set(gca,'Fontsize',fs);
    end
    f2s=[d2plt '34_clst_QPP' num2str(ip) 'SD' num2str(isd) '_1.png'];
    saveas(gcf,f2s); close; % print(gcf,f2s,'-dpng','-r500');    
    MM(ivx,cnt)=M; cnt=cnt+1;
end
end
e.dtseries=MM; e.time=1:n; e.hdr.dim(6)=n; 
ft_write_cifti([d2cft 'clst_QPPs' ],e,'parameter','dtseries');

%% Size of clusters of timecourses per region
n=nvx2clst; n(n<nth)=nan; n(n==0)=nan; 
% ip=1; isd=1; a1=n(:,:,ip,isd)'; % for tables

n=round(log10(n),1); nmx=max(n(:)); 
figure; set(gcf,'Position',ss);
for ip=1:nP
for isd=1:nsd(ip)
    n1=n(:,:,ip,isd);
	subplot(nSD,nP,ip+(isd-1)*nP),im=imagesc(n1,[0 nmx]); colormap('jet'); 
    axis image; title(['QPP' num2str(ip) 'SD' num2str(isd)]); 
    xticks(1:nclst); yticks(1:nrgn);
    yticklabels([]); yticklabels(nmrgn1); if ip==1, ylabel('region'); end
    if isd==nSD, xlabel('QPP cluster'); end
    set(im,'AlphaData',~isnan(n1));
end
end
subplot(nSD,nP,2), title({'log_{10}(cluster size per region)';'QPP2SD1'}); 
subplot(nSD,nP,nP); a=get(gca,'Position'); 
a=[a(1)+a(3)+0.01 a(2)-a(4) 0.02 a(4)*1.75]; colorbar('position',a); 
saveas(gcf,[d2plt '33_clst_size_per-region.png']); close 

%% Plotting timecourses per cluster per region
fs1=18;
for ip=1:nP
for isd=1:nsd(ip)
    r=QPPv_rf(:,1,ip,isd); r=r/max(abs(r));
    figure; set(gcf,'Position',ss);
for ig=1:nrgn
     a=0; for ic=1:nclst, A=max(abs(QPPv2clst{ip,isd}{ig,ic}(:)));
        if any(A), a=max(a,A); end; end; a=ceil(a*20)/20;
for ic=1:nclst
    T=QPPv2clst{ip,isd}{ig,ic}; n=nvx2clst(ig,ic,ip,isd);
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
    set(gca,'Fontsize',fs1);
end
end
    subplot(nrgn,nclst,1), title('cluster1');
    saveas(gcf,[d2plt '34_clst_QPP' num2str(ip) 'SD' num2str(isd) '_3.png']); close
end
end

%% Significant difference of time of peak/dip between clusters & regions 
nrc=nrgn*nclst; il=nclst/2:nclst:nrc; l=cell(nrc,1); 
for i=1:length(il), l{il(i)}=nmrgn1{i}; end
for ip=1:nP
for isd=1:nsd(ip)
for ie=1:2
    figure; set(gcf,'Position',ss); 
    p=max(pttst_te(:,:,ie,ip,isd,:),[],6); 
    im=imagesc(p); set(im,'AlphaData',~isnan(p)); colorbar; axis square; 
    xticks(1:nrc); yticks(1:nrc); xticklabels(l); yticklabels(l); 
    hold on; grid on; set(gca,'fontsize',fs)
    for i=nclst+0.5:nclst:nrc-nclst+0.5, plot([i i],[0 nrc+1],'k'); end
    for i=nclst+0.5:nclst:nrc-nclst+0.5, plot([0 nrc+1],[i i],'k'); end
    a=get(gca,'Position');
    a=[a(1)+a(3)-0.1 a(2)+a(4)/3 0.02 a(4)/3]; colorbar('position',a);
    saveas(gcf,[d2plt '34_zclst_QPP' num2str(ip) 'SD' num2str(isd) ...
        '_pval_' num2str(ie) '_t' nmte(ie) '.png']); close
end
end
end

for i=1:2, m=mste(:,:,:,i,:,:); m(isnan(nte))=nan; mste(:,:,:,i,:,:)=m; end
% isd=1; ip=1; ie=1; a1=mste(:,:,ie,1,ip,isd)'; % for tables

s=[1 ss2(2) ss(3) ss2(4)*3/4];
for isd=1:nSD
for ip=1:nP
    figure; set(gcf,'Position',s);
for ie=1:2
    m=mste(:,:,ie,1,ip,isd);
    subplot(1,2,ie), im=imagesc(m); colormap('jet'); colorbar 
    set(im,'AlphaData',~isnan(m)); title(['QPP' num2str(ip) ' t' nmte(ie)]);
    axis image; xticks(1:nclst); yticks(1:nrgn); yticklabels(nmrgn1)
end
    saveas(gcf,[d2plt '34_zclst_QPP' num2str(ip) 'SD' num2str(isd) ...
        '_med.png']); close
end
end

%% Plotting timecourses per cortical RSNs
s=ss2; s(3)=s(3)*1.25; xtk=[PLc(1) PL PLc(end)]; xtkl=xtk-PLc(1)+1; 
for ip=1:nP
for isd=1:nsd(ip)
    m=QPPv2y(:,:,1,ip,isd); e=QPPv2y(:,:,2,ip,isd);
    figure; set(gcf,'Position',s); hold on; 
    for iy=1:nY
        errorbar(m(iy,:),e(iy,:),'color',cly(iy,:),'linewidth',2); 
    end
    xticks(xtk); xticklabels(xtkl);
    for i=[-ayth ayth], plot([1 PLe],[i i],'k'); end
    a=ceil(max(abs(m(:)))*100)/100;
    ytk=[-a -ayth ayth a]; ytkl={num2str(-a),'','',num2str(a)};
    yticks(ytk); yticklabels(ytkl);
    for i=xtk, plot([i i],[-a a],'k'); end
    axis([1 PLe -a a]); box on; set(gca,'fontsize',fs);
    saveas(gcf,[d2plt '35_QPPv2RSN_' num2str(ip) '_SD' num2str(isd) '.png']); close
end
end

%% Times of peak/dip/zero-crossing for cortical RSNs
s=ss2; s(2)=1; s(4)=ss(4); figure; set(gcf,'Position',s);
for ip=1:nP
for isd=1:nsd(ip)
    t=msyte(:,:,1,ip,isd);
    subplot(nSD,nP,ip+(isd-1)*nP), im=imagesc(t,[3 PL-2]); colormap('jet');
    title(['QPP' num2str(ip) 'SD' num2str(isd)]);
    axis image; xticks(1:3); yticks(1:nY); xticklabels({'tp','td','tz'});
    if ip==1, yticklabels(YLB); else, yticklabels([]); end
end
end
subplot(nSD,nP,nP); a=get(gca,'Position'); 
a=[a(1)+a(3)+0.01 a(2)-a(4) 0.03 a(4)*1.5]; colorbar('position',a);
saveas(gcf,[d2plt '35_QPPv2RSN_te.png']); close

% ip=1; isd=1; a1=msyte(:,[1 3],1,ip,isd); % for tables

s=ss2; s(3)=s(3)*1.25;
for ie=[1 3]
    figure; set(gcf,'Position',s);
for ip=1:nP
for isd=1:nsd(ip)
    p=max(pttsty_te(:,:,ie,ip,isd,:),[],6);
    subplot(nSD,nP,(isd-1)*nP+ip), im=imagesc(p); axis square; 
    set(im,'AlphaData',~isnan(p));
    title(['QPP' num2str(ip) ' SD' num2str(isd)]);
    xticks(1:nY); yticks(1:nY); xticklabels(YLB); yticklabels(YLB); 
end
end
    subplot(nSD,nP,nP); a=get(gca,'Position'); 
    a=[a(1)+a(3)+0.01 a(2)-a(4)/2-0.1 0.02 a(4)*1.5]; colorbar('position',a);
    saveas(gcf,[d2plt '35_QPPv2RSN_t' nmte(ie) '_pval.png']); close
end

% s=ss2; s(3)=s(3)*1.75; figure; set(gcf,'Position',s); % just for paper
% cnt=1; ttl={'peak','zero-crossing'}; ip=1; isd=1;
% for ie=[1 3] 
%     p=max(pttsty_te(:,:,ie,ip,isd,:),[],6);
%     subplot(1,2,cnt), im=imagesc(p); axis square; 
%     set(im,'AlphaData',~isnan(p)); title(['time of ' ttl{cnt}]);
%     xticks(1:nY); yticks(1:nY); xticklabels(YLB); yticklabels(YLB); 
%     set(gca,'fontsize',fs); cnt=cnt+1;
% end
% saveas(gcf,[d2plt '35_QPPv2RSN_tpz_pval_QPP' num2str(ip) ...
%     'SD' num2str(isd) '_4ppr.png']); close

%% Correlation between cortical RSNs
id=find(eye(nY));
figure; s=ss2; s(3)=s(3)*1.5; set(gcf,'Position',s); 
for ip=1:nP
for isd=1:nsd(ip)
    c=corr(QPPv2y(:,PLc,1,ip,isd)'); c(id)=nan;
    subplot(nSD,nP,ip+(isd-1)*nP),PLTFCY7(c,[-1 1],0,1,1);
    title(['QPP' num2str(ip) 'SD' num2str(isd)])
end
end
subplot(nSD,nP,nP); a=get(gca,'Position');
a=[a(1)+a(3)+0.01 a(2)-a(4)/2 0.02 a(4)]; colorbar('position',a);
saveas(gcf,[d2plt '35_QPPv2RSN_z_corr.png']); close

% ip=1; isd=1; c=corr(QPPv2y(:,PLc,1,ip,isd)'); % just for paper
% id=find(eye(nY)); c(id)=nan; 
% figure; PLTFCY7(c,[-1 1],1,3,1); set(gca,'fontsize',fs)
% saveas(gcf,[d2plt '35_QPPv2RSN_z_corr_QPP' 'SD' '_4ppr.png']); close

% a3=nyte(:,3,1); % for tables

%% Number of timecourses per cluster per existing parcels of each region
iPR=[0; cumsum(nPR)]; 
C=nan(iPR(end),nclst,nP,nSD,'single');
for ip=1:nP
for isd=1:nsd(ip)
for ic=1:nclst
for ig=1:nrgn    
    C(iPR(ig)+1:iPR(ig+1),ic,ip,isd)=nvx2clst2prcl{ip,isd}{ig,ic}';
end
end
end
end
C(C<nth)=nan; C(C==0)=nan;
a1=C(:,:,1,1); a2=C(:,:,2,1); a3=C(:,:,3,1); % for tables

cmx=round(log10(max(C(:))),1); 
for isd=1:nSD
    figure; set(gcf,'Position',ss);
for ig=1:nrgn
for ip=1:nP
    nR=nPR(ig); a=zeros(nR,nclst,'single');
    for ic=1:nclst, a(:,ic)=nvx2clst2prcl{ip,isd}{ig,ic}; end
    a(a<nth)=nan; a(a==0)=nan; a=round(log10(a),1);
    subplot(nP,nrgn,ig+(ip-1)*nrgn),im=imagesc(a,[1 cmx]); colormap('jet');
    set(im,'AlphaData',~isnan(a));  
    if ig==1,ylabel(['QPP' num2str(ip)]);end
    if ip==1,title(nmrgn1{ig}); end
    xticks(1:nclst); yticks(1:nR); yticklabels(nmPR{ig});
end
end
    subplot(nP,nrgn,2*nrgn); a=get(gca,'Position'); 
    a=[a(1)+a(3)+0.01 a(2)-a(4)/2 0.02 a(4)*2]; colorbar('position',a);
    saveas(gcf,[d2plt '36_QPPv2Prcls_SD' num2str(isd) '.png']); close
end
