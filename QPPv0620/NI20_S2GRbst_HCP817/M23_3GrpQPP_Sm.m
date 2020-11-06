
%% Obtaining the summary maps within Group QPPs
%%
% The following steps are performed here:
% (2) Finding time-of-peak/-dip (te) of GrpQPP's timecourses, as a fine
% summary of activity within QPP
% > Active vertices/voxels are first found as those corresponding to
% timecourses with root sum square (rss) above rss-threshold & peak or dip
% magnitude above magnitude-threshold (see SA2). Times of peak/dip belong
% to active vertices/voxels with peaks above threshold
% > Timecourses at areas like LPCC & their te are found as refs for timing
% > Histograms of te per 7 brain regions are found (see [p2u 'myHCPcft.m'])
% (3) Clustering GrpQPP's timecourses (as a coarse summary of activity
% within QPP) that involves correlating all pairs of timecourses,
% hierarchical clustering, keeping first few largest clusters & reordering
% the clusters so that 1st cluster is default mode network etc (see below)
% (4) Finding timecourses per cluster per region & their times of peak/dip,
% finding number of timecourses per cluster per existing parcels per region

%% 
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'p2V','nvx','nP','nSD','nsd','PLe','PL','PLc','pth',...
    'avth','apth','nrgn','irgn','irf','nrf','tshclst','nclst','clstth',...
    'PR','nPR','ivx'); 
load(p2V,'QPPv');

%% Finding active vertices/voxels, times-of-peak/dip & refs for timing
%%
QPPv_te=nan(nvx,2,nP,nSD,'single'); QPPv_io=cell(nP,nSD);
QPPv_hte=zeros(PL,nrgn,2,nP,nSD,'single'); % histogram of tp/region
QPPv_rf=zeros(PLe,nrf,nP,nSD,'single'); QPPv_terf=nan(nrf,2,nP,nSD,'single');
fprintf('GrpQPP-te\n');
for ip=1:nP
for isd=1:nsd(ip)

    T=QPPv{ip,isd}(:,PLc);
    p=sqrt(sum(T.^2,2));
    io=single(find(p>pth)); clear p
    
    te=nan(nvx,2,'single'); T=double(T); warning off; tic;
    for iv=io'
        if ~rem(iv,10000), fprintf('QPP%dSD%d-te-vx%d\n',ip,isd,iv); end
        for ie=1:2
            [~,t]=findpeaks((-1)^(ie+1)*T(iv,:),'MinPeakHeight',avth,...
                'MinPeakDistance',PL-2);
            if any(t), te(iv,ie)=t; else, te(iv,ie)=0; end
        end
    end; fprintf('QPP%dSD%d-te %dsec\n\n',ip,isd,round(toc)); 
    warning on; clear T; 
    
    io=find(sum(te,2)>0); QPPv_io{ip,isd}=io;
    te(te==0)=nan; QPPv_te(:,:,ip,isd)=te;
    
    %%
    for ig=1:nrgn
    for ie=1:2
        QPPv_hte(:,ig,ie,ip,isd)=hist(te(irgn{ig},ie),1:PL);
    end
    end; clear te
    
    %%
    for ir=1:nrf
        rf=mean(QPPv{ip,isd}(irf{ir},:)); QPPv_rf(:,ir,ip,isd)=rf;
        rf=double(rf(PLc)); warning off;
        for ie=1:2
            [~,t]=findpeaks((-1)^(ie+1)*rf,'MinPeakHeight',apth,...
                'MinPeakDistance',PL-2);
            if any(t), QPPv_terf(ir,ie,ip,isd)=t; end
        end; warning on
    end; clear rf

end
end
save(p2V,'QPPv_io','QPPv_te','QPPv_hte','QPPv_rf','QPPv_terf','-append'); 
clear QPPv_te QPPv_hte QPPv_rf QPPv_terf

%% Clustering
%% Correlating pairs of timecourses "needs RAM & takes time"
QPPv_Lnkg=cell(nP,nSD); tsh=tshclst; clear tshclst; 
fprintf('GrpQPP-clst\n'); 
for ip=1:nP
for isd=1:nsd(ip)
    io=QPPv_io{ip,isd};
    T=QPPv{ip}(io,:)'; % >> transposed
    no=length(io); clear ipp; 
   
    Tj=zscore(T(PLc,:));
    Ti=zeros(2*tsh+1,PL,no,'single'); Ti(1+tsh,:,:)=Tj;
    for ish=1:tsh
        a=[T(tsh+1-ish+1:end,:); zeros(tsh+1-ish,no,'single')]; 
        Ti(ish,:,:)=zscore(a(PLc,:));  
        a=[zeros(ish,no,'single'); T(1:end-ish,:)];
        Ti(tsh+1+ish,:,:)=zscore(a(PLc,:));
    end; clear a T

    cij=zeros(no*(no-1)/2,1,'single'); tic; %% >> RAM&time consuming part
    for i=1:no-1
        if ~rem(i,10000), fprintf('clst-QPP%dSD%d-vx%d\n',ip,isd,i); end
        c0=Ti(:,:,i)*Tj(:,i+1:end)/PL;
        [c,imx]=max(abs(c0));
        for j=1:length(imx), c(j)=c(j)*sign(c0(imx(j),j)); end
        cij(1+(i-1)*(no-i/2):i*(no-(i+1)/2))=c;
    end; fprintf('clst-QPP%dSD%d %dsec\n\n',ip,isd,round(toc));
    clear c0 c Ti Tj; 
    QPPv_Lnkg{ip,isd}=linkage(1-cij','average'); clear cij; 
end
end
save(p2V,'QPPv_Lnkg','-append');

%% Clustering
fprintf('clst-QPPs-rest\n');
QPPv_clst=nan(nvx,nP,nSD,'single'); 
QPPv_clst_te=nan(2,nclst,nP,nSD,'single'); 
for ip=1:nP
for isd=1:nsd(ip)
    
    %% clustering & getting indices for sorting based on size
    I=cluster(QPPv_Lnkg{ip,isd},'Cutoff',clstth,'Criterion','Distance');
    [~,ih]=sort(hist(I,1:max(I)),'descend');
    
    %% sorting clusters based on size & keeping first nclst 
    nclst1=min(nclst,max(I));
    a=nan(size(I)); for i=1:nclst1, a(I==ih(i))=i; end
    I=nan(nvx,1,'single'); I(QPPv_io{ip,isd})=a; clear ih a

    %% to reorder clusters such that we have DMN 1st, most-anticorrelated-
    %% with-DMN 2nd & rest of clusters with ascending times of peak & dip
    %%
    T=zeros(nclst1,PLe); % average of timecourses per cluster
    t=zeros(2,nclst1); % times of peak (tp) & dip (td)
    for i=1:nclst1
        y=mean(QPPv{ip}(I==i,:)); T(i,:)=y; 
        [~,t(1,i)]=max(y(PLc)); [~,t(2,i)]=min(y(PLc)); 
    end
    iM(1)=mode(I(irf{1})); % getting current cluster index of DMN
    if isnan(iM(1)), iM(1)=mode(I(irf{2})); end % just to avoid error
    [~,iM(2)]=min(corr(T(iM(1),PLc)',T(:,PLc)')); % ~ most anticorr-w-DMN

    [~,i2]=sort(t(2,:)); % getting current cluster indices of ascending td
    [~,i1]=sort(t(1,i2)); ia=i2(i1); % ~ with 1st ascnd. tp then ascnd. td
    clear i1 i2
    
    %% combining current indices such that 1)DMN 2)~DMN 3-end)ascend tp&td
    ib=zeros(1,2); for i=1:2, ib(i)=find(ia==iM(i)); end
    ib=sort(ib); ibr=ib(1)+1:ib(2)-1; nb=length(ibr);
    iM(3:3+nb-1)=ia(ibr);
    iM(3+nb:nclst1)=[ia(ib(2)+1:end) ia(1:ib(1)-1)]; clear ia ib ibr nb

    %% reordering
    M=I; for i=1:nclst1, M(I==iM(i))=i; end; QPPv_clst(:,ip,isd)=M; 
    
    %%
    ic=1:nclst1; QPPv_clst_te(:,ic,ip,isd)=t(:,iM); clear T t I M iM ic 
    
end
end
save(p2V,'QPPv_clst','QPPv_clst_te','-append');

%% Timecourses per cluster per region & their times of peak/dip & 
%% number of timecourses per cluster per existing parcels of each region
%%
c=cell(nP,nSD); c(:)={cell(nrgn,nclst)}; 
QPPv2clst=c; QPPv2clst_te=c; 
nvx2clst=nan(nrgn,nclst,nP,nSD,'single');
nvx2clst2prcl=c; clear c
for ip=1:nP
for isd=1:nsd(ip)
for ig=1:nrgn
    I=nan(nvx,1,'single'); I(irgn{ig})=QPPv_clst(irgn{ig},ip,isd);
    R=PR{ig}(ivx); nR=nPR(ig); 
    for ic=1:nclst
        io=find(I==ic); nvx2clst2prcl{ip,isd}{ig,ic}=nan(nR,1,'single');
        if any(io)
            QPPv2clst{ip,isd}{ig,ic}=QPPv{ip,isd}(io,:);
            QPPv2clst_te{ip,isd}{ig,ic}=QPPv_te(io,:,ip,isd);
            nvx2clst(ig,ic,ip,isd)=length(io);
            nvx2clst2prcl{ip,isd}{ig,ic}=hist(R(io),1:nR);
        end
    end
end
end
end
save(p2V,'QPPv2clst','QPPv2clst_te','nvx2clst','nvx2clst2prcl','-append');
