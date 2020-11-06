
%% Obtaining group QPP in grayordinate, finding active vertices/voxels 
%% within group QPP & obtaining the summary maps 
%%
% The following steps are performed here:
% (1) Building GrpQPP by averaging its contributing segments (identified
% based on parcel-space analysis) over grayordinate scans 
% > Null patterns are also built in grayordinates (by averaging random
% segments selected previously), root sum square (rss) & absolute value of
% amplitude (magnitude) of timecourses are found & 99th quantiles of all
% rss values & magnitudes are taken as rss & magnitude threshold for GrpQPP
% (2) Finding time-of-peak (tp) of GrpQPP's timecourses, as a fine
% summary of activity within QPP
% > Active vertices/voxels are first found as those corresponding to
% timecourses with rss above rss-threshold & peak or dip magnitude above
% magnitude-threshold. Times of peak belong to active vertices/voxels with
% peaks above threshold
% > Timecourses at areas like LPCC & their tp are found as refs for timing
% > Histograms of tp per 7 brain regions are found (see [p2u 'myHCPcft.m'])
% (3) Clustering GrpQPP's timecourses (as a coarse summary of activity
% within QPP) that involves correlating all pairs of timecourses,
% hierarchical clustering, keeping first few largest clusters & reordering
% the clusters so that 1st cluster is default mode network etc (see below)
% (4) Finding timecourses per cluster per region & their times of peak,
% finding number of timecourses per cluster per existing parcels per region

%% 
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nsbj','nscn','nvx','nt','nT','p2bpr','p2O','p2V','p2NV',...
    'p2qppf','p2u','nP','PLe','PL','PLh','PLc','nnll','nrgn','irgn',...
    'irf','nrf','tshclst','nclst','clstth','PR','nPR','ivx'); isd=1;
load([p2u 'myGlssr.mat'],'nG','ixG'); nX=1+nG; iX=[irf{1}; ixG];
load(p2O,'TMXa','TMX','TMXN'); addpath(p2qppf);

%% GrpQPP in grayordinate
%% Dividing segments per each subject (note the next section for reason)
TMX1=cell(nP,nsbj); nmx1=zeros(nP,1,'single'); 
TMXN1=cell(nP,nnll,nsbj); nmxN1=zeros(nP,1,'single');
for ip=1:nP
    if any(TMXa{ip,isd}), t=TMXa{ip,isd}; else, t=TMX{ip}; end 
    for is=1:nsbj
        t1=t; t1( t1<(is-1)*nT+1 | t1>is*nT )=[]; t1=t1-(is-1)*nT;
        TMX1{ip,is}=t1;
    end
    nmx1(ip)=length(t);
    for in=1:nnll
        t=TMXN{ip,in};
        for is=1:nsbj
            t1=t; t1( t1<(is-1)*nT+1 | t1>is*nT )=[]; t1=t1-(is-1)*nT;
            TMXN1{ip,in,is}=t1;
        end
    end
    nmxN1(ip)=length(TMXN{ip,1}); 
end; clear TMX TMXa TMXN t t1

%% Averaging segments across subjects
QPPv=cell(nP,1); QPPv(:)={zeros(nvx,PLe,'single')};
QPPNv=zeros(nvx,PL,nP,nnll,'single');
fprintf('GrpQPP-Vx\n');
for is=1:nsbj
    fprintf('Sbj%d ',is); tic;
    D=zeros(nvx,nT,'single'); 
    for i=1:nscn
        load(p2bpr{is,i},'bpr'); 
        D(:,(i-1)*nt+(1:nt))=zscore(bpr,[],2);
    end; clear bpr
    for ip=1:nP
        t=TMX1{ip,is};
        if ~isempty(t), QPPv{ip}=QPPv{ip}+Tbld(D,t,PL,PLh,0); end
        for in=1:nnll
            t=TMXN1{ip,in,is};
            if ~isempty(t)
                QPPNv(:,:,ip,in)=QPPNv(:,:,ip,in)+Tbld(D,t,PL,[0 0],0); 
            end
        end
    end
    fprintf('%dsec\n',round(toc));
end; clear D TMX1 TMXN1 T t
for ip=1:nP, QPPv{ip}=QPPv{ip}/nmx1(ip); end; save(p2V,'QPPv');

%% Finding root-sum-square & magnitude thresholds based on null patterns
qvN=zeros(nP,3,'single'); rssN=zeros(nvx,nP,nnll,'single'); 
QPPNp=zeros(nX,PL,nP,nnll,'single');
for ip=1:nP
    QPPNv(:,:,ip,:)=QPPNv(:,:,ip,:)/nmxN1(ip);
    for in=1:nnll
        rssN(:,ip,in)=sqrt(sum(QPPNv(:,:,ip,in).^2,2)); 
        for i=1:nX, QPPNp(i,:,ip,in)=mean(QPPNv(iX{i},:,ip,in),1); end
    end
    qvN(ip,1)=quantile(reshape(rssN(:,ip,:),[],1),0.99);
    qvN(ip,2)=quantile(reshape(abs(QPPNv(:,:,ip,:)),[],1),0.99);
    qvN(ip,3)=quantile(reshape(abs(QPPNp(:,:,ip,:)),[],1),0.99);
end
save(p2V,'qvN','-append'); % save(p2NV,'QPPNv','-v7.3'); 
clear QPPNv QPPNp rssN

%%  Finding active vertices/voxels, times of peak & references for timing
%%
QPPv_tp=nan(nvx,nP,'single'); QPPv_io=cell(nP,1); 
QPPv_htp=zeros(PL,nrgn,nP,'single'); % histogram of tp/region
QPPv_rf=zeros(PLe,nrf,nP,'single'); QPPv_tprf=nan(nrf,nP,'single'); 
fprintf('\nGrpQPP-tp\n');
for ip=1:nP
    
    T=QPPv{ip}(:,PLc);
    rss=sqrt(sum(T.^2,2));
    io=single(find(rss>qvN(ip,1))); clear rss
    
    te=nan(nvx,2,'single'); T=double(T); warning off; tic;
    for iv=io'
        if ~rem(iv,10000), fprintf('QPP%d-tp-vx%d\n',ip,iv); end
        for ie=1:2
            [~,t]=findpeaks((-1)^(ie+1)*T(iv,:),...
                'MinPeakHeight',qvN(ip,2),'MinPeakDistance',PL-2);
            if any(t), te(iv,ie)=t; else, te(iv,ie)=0; end
        end
    end; fprintf('QPP%d-tp %dsec\n\n',ip,round(toc)); 
    warning on; clear T; 
    
    io=find(sum(te,2)>0); QPPv_io{ip}=io;
    te(te==0)=nan; QPPv_tp(:,ip)=te(:,1);
    
    %%
    for i=1:nrgn, QPPv_htp(:,i,ip)=hist(te(irgn{i},1),1:PL); end; clear te
    
    %%
    for ir=1:nrf
        rf=mean(QPPv{ip}(irf{ir},:)); QPPv_rf(:,ir,ip)=rf;
        rf=double(rf(PLc)); warning off;
        [~,t]=findpeaks(rf,'MinPeakHeight',qvN(ip,3),'MinPeakDistance',PL-2);
        if any(t), QPPv_tprf(ir,ip)=t; end
        warning on
    end; clear rf
      
end
save(p2V,'QPPv_io','QPPv_tp','QPPv_htp','QPPv_rf','QPPv_tprf','-append'); 
clear QPPv_htp QPPv_rf QPPv_tprf; 

%% Clustering
%% Correlating pairs of timecourses "needs RAM & takes time"
QPPv_Lnkg=cell(nP,1); tsh=tshclst; clear tshclst
fprintf('GrpQPP-clst\n'); 
for ip=1:nP    
    io=QPPv_io{ip};
    T=QPPv{ip}(io,:)'; % >> transposed
    no=length(io); clear io; 
   
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
        if ~rem(i,10000), fprintf('clst-QPP%d-vx%d\n',ip,i); end
        c0=Ti(:,:,i)*Tj(:,i+1:end)/PL;
        [c,imx]=max(abs(c0));
        for j=1:length(imx), c(j)=c(j)*sign(c0(imx(j),j)); end
        cij(1+(i-1)*(no-i/2):i*(no-(i+1)/2))=c;
    end; fprintf('clst-QPP%d %dsec\n\n',ip,round(toc));
    clear c0 c Ti Tj;
    QPPv_Lnkg{ip}=linkage(1-cij','average'); clear cij;  
end
save(p2V,'QPPv_Lnkg','-append');

%% Clustering
fprintf('clst-QPPs-rest\n');
QPPv_clst=nan(nvx,nP,'single'); 
QPPv_clst_te=zeros(2,nclst,nP,'single');
for ip=1:nP
    
    %% clustering & getting indices for sorting based on size
    I=cluster(QPPv_Lnkg{ip},'Cutoff',clstth(ip),'Criterion','Distance');
    [~,ih]=sort(hist(I,1:max(I)),'descend');
    
    %% sorting clusters based on size & keeping first nclst 
    nclst1=min(nclst,max(I));
    a=nan(size(I)); for i=1:nclst1, a(I==ih(i))=i; end
    I=nan(nvx,1,'single'); I(QPPv_io{ip})=a; clear ih a

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
    clear T 

    [~,i2]=sort(t(2,:)); % getting current cluster indices of ascending td
    [~,i1]=sort(t(1,i2)); ia=i2(i1); % ~ with 1st ascnd. tp then ascnd. td
    clear i1 i2
    
    %% combining current indices such that 1)DMN 2)~DMN 3-end)ascend tp&td
    ib=zeros(1,2); for i=1:2, ib(i)=find(ia==iM(i)); end
    ib=sort(ib); ibr=ib(1)+1:ib(2)-1; nb=length(ibr);
    iM(3:3+nb-1)=ia(ibr);
    iM(3+nb:nclst1)=[ia(ib(2)+1:end) ia(1:ib(1)-1)]; clear ia ib ibr nb

    %% reordering
    M=I; for i=1:nclst1, M(I==iM(i))=i; end; QPPv_clst(:,ip)=M; 
    
    %%
    ic=1:nclst1; QPPv_clst_te(:,ic,ip)=t(:,iM); clear T t I M iM ic
    
end
save(p2V,'QPPv_clst','QPPv_clst_te','-append');

%% Timecourses per cluster per region & their times of peak & 
%% number of timecourses per cluster per existing parcels of each region
%%
c=cell(nP,1); c(:)={cell(nrgn,nclst)}; 
QPPv2clst=c; QPPv2clst_tp=c; 
nvx2clst=nan(nrgn,nclst,nP,'single');
nvx2clst2prcl=c; clear c
for ip=1:nP
for ig=1:nrgn
    I=nan(nvx,1,'single'); I(irgn{ig})=QPPv_clst(irgn{ig},ip);
    R=PR{ig}(ivx); nR=nPR(ig); 
    for ic=1:nclst
        io=find(I==ic); nvx2clst2prcl{ip}{ig,ic}=nan(nR,1,'single');
        if any(io)
            QPPv2clst{ip}{ig,ic}=QPPv{ip}(io,:);
            QPPv2clst_tp{ip}{ig,ic}=QPPv_tp(io,ip);
            nvx2clst(ig,ic,ip)=length(io);
            nvx2clst2prcl{ip}{ig,ic}=hist(R(io),1:nR);
        end
    end
end
end
save(p2V,'QPPv2clst','QPPv2clst_tp','nvx2clst','nvx2clst2prcl','-append');
