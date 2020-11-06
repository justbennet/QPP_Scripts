
%% Obtaining group QPP in grayordinate, finding active vertices/voxels 
%% within group QPP & finding their times of peak (fine summary map)
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

%% 
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nsbj','nscn','nvx','nt','nT','p2bpr','p2O','p2V','p2NV',...
    'p2qppf','nP','PLe','PL','PLh','PLc','nnll'); isd=1;
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

qvN=zeros(nP,2,'single'); rssN=zeros(nvx,nP,nnll,'single'); 
for ip=1:nP
    QPPNv(:,:,ip,:)=QPPNv(:,:,ip,:)/nmxN1(ip);
    for in=1:nnll
        rssN(:,ip,in)=sqrt(sum(QPPNv(:,:,ip,in).^2,2)); 
    end
    qvN(ip,1)=quantile(reshape(rssN(:,ip,:),[],1),0.99);
    qvN(ip,2)=quantile(reshape(abs(QPPNv(:,:,ip,:)),[],1),0.99);
end
save(p2V,'qvN','-append'); % save(p2NV,'QPPNv','-v7.3'); 

%% Finding active vertices/voxels & times of peak 
%%
QPPv_tp=nan(nvx,nP,'single'); QPPv_io=cell(nP,1); 
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
      
end
save(p2V,'QPPv_io','QPPv_tp','-append');
