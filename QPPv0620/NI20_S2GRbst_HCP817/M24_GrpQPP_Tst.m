
%% Testing significance of time of peak/dip differences between: (1) seven 
%% brain regions, (2) QPP clusters per region, (3) cortical RSNs within QPP
%%
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'p2V','nP','nSD','nsd','PL','nvx','nrgn','irgn','nclst',...
    'p2u','nY','ivx','PLe','PLc','ntst'); nth=15; nrc=nrgn*nclst; 
benf1=nrgn*(nrgn-1)/2; alph1=0.01/benf1; 
benf2=nrc*(nrc-1)/2; alph2=0.01/benf2; 
benf3=nY*(nY-1)/2; alph3=0.01/benf3;
load([p2u 'myGlssr.mat'],'Y','YLB'); Y=Y(ivx);
load(p2V,'QPPv_te','QPPv_hte','QPPv2clst_te','QPPv','QPPv_io'); 

%% (1) Timing differences between seven brain regions
%%
fprintf('PreTest1\n');
QPPv_tm=zeros(nP,nSD,'single'); ntp2=nan(nrgn,nP,nSD,'single');
QPPv_tp2=cell(nrgn,nP,nSD); QPPv_tp2(:)={single([])};
mstp2=nan(nrgn,2,nP,nSD,'single');
for ip=1:nP
for isd=1:nsd(ip)
    [~,tm]=findpeaks(-smooth(QPPv_hte(:,1,1,ip,isd)));
    [~,im]=min(abs(tm-PL/2)); tm=tm(im); QPPv_tm(ip,isd)=tm;
    tp=QPPv_te(:,1,ip,isd);
    for ig=1:nrgn
        tp2=nan(nvx,1,'single'); tp2(irgn{ig})=tp(irgn{ig});
        tp2(isnan(tp2))=[]; tp2(tp2<tm)=[]; n=length(tp2); 
        if n>=nth
            ntp2(ig,ip,isd)=n; 
            QPPv_tp2{ig,ip,isd}=tp2; 
            mstp2(ig,:,ip,isd)=[median(tp2) std(tp2)/sqrt(n)];
        end
    end
end
end; clear tm tp tp2 n
save(p2V,'QPPv_tm','ntp2','mstp2','-append');

tic; fprintf('Test1 ');
pttst_tp2=nan(nrgn,nrgn,nP,nSD,ntst,'single'); 
for it=1:ntst
for ip=1:nP
for isd=1:nSD 
for ig1=1:nrgn
    tp1=QPPv_tp2{ig1,ip,isd}; l1=length(tp1);
    if l1
    for ig2=1:nrgn
        tp2=QPPv_tp2{ig2,ip,isd}; l2=length(tp2);
        if l2
            l=min(l1,l2);
            if l==l1, t1=tp1; t2=tp2(randperm(l2)); t2=t2(1:l);
            else, t2=tp2; t1=tp1(randperm(l1)); t1=t1(1:l); end
            [h,p]=ttest(t1,t2,'Alpha',alph1);
            if ~h, pttst_tp2(ig1,ig2,ip,isd,it)=1; 
            else, pttst_tp2(ig1,ig2,ip,isd,it)=p; end
        end
    end
    end
end
end
end
end; fprintf('%dsec\n',round(toc)); clear tp1 tp2 t1 t2 l1 l2 l h p
save(p2V,'pttst_tp2','-append');

%% (2) Timing differences between QPP clusters per region
%%
fprintf('PreTest2\n');
nte=nan(nrgn,nclst,2,nP,nSD,'single');
mste=nan(nrgn,nclst,2,2,nP,nSD,'single'); 
for ip=1:nP
for isd=1:nSD
for ig=1:nrgn
for ic=1:nclst
    te=QPPv2clst_te{ip,isd}{ig,ic};
    if ~isempty(te)
        for ie=1:2
            t=te(:,ie); t(isnan(t))=[];
            nte(ig,ic,ie,ip,isd)=length(t);
            mste(ig,ic,ie,:,ip,isd)=[median(t) std(t)/sqrt(length(t))];
        end
    end
end
end
end
end; clear te t
nte(nte<nth)=nan;
save(p2V,'nte','mste','-append');

icR=cell(nrgn,2,nP,nSD);
for ip=1:nP
for isd=1:nsd(ip)
for ie=1:2
for ig=1:nrgn
    icR{ig,ie,ip,isd}=find(~isnan(nte(ig,:,ie,ip,isd)));
end
end
end
end

pttst_te=nan(nrc,nrc,2,nP,nSD,ntst,'single'); tic; 
for it=1:ntst
    fprintf('Test2 rep%d\n',it);
for ip=1:nP
for isd=1:nsd(ip)
for ie=1:2
for ig1=1:nrgn
    icr1=icR{ig1,ie,ip,isd}; 
    if ~isempty(icr1)
    for ic1=icr1
        te1=QPPv2clst_te{ip,isd}{ig1,ic1}(:,ie); 
        te1(isnan(te1))=[]; le1=length(te1);
        for ig2=1:nrgn
            icr2=icR{ig2,ie,ip,isd};
            if ~isempty(icr2)
            for ic2=icr2
                te2=QPPv2clst_te{ip,isd}{ig2,ic2}(:,ie);
                te2(isnan(te2))=[]; le2=length(te2);
                le=min(le1,le2);
                if le==le1, t1=te1; t2=te2(randperm(le2)); t2=t2(1:le);
                else, t2=te2; t1=te1(randperm(le1)); t1=t1(1:le); 
                end 
                [h,p]=ttest(t1,t2,'Alpha',alph2);
                i=(ig2-1)*nclst+ic2; j=(ig1-1)*nclst+ic1;
                if ~h, pttst_te(i,j,ie,ip,isd,it)=1;
                else, pttst_te(i,j,ie,ip,isd,it)=p;
                end
            end
            end
        end
    end
    end
end
end
end
end
end; fprintf('Test2 %dsec\n',round(toc)); clear te1 te2 t1 t2 le1 le2 le h p
save(p2V,'pttst_te','-append');

%% (3) Timing differences between cortical RSNs within QPP
%%
fprintf('PreTest3\n');
QPPv2y=nan(nY,PLe,2,nP,nSD,'single'); 
QPPv2y_te=cell(nP,nSD); QPPv2y_te(:)={cell(nY,1)};
msyte=nan(nY,3,2,nP,nSD,'single'); nyte=nan(nY,nP,nSD,'single'); 
for ip=1:nP
for isd=1:nsd(ip)
    
    T=QPPv{ip,isd}; io=QPPv_io{ip,isd};
    ix=1:nvx; ix(io)=[]; T(ix,:)=nan;
    
    tz=nan(nvx,1,'single'); % time of zero-crossing
    for iv=io'
        z=T(iv,PLc); t=find(z.*circshift(z,-1)<=0); t(t==PL)=[];
        if length(t)==1, tz(iv)=t; end
    end
    
    for iy=1:nY
        T1=T(Y==iy,:); T1(isnan(sum(T1,2)),:)=[]; 
        l1=size(T1,1); nyte(iy,ip,isd)=l1;
        QPPv2y(iy,:,1,ip,isd)=mean(T1);
        QPPv2y(iy,:,2,ip,isd)=std(T1)/sqrt(l1);

        tez=[QPPv_te(Y==iy,:,ip,isd) tz(Y==iy)]; 
        QPPv2y_te{ip,isd}{iy}=tez;
        msyte(iy,:,1,ip,isd)=median(tez,'omitnan');
        msyte(iy,:,2,ip,isd)=std(tez,'omitnan')/sqrt(l1);
    end
end
end; clear T io ix t tz z tez T1 l1
save(p2V,'QPPv2y','msyte','nyte','-append');

tic; fprintf('Test3 ');
pttsty_te=nan(nY,nY,3,nP,nSD,ntst,'single'); 
for it=1:ntst
for ip=1:nP
for isd=1:nsd(ip)
for ie=1:3
for iy1=1:nY
    te1=QPPv2y_te{ip,isd}{iy1}(:,ie);
    te1(isnan(te1))=[]; l1=length(te1);
    for iy2=1:nY
        te2=QPPv2y_te{ip,isd}{iy2}(:,ie);
        te2(isnan(te2))=[]; l2=length(te2);
        l=min(l1,l2);
        if l==l1, t1=te1; t2=te2(randperm(l2)); t2=t2(1:l);
        else, t2=te2; t1=te1(randperm(l1)); t1=t1(1:l);
        end
        [h,p]=ttest(t1,t2,'Alpha',alph3);
        if ~h, pttsty_te(iy1,iy2,ie,ip,isd,it)=1;
        else, pttsty_te(iy1,iy2,ie,ip,isd,it)=p;
        end
    end
end
end
end
end
end; fprintf('%dsec\n',round(toc)); clear te1 te2 t1 t2 l l1 l2 h p
save(p2V,'pttsty_te','-append');

