
%% Testing significance of time of peak differences between: (1) seven
%% brain regions, (2) QPP clusters per region, (3) cortical RSNs within QPP
%%
clear; clc; p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'p2V','nP','nSD','nsd','PL','nvx','nrgn','irgn','nclst',...
    'p2u','nY','ivx','PLe','PLc','ntst'); isd=1; nth=15; nrc=nrgn*nclst; 
benf1=nrgn*(nrgn-1)/2; alph1=0.01/benf1; 
benf2=nrc*(nrc-1)/2; alph2=0.01/benf2; 
benf3=nY*(nY-1)/2; alph3=0.01/benf3;
load([p2u 'myGlssr.mat'],'Y','YLB'); Y=Y(ivx);
load(p2V,'QPPv_tp','QPPv_htp','QPPv2clst_tp','QPPv','QPPv_io'); 

%% (1) Timing differences between seven brain regions
%%
fprintf('PreTest1\n');
QPPv_tm=zeros(nP,1,'single'); ntp2=nan(nrgn,nP,'single');
QPPv_tp2=cell(nrgn,nP); QPPv_tp2(:)={single([])};
mstp2=nan(nrgn,2,nP,'single');
for ip=1:nP
    [~,tm]=findpeaks(-smooth(QPPv_htp(:,1,ip)));
    [~,im]=min(abs(tm-PL/2)); tm=tm(im); QPPv_tm(ip)=tm;
    tp=QPPv_tp(:,ip);
    for ig=1:nrgn
        tp2=nan(nvx,1,'single'); tp2(irgn{ig})=tp(irgn{ig});
        tp2(isnan(tp2))=[]; tp2(tp2<tm)=[]; n=length(tp2); 
        if n>=nth
            ntp2(ig,ip)=n; 
            QPPv_tp2{ig,ip}=tp2; 
            mstp2(ig,:,ip)=[median(tp2) std(tp2)/sqrt(n)];
        end
    end
end; clear tm tp tp2 n
save(p2V,'QPPv_tm','ntp2','mstp2','-append');

tic; fprintf('Test1 ');
pttst_tp2=nan(nrgn,nrgn,nP,ntst,'single'); 
for it=1:ntst
for ip=1:nP
for ig1=1:nrgn
    tp1=QPPv_tp2{ig1,ip}; l1=length(tp1);
    if l1
    for ig2=1:nrgn
        tp2=QPPv_tp2{ig2,ip}; l2=length(tp2);
        if l2
            l=min(l1,l2);
            if l==l1, t1=tp1; t2=tp2(randperm(l2)); t2=t2(1:l);
            else, t2=tp2; t1=tp1(randperm(l1)); t1=t1(1:l); end
            [h,p]=ttest(t1,t2,'Alpha',alph1);
            if ~h, pttst_tp2(ig1,ig2,ip,it)=1; 
            else, pttst_tp2(ig1,ig2,ip,it)=p; end
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
ntp=nan(nrgn,nclst,nP,'single');
mstp=nan(nrgn,nclst,2,nP,'single'); 
for ip=1:nP
for ig=1:nrgn
for ic=1:nclst
    tp=QPPv2clst_tp{ip}{ig,ic};
    if ~isempty(tp)
        t=tp; t(isnan(t))=[];
        ntp(ig,ic,ip)=length(t);
        mstp(ig,ic,:,ip)=[median(t) std(t)/sqrt(length(t))];
    end
end
end
end; clear tp t
ntp(ntp<nth)=nan;
save(p2V,'ntp','mstp','-append');

icR=cell(nrgn,nP);
for ip=1:nP
for ig=1:nrgn
    icR{ig,ip}=find(~isnan(ntp(ig,:,ip)));
end
end

tic; fprintf('Test2 ');
pttst_tp=nan(nrc,nrc,nP,ntst,'single'); 
for it=1:ntst
for ip=1:nP
for ig1=1:nrgn
    icr1=icR{ig1,ip}; 
    if ~isempty(icr1)
    for ic1=icr1
        tp1=QPPv2clst_tp{ip}{ig1,ic1}; 
        tp1(isnan(tp1))=[]; l1=length(tp1);
        for ig2=1:nrgn
            icr2=icR{ig2,ip};
            if ~isempty(icr2)
            for ic2=icr2
                tp2=QPPv2clst_tp{ip}{ig2,ic2};
                tp2(isnan(tp2))=[]; l2=length(tp2);
                l=min(l1,l2);
                if l==l1, t1=tp1; t2=tp2(randperm(l2)); t2=t2(1:l);
                else, t2=tp2; t1=tp1(randperm(l1)); t1=t1(1:l); 
                end 
                [h,p]=ttest(t1,t2,'Alpha',alph2);
                i=(ig2-1)*nclst+ic2; j=(ig1-1)*nclst+ic1;
                if ~h, pttst_tp(i,j,ip,it)=1;
                else, pttst_tp(i,j,ip,it)=p;
                end
            end
            end
        end
    end
    end
end
end
end; fprintf('%dsec\n',round(toc)); clear tp1 tp2 t1 t2 l1 l2 l h p
save(p2V,'pttst_tp','-append');

%% (3) Timing differences between cortical RSNs within QPP
%%
fprintf('PreTest3\n');
QPPv2y=nan(nY,PLe,2,nP,'single'); 
QPPv2y_tp=cell(nP,1); QPPv2y_tp(:)={cell(nY,1)};
msytp=nan(nY,2,nP,'single'); nytp=nan(nY,nP,'single'); 
for ip=1:nP
    
    T=QPPv{ip}; io=QPPv_io{ip};
    ix=1:nvx; ix(io)=[]; T(ix,:)=nan;
    
    for iy=1:nY
        T1=T(Y==iy,:); T1(isnan(sum(T1,2)),:)=[]; 
        l1=size(T1,1); nytp(iy,ip)=l1;
        QPPv2y(iy,:,1,ip)=mean(T1);
        QPPv2y(iy,:,2,ip)=std(T1)/sqrt(l1);

        tp=QPPv_tp(Y==iy,ip); 
        QPPv2y_tp{ip}{iy}=tp;
        msytp(iy,:,ip)=[median(tp,'omitnan') std(tp,'omitnan')/sqrt(l1)];
    end
end; clear T io ix tp T1 l1
save(p2V,'QPPv2y','msytp','nytp','-append');

tic; fprintf('Test3 ');
pttsty_tp=nan(nY,nY,nP,ntst,'single'); 
for it=1:ntst
for ip=1:nP
for iy1=1:nY
    tp1=QPPv2y_tp{ip}{iy1};
    tp1(isnan(tp1))=[]; l1=length(tp1);
    for iy2=1:nY
        tp2=QPPv2y_tp{ip}{iy2};
        tp2(isnan(tp2))=[]; l2=length(tp2);
        l=min(l1,l2);
        if l==l1, t1=tp1; t2=tp2(randperm(l2)); t2=t2(1:l);
        else, t2=tp2; t1=tp1(randperm(l1)); t1=t1(1:l);
        end
        [h,p]=ttest(t1,t2,'Alpha',alph3);
        if ~h, pttsty_tp(iy1,iy2,ip,it)=1;
        else, pttsty_tp(iy1,iy2,ip,it)=p;
        end
    end
end
end
end; fprintf('%dsec\n',round(toc)); clear tp1 tp2 t1 t2 l l1 l2 h p
save(p2V,'pttsty_tp','-append');
