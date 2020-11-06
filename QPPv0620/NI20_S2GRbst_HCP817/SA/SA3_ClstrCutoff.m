
%% Increasing cut-off when clustering QPP's timecourses
%%
clear; clc; close all; p1='../'; clstth=[0.25 0.25 0.4]; isd=1;
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2V','p2qppf','p2u','nvx','nP','nSD','nsd','PLe','PL'...
    ,'PLc','irf','nclst','clstclr','nVX','ivx','d2SAplt','d2SAcft');  
load([p1 p2V],'QPPv','QPPv_io','QPPv_Lnkg');
addpath(genpath([p1 p2u])); addpath([p1 p2qppf]);

alm=repmat(0.75*[-1 1],nP,1);
ss=get(0,'Screensize'); ssw=ss; ssw(2)=ss(4)/4; ssw(4)=ss(4)/4;
set(0,'DefaultAxesTitleFontWeight','normal');
e=ft_read_cifti('empty.dtseries.nii');

%% 
QPPv_clst=nan(nvx,nP,'single'); h=cell(nP,1); 
for ip=1:nP 
    I=cluster(QPPv_Lnkg{ip,isd},'Cutoff',clstth(ip),'Criterion','Distance');
    [h{ip},ih]=sort(hist(I,1:max(I)),'descend');
      
    nclst1=min(nclst,max(I));
    a=nan(size(I)); for i=1:nclst1, a(I==ih(i))=i; end
    I=nan(nvx,1,'single'); I(QPPv_io{ip,isd})=a; clear ih a

    T=zeros(nclst1,PLe); 
    t=zeros(2,nclst1);
    for i=1:nclst1
        y=mean(QPPv{ip}(I==i,:)); T(i,:)=y; 
        [~,t(1,i)]=max(y(PLc)); [~,t(2,i)]=min(y(PLc)); 
    end
    iM(1)=mode(I(irf{1}));
    [~,iM(2)]=min(corr(T(iM(1),:)',T')); clear T 

    [~,i2]=sort(t(2,:));
    [~,i1]=sort(t(1,i2)); ia=i2(i1); clear i1 i2
     
    ib=zeros(1,2); for i=1:2, ib(i)=find(ia==iM(i)); end
    ib=sort(ib); ibr=ib(1)+1:ib(2)-1; nb=length(ibr);
    iM(3:3+nb-1)=ia(ibr);
    iM(3+nb:nclst1)=[ia(ib(2)+1:end) ia(1:ib(1)-1)]; clear ia ib ibr nb
 
    M=I; for i=1:nclst1, M(I==iM(i))=i; end
    QPPv_clst(:,ip)=M; 
    
end

%%
figure; set(gcf,'Position',ss);
for ip=1:nP
    subplot(2,nP,ip), stem(h{ip}); xlim([0 length(h{ip})+1]);
    title(['QPP' num2str(ip) 'SD' num2str(isd)]);
    if ip==1, ylabel('cluster size'); end
    subplot(2,nP,nP+ip), stem(log10(h{ip})); xlim([0 length(h{ip})+1]);
    if ip==1, ylabel('log_{10} (cluster size)'); end
end; saveas(gcf,[d2SAplt 'SA3_clst_size_HiTh.png' ]); close;

MM=nan(nVX,nP);
for ip=1:nP
    M=QPPv_clst(:,ip); MM(ivx,ip)=M;
    figure; set(gcf,'Position',ssw);
    for i=1:nclst
        T=QPPv{ip,isd}(M==i,:); 
        if any(T)
            iT=randperm(size(T,1)); iT=iT(1:min(50,length(iT))); T=T(iT,:);
            subplot(1,nclst,i), plot(1:PLe,T,'color',clstclr(i,:));
            PLTUSD2(PLc,alm(ip,:),~(i-1),~(i-1)); 
        end
    end; saveas(gcf,[d2SAplt 'SA3_clst_HiTh_QPP' num2str(ip) '_1.png']); close 
end
e.dtseries=MM; e.time=1:nP; e.hdr.dim(6)=nP; 
ft_write_cifti([d2SAcft 'clst_HiTh_QPPs' ],e,'parameter','dtseries');
