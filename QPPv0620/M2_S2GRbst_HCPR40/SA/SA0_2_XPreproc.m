
%% Preprocessing: averaging over non-preprocessed scans
%%
clear; clc; p1='../'; 
p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load([p1 p2p],'p2prep','p2O','p2qppf','p2u','d2SA','d2SAcft',...
    'nsbj','nscn','nt','nT','nvx','ivx','nVX','nP',...
    'nSD','nsd','PL','PLh','PLe'); isd=1;
p2s=[d2SA 'SA02.mat'];
load([p1 p2prep],'p2b'); load([p1 p2O],'TMX');
addpath([p1 p2qppf]); addpath(genpath([p1 p2u]));

%%
QPPvxp=cell(nP,1); QPPvxp(:)={zeros(nvx,PLe,'single')}; 
nmx=zeros(nP,nsbj,'single');
for is=1:nsbj
    fprintf('Sbj%d ',is); tic;
    D=zeros(nvx,nT,'single'); 
    for iscn=1:nscn
        b=ft_read_cifti(p2b{is,iscn});
        b=b.dtseries(ivx,:);
        b=b-repmat(mean(b,2),1,nt);
        D(:,(iscn-1)*nt+(1:nt))=b;
    end; clear b   
    for ip=1:nP
        tmx=TMX{ip,isd}{is}; n=length(tmx); nmx(ip,is)=n;
        if n, T=Tbld(D,tmx,PL,PLh,0); QPPvxp{ip}=QPPvxp{ip}+T; end
    end; fprintf('%dsec\n',round(toc));
end; clear T D
for i=1:nP, QPPvxp{i}=QPPvxp{i}/sum(nmx(i,:)); end
for i=1:nP
    P=QPPvxp{i}; mx=max(abs(P),[],2); P(mx>1000,:)=nan; QPPvxp{i}=P;
end; clear P; save(p2s,'QPPvxp')

%%
e=ft_read_cifti('empty.dtseries.nii');
e.time=1:PLe; e.hdr.dim(6)=PLe; 
for i=1:nP
    T=nan(nVX,PLe,'single'); T(ivx,:)=QPPvxp{i}; e.dtseries=T;
    ft_write_cifti([d2SAcft 'QPPxp' num2str(i)],e,'parameter','dtseries'); 
end
