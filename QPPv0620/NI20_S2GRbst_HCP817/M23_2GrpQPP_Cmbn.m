
%% Combining mat files
clear; clc; prfx='f';

%%
p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nP','nSD','nsd','nvx','PLe','nsbj','p2O','p2V','p2VG2S'); 
load(p2O,'NMX');

%%
QPPv=cell(nP,nSD); QPPv(:)={single([])};
for ip=1:nP, for isd=1:nsd(ip), QPPv{ip,isd}=zeros(nvx,PLe,'single'); end; end

a=dir([prfx '*' p2V]); n=size(a,1); 
for i=1:n
    b=load(a(i).name); 
    for ip=1:nP, for isd=1:nsd(ip)
            QPPv{ip,isd}=QPPv{ip,isd}+b.QPPv{ip,isd}; end; end
end

for ip=1:nP, for isd=1:nsd(ip)
    QPPv{ip,isd}=QPPv{ip,isd}/sum(NMX(ip,isd,:)); end; end
save(p2V,'QPPv');

%%
s=strfind(p2VG2S{1,1},'/'); aa=p2VG2S{1,1}(1:s(end)); ind=cell(n,1);
for ip=1:nP
for isd=1:nsd(ip) 
    a=dir([aa prfx '*' p2VG2S{ip,isd}(s(end)+1:end)]);
    QPPvg2s=cell(nsbj,1);
    for i=1:n
        tic; fprintf('Loading G2SQPP%dSD%d_%i',ip,isd,i); 
        b=load([aa a(i).name]); 
        fprintf(' %ds\n',round(toc)); 
        IS=b.IS; b=b.QPPvg2s; 
        QPPvg2s(IS)=b(IS);
    end
    tic; fprintf('Saving G2SQPP%dSD%d',ip,isd); 
    save(p2VG2S{ip,isd},'QPPvg2s','-v7.3'); fprintf(' %ds\n',round(toc));
end
end
