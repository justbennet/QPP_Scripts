
function fXPreproc(IS,p2s)
p1='../'; p2p=dir([p1 'Params_*.mat']); p2p=[p1 p2p.name];
load(p2p,'p2prep','p2O','p2qppf','p2u','nscn','nt','nT','nvx',...
    'ivx','nP','nSD','nsd','PL','PLh','PLe'); 
load([p1 p2prep],'p2b');
load([p1 p2O],'TMX','NMX');
addpath([p1 p2qppf]); addpath(genpath([p1 p2u]));

%%
QPPvxp=cell(nP,nSD); QPPvxp(:)={single([])};
for ip=1:nP, for isd=1:nsd(ip)
        QPPvxp{ip,isd}=zeros(nvx,PLe,'single'); end; end

for is=IS
    fprintf('Sbj%d ',is); tic;
    D=zeros(nvx,nT,'single'); 
    for iscn=1:nscn
        b=ft_read_cifti(p2b{is,iscn});
        b=b.dtseries(ivx,:);
        b=b-repmat(mean(b,2),1,nt);
        D(:,(iscn-1)*nt+(1:nt))=b;
    end; clear b   
    for ip=1:nP
    for isd=1:nsd(ip)
    if NMX(ip,isd,is)
        T=Tbld(D,TMX{ip,isd}{is},PL,PLh,0);
        QPPvxp{ip,isd}=QPPvxp{ip,isd}+T;
    end
    end
    end; fprintf('%ds\n',round(toc));
end; clear T D
save(p2s,'QPPvxp');

