
function fAmpNull(IS,p2s,TN)
p2p=dir('../Params_*.mat'); p2p=p2p.name;
load(['../' p2p],'nscn','nvx','nt','nT','p2bpr','p2qppf','PLe','PL','PLh'); 
addpath(['../' p2qppf]); [nrep,n]=size(TN); 

PN=cell(nrep,n); PN(:)={zeros(nvx,PLe,'single')};
for is=IS
    fprintf('Sbj%d ',is); tic;
    D=zeros(nvx,nT,'single'); 
    for iscn=1:nscn
        load(['../' p2bpr{is,iscn}],'bpr'); 
        D(:,(iscn-1)*nt+(1:nt))=zscore(bpr,[],2);
    end; clear bpr
    for ir=1:nrep
    for in=1:n
        if ~isempty(TN{ir,in}{is})
            T=Tbld(D,TN{ir,in}{is},PL,PLh,0); 
            PN{ir,in}=PN{ir,in}+T;
        end
    end
    end; fprintf('%dsec\n',round(toc));
end
save(p2s,'PN','-v7.3');
