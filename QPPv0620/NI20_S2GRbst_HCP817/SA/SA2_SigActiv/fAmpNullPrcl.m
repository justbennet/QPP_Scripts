
function PN=fAmpNullPrcl(TN,nmxn,p2B)
p2p=dir('../Params_*.mat'); p2p=p2p.name;
load(['../' p2p],'nsbj','nscn','nt','nT','nX','p2qppf','PLe','PL','PLh'); 
addpath(['../' p2qppf]); [nrep,n]=size(TN); 
fprintf('Loading Data\n'); load(p2B,'B');

PN=cell(nrep,n); PN(:)={zeros(nX,PLe,'single')};
for is=1:nsbj
    fprintf('Sbj%d\n',is);
    D=zeros(nX,nT,'single'); 
    for iscn=1:nscn, D(:,(iscn-1)*nt+(1:nt))=B{is,iscn}; end
    for ir=1:nrep
    for in=1:n
        if ~isempty(TN{ir,in}{is})
            T=Tbld(D,TN{ir,in}{is},PL,PLh,0);
            PN{ir,in}=PN{ir,in}+T;
        end
    end
    end
end
for ir=1:nrep, for in=1:n, PN{ir,in}=PN{ir,in}/nmxn(in); end; end

