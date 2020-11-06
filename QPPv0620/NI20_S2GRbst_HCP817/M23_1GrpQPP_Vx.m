
function M23_1GrpQPP_Vx(IS,prfxin)
%%
% Averaging the contributing segments of GrpQPP (identified based on
% parcel-space analysis) over grayordinate scans
 
%% 
p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nsbj','nscn','nvx','nt','nT','nSD','nsd',... 
    'p2bpr','p2O','p2V','p2VG2S','p2qppf','nP','PLe','PL','PLh'); 
load(p2O,'TMX','NMX'); addpath(p2qppf);
prfx=''; if nargin>1, prfx=prfxin; end
a=strfind(p2VG2S{1,1},'/'); 
for ip=1:nP, for i=1:nsd(ip), p2VG2S{ip,i}=...
	[p2VG2S{ip,i}(1:a(end)) prfx p2VG2S{ip,i}(a(end)+1:end)]; end; end

%% 
QPPv=cell(nP,nSD); QPPv(:)={single([])};
for ip=1:nP, for isd=1:nsd(ip), QPPv{ip,isd}=zeros(nvx,PLe,'single'); end; end
qppvg2s=cell(nP,nSD); qppvg2s(:)={single([])};
for ip=1:nP, for isd=1:nsd(ip), qppvg2s{ip,isd}=cell(nsbj,1); end; end
for is=IS
    fprintf('Sbj%d ',is); tic;
    D=zeros(nvx,nT,'single'); 
    for iscn=1:nscn
        load(p2bpr{is,iscn},'bpr'); 
        D(:,(iscn-1)*nt+(1:nt))=zscore(bpr,[],2);
    end; clear bpr
    for ip=1:nP
    for isd=1:nsd(ip)
        if NMX(ip,isd,is)
            T=Tbld(D,TMX{ip,isd}{is},PL,PLh,0); 
            QPPv{ip,isd}=QPPv{ip,isd}+T; 
            qppvg2s{ip,isd}{is}=T/NMX(ip,isd,is);
        end
    end
    end; fprintf('%dsec\n',round(toc));
end; clear D TMX T
save([prfx p2V],'QPPv');

fprintf('\nSaving GrpQPP-proj2-Sbj\n');
for i=1:nP, for j=1:nsd(i)
    QPPvg2s=qppvg2s{i,j}; save(p2VG2S{i,j},'QPPvg2s','IS','-v7.3'); end; end

