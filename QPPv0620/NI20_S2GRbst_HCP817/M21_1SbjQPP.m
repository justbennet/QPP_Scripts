function M21_1SbjQPP(IS)

%% Subject-level QPPs
% Per subject, the following steps are performed here:
% (1) Detecting the dominant QPP (QPPf1detect)
% (2) Phase-adjusting this QPP (QPPf2phadj) based on 1 or 2 seed parcels
% (3) Building QPP2 & above by averaging the contributing segments over the
% original scans NOT the residual scans (QPPf3xresid)
% (4) After detection of QPPi, regressing QPPs 1 to i scanwise(QPPf4regscn)
% to detect the next dominant QPP or QPPi+1 using the same steps
% > In the QPPf4regscn, FC after regression of QPPs 1-i is also calculated

%%
p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nsbj','nscn','nX','nt','nT','tres','p2B','p2S','p2qppf',...
    'nP','PL','cth','ncth1','nitr','ssg','ITP','ITPstp',...
    'PLh','tsh','PLc','cthph','sdph','ibY','iG2Y'); 
load(p2B,'B'); addpath(p2qppf); if nargin<1, IS=1:nsbj; end

%%
for is=IS
    
    D1=zeros(nX,nT,'single');
    for iscn=1:nscn, D1(:,(iscn-1)*nt+(1:nt))=B{is,iscn}; end; D=D1;
    TT=cell(nP,1); CT=zeros(nP,nT,'single');
    
for ip=1:nP  
    
    tic; if ip~=1, D=load(p2S{is,ip-1},'Dr'); D=D.Dr; end
    isip=sprintf('Sbj%d-QPP%d-',is,ip); 
    
    [QPP,TMX,C,MET,ITER,TMPL,TMXTMPL,CTMPL,SCMX]=QPPf1detect...
        (D,nscn,PL,cth{ip},ncth1(ip),nitr,ssg(ip),ITP{ip},PLh,tres,...
        [isip 'f1detect'],ITPstp);
    
    fprintf([isip 'f2phadj\n']);
    [QPPa,TMXa,Ca,METa,SDa,SD,flga,cT1Tj,nsim]=QPPf2phadj...
        (QPP,TMPL,SCMX,TMXTMPL,CTMPL,tsh,PLc,cthph,sdph{ip},1,tres);

    fprintf([isip 'f3xresid\n']);
    if ip~=1
        QPP0=QPP; QPPa0=QPPa; SD0=SD; SD=nan; SDa0=SDa; SDa=nan;
        [QPP,C,MET(1),QPPa,Ca,METa(:,1)]=QPPf3xresid...
            (D1,TMX,TMXa,PLh,PLc,nscn,ssg(ip));
    end

    fprintf([isip 'f4regscn\n']);
    TT{ip}=QPP; CT(ip,:)=C;
    [Dr,Cr,FCr]=QPPf4regscn...
        (D1,TT(1:ip),CT(1:ip,:),nscn,PL,PLc,ibY,iG2Y,0);   

    save(p2S{is,ip},'QPP','TMX','C','MET','ITER','TMPL','TMXTMPL',...
        'CTMPL','SCMX','QPPa','TMXa','Ca','METa','SDa','SD','flga',...
        'cT1Tj','nsim','Dr','Cr','FCr'); 
    if ip~=1, save(p2S{is,ip},'QPP0','QPPa0','SD0','SDa0','-append'); end
    fprintf('Sbj%d-QPP%d %dsec\n\n',is,ip,round(toc));
        
end
end; clear D1 D
