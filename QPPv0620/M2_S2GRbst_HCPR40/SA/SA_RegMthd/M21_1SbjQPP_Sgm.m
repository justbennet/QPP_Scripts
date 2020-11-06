
function M21_1SbjQPP_Sgm(IS)
%%
p2p=dir('Params_*.mat'); p2p=p2p.name;
load(p2p,'nscn','nX','nt','nT','tres','p2B','p2S','p2qppf',...
    'nP','PL','cth','ncth1','nitr','ssg','ITP','ITPstp',...
    'PLh','tsh','PLc','cthph','sdph','ibY','iG2Y');
addpath(p2qppf);

%%
load(p2B,'B'); 
for is=IS
    D1=zeros(nX,nT,'single');
    for iscn=1:nscn, D1(:,(iscn-1)*nt+(1:nt))=B{is,iscn}; end; D=D1;
    TT=cell(nP,1); 
    
for ip=1:nP  
    
    tic; if ip~=1, D=load(p2S{is,ip-1},'Dr'); D=D.Dr; end
    isip=sprintf('Sbj%d-QPP%d-',is,ip); 
    
    [QPP,TMX,C,MET,~,TMPL,TMXTMPL,CTMPL,SCMX]=QPPf1detect...
        (D,nscn,PL,cth{ip},ncth1(ip),nitr,ssg(ip),ITP{ip},PLh,tres,...
        [isip 'f1detect'],ITPstp);
    
    fprintf([isip 'f2phadj\n']);
    [QPPa,TMXa,Ca,METa]=QPPf2phadj...
        (QPP,TMPL,SCMX,TMXTMPL,CTMPL,tsh,PLc,cthph,sdph{ip},1,tres);

    fprintf([isip 'f3xresid\n']);
    if ip~=1
        [QPP,C,MET(1),QPPa,Ca,METa(:,1)]=QPPf3xresid...
            (D1,TMX,TMXa,PLh,PLc,nscn,ssg(ip));
    end

    fprintf([isip 'f4regsgm\n']);
    TT{ip}=QPP; 
    [Dr,Cr,FCr]=QPPf4regsgm...
        (D1,TT(1:ip),nscn,PL,PLc,cth{ip}(2),tsh,ibY,iG2Y,0);

    save(p2S{is,ip},'QPP','TMX','C','MET','QPPa','TMXa','Ca','METa',...
        'Dr','Cr','FCr'); 
    fprintf('Sbj%d-QPP%d %dsec\n\n',is,ip,round(toc));
        
end
end; clear D1 D
