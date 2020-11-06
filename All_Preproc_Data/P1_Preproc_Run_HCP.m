
function P1_Preproc_Run_HCP(nmD,prfxin,ISin)

%% Preprocessing
% Per scan the following steps are performed: demeaning & filtering,
% regressing WM&CSF & parcellating , regressing GM&WM&CSF & parcellating 
% > after regressing WM&CSF or GM&WM&CSF, each scan is saved seperately
% > after parcellating, all scans are saved in a mat file, seperately for
% WCR & GWCR
% > Preprocessing of subjects can be divided across cpus/servers and the
% cell arrays containing parcellated scans later combined by P2_Bcmbn
% IS: indices of sbjs to be preprocessed per cpu/server, [] for all sbjs
% prfx: temporary prefix to save the mat file of IS subjects
% nmD: dataset name
% > Example: P1_Preproc_Run_HCP('HCPR40','f1_',1:10) or ('HCPR3')

%%
p2prep=['Prep_Params_' nmD '.mat'];
load(p2prep,'nsbj','nscn','ivx','nvx','nt','nzFltr','Fltr','nX','ixG',...
    'p2b','p2rg','p2bwcr','p2bgwcr','p2BWCR','p2BGWCR','p2u','INFO');
addpath(genpath(p2u)); 
prfx=''; IS=1:nsbj; if nargin>1, prfx=prfxin; IS=ISin; end

BWCR=cell(nsbj,nscn); BGWCR=cell(nsbj,nscn);
for is=IS
for iscn=1:nscn, fprintf('Sbj%dScn%d ',is,iscn), tic;
    
    b=ft_read_cifti(p2b{is,iscn}); 
    b=b.dtseries(ivx,:);
    b=b-repmat(mean(b,2),1,nt); 
    b0=[zeros(nvx,nzFltr) b zeros(nvx,nzFltr)]; clear b 
    bft=zeros(nt,nvx,'single');
    for i=1:nvx
        y0=filtfilt(Fltr{1},Fltr{2},b0(i,:));
        y=y0(nzFltr+1:end-nzFltr);
        bft(:,i)=y-sum(y)/nt;
    end; clear b0 
    
    rg=[dlmread(p2rg{is,iscn}{1}) dlmread(p2rg{is,iscn}{2})];
    rg=rg-repmat(mean(rg),nt,1); 
    rg0=[zeros(nzFltr,2); rg; zeros(nzFltr,2)];   
    rg=zeros(nt,2,'single');
    for i=1:2
        y0=filtfilt(Fltr{1},Fltr{2},rg0(:,i));
        rg(:,i)=y0(nzFltr+1:end-nzFltr);     
    end 
    
    rg=zscore(rg); beta=(rg'*rg)\rg'*bft; bpr=(bft-rg*beta)'; 
    save(p2bwcr{is,iscn},'bpr','-v7.3'); 
    BWCR{is,iscn}=zeros(nX,nt,'single');
    for i=1:nX, BWCR{is,iscn}(i,:)=zscore(mean(bpr(ixG{i},:))); end
    
    rg=[zscore(mean(bft,2)) rg]; beta=(rg'*rg)\rg'*bft; bpr=(bft-rg*beta)'; 
    save(p2bgwcr{is,iscn},'bpr','-v7.3'); 
    BGWCR{is,iscn}=zeros(nX,nt,'single');
    for i=1:nX, BGWCR{is,iscn}(i,:)=zscore(mean(bpr(ixG{i},:))); end
    
    clear bpr bft; t=toc; t=round(t); fprintf('%d sec\n',t)
end
end
B=BWCR; save([prfx p2BWCR],'B','INFO','-v7.3');
B=BGWCR; save([prfx p2BGWCR],'B','INFO','-v7.3');

