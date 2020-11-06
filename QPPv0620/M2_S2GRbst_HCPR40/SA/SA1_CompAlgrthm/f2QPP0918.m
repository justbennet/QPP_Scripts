function [FTP,TMPL,iter]=f2QPP0918(B,nd,WL,ITP,cth,n_iter_th1,mx_iter)
[nRG,nT]=size(B); nt=nT/nd; nch=nt-WL+1; nTf=nRG*WL;
    
bchf=cell(1,nT); bchfn=cell(1,nT);
for id=1:nd
for ich=1:nch
    T=B(:,(id-1)*nt+ich:(id-1)*nt+WL+ich-1); T=T(:); bchf{(id-1)*nt+ich}=T;
    T=T-sum(T)/nTf; T=T/sqrt(T'*T); bchfn{(id-1)*nt+ich}=T;
end
end
 
c=zeros(nT,1,'single');
for id=1:nd
for ich=1:nch
	c((id-1)*nt+ich)=bchfn{ITP}'*bchfn{(id-1)*nt+ich};
end
end
[~,icmx]=findpeaks(double(c),'MinPeakHeight',cth(1));
for id=1:nd, icmx( icmx==(id-1)*nt+1 | icmx==(id-1)*nt+nch )=[]; end

c0=c; c00=c; c000=c; iter=1;
while iter<=mx_iter
    c=smooth(c);
    if iter<=n_iter_th1, ith=1; else, ith=2; end; th=cth(ith);
    
    tpsgth=icmx; n_tpsgth=length(tpsgth);
    if n_tpsgth<=1, break; end
    
    T=bchf{tpsgth(1)};
    for i=2:n_tpsgth
        T=T+bchf{tpsgth(i)}; end; T=T/n_tpsgth;
    T=T-sum(T)/nTf; T=T/sqrt(T'*T);
    
    for id=1:nd
    for ich=1:nch
    	c((id-1)*nt+ich)=T'*bchfn{(id-1)*nt+ich};
    end
    end
    [~,icmx]=findpeaks(double(c),'MinPeakHeight',th);
    for id=1:nd
        icmx( icmx==(id-1)*nt+1 | icmx==(id-1)*nt+nch )=[];
        icmx( icmx>=(id-1)*nt+nch+1 & icmx<=id*nt )=[];
    end
    
    if (ncc(c0,c)>.9999) || (ncc(c00,c)>.9999) || (ncc(c000,c)>.9999)
        break; end
    c000=c00; c00=c0; c0=c; iter=iter+1;
end
if n_tpsgth>1, FTP=tpsgth; T=zeros(nRG,WL,'single');
    for i=1:n_tpsgth, T=T+B(:,tpsgth(i):tpsgth(i)+WL-1);end; TMPL=T/n_tpsgth;
end
end
function z =ncc(X,Y)
    X =X(:)-mean(X(:));
    Y =Y(:)-mean(Y(:));
    if(norm(X)==0 || norm(Y)==0), z = nan; return; end
    z = (X'*Y)/norm(X)/norm(Y);
end
