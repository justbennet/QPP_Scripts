function [TMXN,qN,rssN,QPPN]=QPPf6null(D,nmx,nscn,nnll,PL)
[nX,nT]=size(D); nt=nT/nscn; nP=length(nmx);

a=repmat((1:nt-PL+1)',1,nscn);
for i=1:nscn, a(:,i)=(i-1)*nt+a(:,i); end; a=single(a(:));

TMXN=cell(nP,nnll); QPPN=zeros(nX,PL,nP,nnll,'single');
rssN=zeros(nX,nP,nnll,'single'); qN=zeros(nP,2,'single');
for ip=1:nP
    for in=1:nnll
        b=a(randperm(length(a)));
        d=abs(diff(b)); I=find(d<=PL); ix=I;
        while any(I), d(I)=nan; I=find(d<=PL); ix=[ix;I]; end
        b(ix)=[];
        tmx=b(1:nmx(ip)); 
        tmx=sort(tmx); 
        TMXN{ip,in}=tmx;
        T=Tbld(D,tmx,PL,[0 0],1); QPPN(:,:,ip,in)=T;
        rssN(:,ip,in)=sqrt(sum(T.^2,2));
    end
    z=rssN(:,ip,:); qN(ip,1)=quantile(z(:),0.99);
    z=abs(QPPN(:,:,ip,:)); qN(ip,2)=quantile(z(:),0.99);
end
