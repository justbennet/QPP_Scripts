function [TMXN,qN,rssN,QPPN]=QPPsnull(B,nmx,nnll,PL)
[nsbj,nscn]=size(B); [nX,nt]=size(B{1,1}); nP=length(nmx);
nscng=nsbj*nscn; nT=nscn*nt;

a=repmat((1:nt-PL+1)',1,nscng);
for i=1:nscng, a(:,i)=(i-1)*nt+a(:,i); end; a=single(a(:));

TMXN=cell(nP,nnll); TMXN(:)={cell(nsbj,1)};
QPPN=zeros(nX,PL,nP,nnll,'single');
rssN=zeros(nX,nP,nnll,'single'); qN=zeros(nP,2,'single');
for ip=1:nP
    for in=1:nnll
        b=a(randperm(length(a)));
        d=abs(diff(b)); I=find(d<=PL); ix=I;
        while any(I), d(I)=nan; I=find(d<=PL); ix=[ix;I]; end
        b(ix)=[];
        tmx=b(1:nmx(ip)); 
        tmx=sort(tmx); 
        for is=1:nsbj
            t=tmx; t( t<(is-1)*nT+1 | t>is*nT )=[]; t=t-(is-1)*nT;
            TMXN{ip,in}{is}=t;
            D=zeros(nX,nT,'single'); 
            for i=1:nscn, D(:,(i-1)*nt+1:i*nt)=B{is,i}; end
            QPPN(:,:,ip,in)=QPPN(:,:,ip,in)+Tbld(D,t,PL,[0 0],0);
        end
        QPPN(:,:,ip,in)=QPPN(:,:,ip,in)/nmx(ip);
        rssN(:,ip,in)=sqrt(sum(QPPN(:,:,ip,in).^2,2));
    end
    z=rssN(:,ip,:); qN(ip,1)=quantile(z(:),0.99);
    z=abs(QPPN(:,:,ip,:)); qN(ip,2)=quantile(z(:),0.99);
end
