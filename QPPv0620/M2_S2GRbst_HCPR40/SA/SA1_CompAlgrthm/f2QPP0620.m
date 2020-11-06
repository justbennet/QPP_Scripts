function [TMXTMPL,TMPL,itr]=f2QPP0620(D,nscn,PL,ITP,cth,ncth1,nitr)
[nX,nT]=size(D); nt=nT/nscn; esg=nt-PL+1; nXL=nX*PL; ssg=1; e=0.9999;

SGf=zeros(nXL,nT,'single'); SGfn=zeros(nXL,nT,'single');
for iscn=1:nscn
for isg=ssg:esg
    t=(isg:isg+PL-1)+(iscn-1)*nt; 
    S=D(:,t); S=S(:); SGf(:,(iscn-1)*nt+isg)=S;
    S=S-sum(S)/nXL; S=S/sqrt(S'*S); SGfn(:,(iscn-1)*nt+isg)=S; 
end
end; clear S

c=SGfn(:,ITP)'*SGfn;
[~,tmx]=findpeaks(c,'MinPeakHeight',cth(1));
for iscn=1:nscn
    tmx( tmx==ssg+(iscn-1)*nt | tmx==esg+(iscn-1)*nt)=[]; end
nmx=length(tmx);
cn=c-sum(c)/nT; cn=cn'/sqrt(cn*cn');
    
cn1=cn; cn2=cn; cn3=cn; itr=1;
while itr<=nitr
    if nmx<=1, break; end
    
    T=sum(SGf(:,tmx),2)/nmx; T=T-sum(T)/nXL; T=T/sqrt(T'*T);
    c=T'*SGfn;
    if itr<=ncth1-1, ith=1; else, ith=2; end
    [~,tmx]=findpeaks(c,'MinPeakHeight',cth(ith));
    for iscn=1:nscn
        tmx( tmx==ssg+(iscn-1)*nt | tmx==esg+(iscn-1)*nt)=[]; end
    nmx=length(tmx);
    cn=c-sum(c)/nT; cn=cn/sqrt(cn*cn');
    
    if cn*cn1>e || cn*cn2>e || cn*cn3>e, break; end
    cn3=cn2; cn2=cn1; cn1=cn'; itr=itr+1;
end; clear cn cn1 cn2 cn3

TMPL=[]; TMXTMPL=[];
if nmx>1
    T=zeros(nX,PL,'single');
    for i=1:nmx, T=T+D(:,tmx(i):tmx(i)+PL-1); end; TMPL=T/nmx;
    TMXTMPL=single(tmx)';
end