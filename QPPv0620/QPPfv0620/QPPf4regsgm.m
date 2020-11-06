function [Dr,Cr,FCr]=QPPf4regsgm(D,TT,nscn,PL,PLc,cth,tsh,ibY,iG2Y,nz)
[nX,nT]=size(D); nt=nT/nscn; nXL=nX*PL; npr=length(TT); nsh=2*tsh+1; nitr=5;

%%
Tsh=cell(npr,1); Tsh(:)={zeros(nXL,nsh,'single')}; 
for ir=1:npr    
    cnt=1;
    for s=-tsh:tsh
        if s>0, T=[zeros(nX,s,'single') TT{ir}(:,1:end-s)];
        elseif s<0, T=[TT{ir}(:,-s+1:end) zeros(nX,-s,'single')];
        else, T=TT{ir}; 
        end
        T=T(:,PLc); T=T(:); T=T-sum(T)/nXL; 
        Tsh{ir}(:,cnt)=T/sqrt(T'*T); cnt=cnt+1;
    end
end; clear T cnt s

%%
Dr0=cell(npr+1,nscn); for i=1:nscn, Dr0{1,i}=D(:,(i-1)*nt+1:i*nt); end
Dr=zeros(nX,nT,'single'); Cr=zeros(npr,nT,'single'); FCr=zeros(nX);

for iscn=1:nscn
%% Regressing QPPs 1 to npr
for ir=1:npr
    itr=1; Dr0{ir+1,iscn}=Dr0{ir,iscn};
    while itr<=nitr  
        Csh=zeros(nsh,nt,'single');
        for i=1:nsh
        for isg=1:nt-PL+1
            S=Dr0{ir+1,iscn}(:,isg:isg+PL-1); S=S(:); S=S-sum(S)/nXL; 
            S=S'/sqrt(S'*S); Csh(i,isg)=S*Tsh{ir}(:,i);
        end
        end
        
        warning off; 
        [p,ip]=findpeaks(double(abs(Csh(:))),'MinPeakHeight',cth); warning on; 
        if ~any(p), break; end
        [~,is]=sort(p,'descend'); ips=ip(is); [ISH,TP]=ind2sub([nsh,nt],ips);
        clear A; cnt=1;
        while any(TP)  
            A(cnt,1)=TP(1); A(cnt,2)=ISH(1);
            [~,ix,~]=intersect(TP,TP(1)-PL:TP(1)+PL-1);
            TP(ix)=[]; ISH(ix)=[]; cnt=cnt+1;
        end; A=sortrows(A);

        for i=1:size(A,1)
            tp=A(i,1); ish=A(i,2);
            y=Dr0{ir+1,iscn}(:,tp:tp+PL-1); y=y(:); y=detrend(y);
            x=Tsh{ir}(:,ish); x=detrend(x);
            beta=(x'*x)\x'*y; yr=y-x*beta; 
            Dr0{ir+1,iscn}(:,tp:tp+PL-1)=reshape(yr,nX,PL);
        end    
        itr=itr+1;
    end
end
Dr(:,(iscn-1)*nt+1:iscn*nt)=Dr0{npr+1,iscn};

%% Correlating QPPs 1-npr with residual scan (goodness of regression)
d=Dr0{npr+1,iscn}; 
for ir=1:npr
    T=TT{ir}(:,PLc); T=T(:); T=T-sum(T)/nXL; T=T'/sqrt(T'*T);
    for isg=1:nt-PL+1
        S=d(:,isg:isg+PL-1); S=S(:); S=S-sum(S)/nXL; S=S/sqrt(S'*S);
        Cr(ir,(iscn-1)*nt+isg)=T*S;
    end
end

%% FC of residuals (FCr), averaged across scans using Fisher transform
if nargin>6
    dy=zeros(nX,nt);
    for i=1:length(iG2Y), dy(ibY(i)+1:ibY(i+1),:)=d(iG2Y{i},:); end
    c=corr(dy'); c(c==1)=0.9999; FCr=FCr+0.5*log((1+c)./(1-c));
end
end
if nargin>6, FCr=FCr/nscn;
    if nz, FCr=(exp(2*FCr)-1)./(exp(2*FCr)+1); FCr(FCr>=0.9999)=1; end
end

