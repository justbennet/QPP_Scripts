function [CNTtrns,CNTovlp,CNTnp]=QPPstrns(TMXa,nT,nP1)
[nsbj,~]=size(TMXa);

IS=zeros(nsbj,nP1,'single'); 
for is=1:nsbj, for ip=1:nP1, if any (TMXa{is,ip}), IS(is,ip)=1; end; end; end
IS1=find(sum(IS,2)==nP1); nsbj1=length(IS1);

A=zeros(nsbj1,nT,'single'); CNTnp=zeros(nsbj1,nP1,'single');
z=zeros(nP1,'single'); CNTovlp=z; CNTtrns=z; clear z
for is=1:nsbj1
    for ip=1:nP1
        tmx=TMXa{IS1(is),ip}; CNTnp(is,ip)=length(tmx);
        for i=1:CNTnp(is,ip)
            a=A(is,tmx(i));
            if ~a, A(is,tmx(i))=ip;
            else, CNTovlp(a,ip)=CNTovlp(a,ip)+1;
            end
        end
    end
    a=A(is,:); a(~a)=[];
    for i=2:length(a)
        CNTtrns(a(i-1),a(i))=CNTtrns(a(i-1),a(i))+1;
    end
end
