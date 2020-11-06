function [QPPs,TMXs,METs,QPPsa,TMXsa,METsa,FCsr]=QPPscmbn(p2S)
[nsbj,nP]=size(p2S); load(p2S{1,1},'QPP'); nX=size(QPP,1);
nsd=zeros(1,nP); for i=1:nP, load(p2S{1,i},'QPPa'); 
    nsd(i)=length(QPPa); end; nSD=max(nsd);

c=cell(nsbj,nP); QPPs=c; TMXs=c; METs=zeros(nsbj,3,nP,'single'); 
c=cell(nsbj,nP,nSD); QPPsa=c; TMXsa=c; METsa=zeros(nsbj,3,nP,nSD,'single');
FCsr=cell(nP,1);
for ip=1:nP, c=zeros(nX);
for is=1:nsbj
    load(p2S{is,ip},'QPP','TMX','MET','QPPa','TMXa','METa','FCr');
    QPPs{is,ip}=QPP; TMXs{is,ip}=TMX; METs(is,:,ip)=MET; c=c+FCr;
    for i=1:nsd(ip)
        QPPsa{is,ip,i}=QPPa{i}; TMXsa{is,ip,i}=TMXa{i}; 
        METsa(is,:,ip,i)=METa(i,:);
    end
end; c=c/nsbj; c=(exp(2*c)-1)./(exp(2*c)+1); c(c>=0.9999)=1; FCsr{ip}=c;
end
