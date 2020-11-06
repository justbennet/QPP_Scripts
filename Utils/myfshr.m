function [mc,md,mn,mx]=myfshr(c,n)
cz=0.5*log((1+c)./(1-c)); cz=mean(cz,n,'omitnan');
mc=(exp(2*cz)-1)./(exp(2*cz)+1);
md=round(median(mc(:),'omitnan')*100)/100;
mn=floor(min(mc(:))*100)/100;
mx=ceil(max(mc(:))*100)/100;
