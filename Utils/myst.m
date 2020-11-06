function [md,mn,mx]=myst(a,r); r=10^r;
md=round(median(a,'omitnan')*r)/r;
mn=floor(min(a)*r)/r;
mx=ceil(max(a)*r)/r;