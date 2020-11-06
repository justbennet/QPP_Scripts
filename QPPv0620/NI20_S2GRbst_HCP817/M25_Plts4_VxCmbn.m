
%% Combining screenshots of QPPs (captured by workbench)
%%
clear; clc; close all;
ps='./Plts/cfts2/mv/';
pe='.png'; fmv='zmv_';
v='_s_'; o=[2 2; 2 1; 1 2]; no=size(o,1); v0='_C_';
nP=3; PL=30; PLh=PL/2;
zy=40; zx=20;
t=PLh+(1:PL); 
t1=PLh+[3 11 13 15:19 21 28]; tf={t1;t1;t1};
fs=45; ss=get(0,'screensize'); 

%%
ip=1; it=16;
f=[ps num2str(ip) v num2str(it) pe]; 
I=imread(f); I=I(:,2:end,:); 
[y,x,~]=size(I); x3=x/3; y2=y/2; x32=round(x3/2);

f=[ps num2str(ip) v0 num2str(it) pe]; 
I=imread(f);
d=80; I=I(d+1:end-d,:,:); I=imresize(I,[y2 nan]);

x0=size(I,2); x02=round(x0/2); x04=round(x0/4);
ZY=zeros(zy,x0+x,3,'uint8');
ZX=zeros((zy+y2)*nP+zy*5,zx,3,'uint8');

cb1=imread([ps 'clr1.png']); cb2=imread([ps 'clr2.png']); % ps(1:end-1)
[by,bx,~]=size(cb1); bx2=round(bx/2);

%% For Movie
for it=t
    M=[];
for ip=1:nP 
    f=[ps num2str(ip) v num2str(it) pe]; 
    I=imread(f); I=I(:,2:end,:);
    J=zeros(y2,x,3,'uint8');
    for io=1:no
        i=o(io,1); j=o(io,2);
        J(:,(io-1)*x3+1:io*x3,:)=I((i-1)*y2+1:i*y2,(j-1)*x3+1:j*x3,:); 
    end  
    f=[ps num2str(ip) v0 num2str(it) pe]; 
    I=imread(f); I=I(d+1:end-d,:,:); I=imresize(I,[y2 nan]); 
    K=[ZY; I J];
    K=insertText(K,[x02,zy+10],['QPP' num2str(ip)],'fontsize',fs+15,...
        'TextColor','white','BoxOpacity',0,'AnchorPoint','Center');
    M=[M; K];
end
    M=[ZY;ZY;M;ZY;ZY;ZY]; M=[ZX M ZX];
    M(end-zy-by+1:end-zy,zx+x02+x04-bx2+1:zx+x02+x04+bx2,:)=cb1;
    M(end-zy-by+1:end-zy,zx+x0+x3+x32-bx2+1:zx+x0+x3+x32+bx2,:)=cb2;
    [m1,m2,~]=size(M);
    M=insertText(M,[zx+x0+x32,m1-zy-20],'QPP`s amplitude','fontsize',fs,...
        'TextColor','white','BoxOpacity',0,'AnchorPoint','Center');
    M=insertText(M,[zx+x0+x32,zy+20],['QPP`s timepoint '...
        num2str(it-PLh)],'fontsize',fs,'TextColor','white',...
        'BoxOpacity',0,'AnchorPoint','Center');
%     figure; set(gcf,'Position',ss); imagesc(M); axis image
    imwrite(M,[ps fmv num2str(it) pe],'Mode','lossless');
end

%% For figure
for ip=1:nP 
    ps2=[ps num2str(ip)];
    F=[];
for it=tf{ip}
    f=[ps2 v num2str(it) pe]; 
    I=imread(f); I=I(:,2:end,:);
    J=zeros(y2,x,3,'uint8');
    for io=1:no
        i=o(io,1); j=o(io,2);
        J(:,(io-1)*x3+1:io*x3,:)=I((i-1)*y2+1:i*y2,(j-1)*x3+1:j*x3,:); 
    end    
    f=[ps2 v0 num2str(it) pe]; 
    I=imread(f); I=I(d+1:end-d,:,:); I=imresize(I,[y2 nan]);
    F=[F; I J];
end
    imwrite(F,[ps2 pe],'Mode','lossless');
end
