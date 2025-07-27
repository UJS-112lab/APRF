function [RTVx,RTVy,RTVl,RTVr] = MyRTV(img,sigma)
%     img=im2single(img);
    
    dx = [img(:,1,:) - img(:,end,:), diff(img,1,2)];
    dy = [img(1,:,:) - img(end,:,:); diff(img,1,1)];
    
    Sx1=[img(:,end,:), img(:,1:end-1,:)];%右移1
    S1=[Sx1(end,:,:); Sx1(1:end-1,:,:)];%下移1，左上角
    
    Sx2=[img(:,2:end,:), img(:,1,:)];%左移1
    S2=[Sx2(end,:,:); Sx2(1:end-1,:,:)];%下移1，右上角

    d1=img-S1;%与左上角像素差
    d2=img-S2;%与右上角像素差
    
    Px(:,:,1)=imgaussfilt(dx(:,:,1),sigma);%高斯滤波
    Px(:,:,2)=imgaussfilt(dx(:,:,2),sigma);%高斯滤波
    Px(:,:,3)=imgaussfilt(dx(:,:,3),sigma);%高斯滤波
    
    Py(:,:,1)=imgaussfilt(dy(:,:,1),sigma);%高斯滤波
    Py(:,:,2)=imgaussfilt(dy(:,:,2),sigma);%高斯滤波
    Py(:,:,3)=imgaussfilt(dy(:,:,3),sigma);%高斯滤波
    
    P1(:,:,1)=imgaussfilt(d1(:,:,1),sigma);%高斯滤波
    P1(:,:,2)=imgaussfilt(d1(:,:,2),sigma);%高斯滤波
    P1(:,:,3)=imgaussfilt(d1(:,:,3),sigma);%高斯滤波
    
    P2(:,:,1)=imgaussfilt(d2(:,:,1),sigma);%高斯滤波
    P2(:,:,2)=imgaussfilt(d2(:,:,2),sigma);%高斯滤波
    P2(:,:,3)=imgaussfilt(d2(:,:,3),sigma);%高斯滤波
    
    Lx=abs(Px);
    Ly=abs(Py);
    
    L1=abs(P1);
    L2=abs(P2);
    
    absDx=abs(dx);
    absDy=abs(dy);
    
    absD1=abs(d1);
    absD2=abs(d2);
    
    Dx(:,:,1)=imgaussfilt(absDx(:,:,1),sigma);%高斯滤波
    Dx(:,:,2)=imgaussfilt(absDx(:,:,2),sigma);%高斯滤波
    Dx(:,:,3)=imgaussfilt(absDx(:,:,3),sigma);%高斯滤波
    
    Dy(:,:,1)=imgaussfilt(absDy(:,:,1),sigma);%高斯滤波
    Dy(:,:,2)=imgaussfilt(absDy(:,:,2),sigma);%高斯滤波
    Dy(:,:,3)=imgaussfilt(absDy(:,:,3),sigma);%高斯滤波
    
    D1(:,:,1)=imgaussfilt(absD1(:,:,1),sigma);%高斯滤波
    D1(:,:,2)=imgaussfilt(absD1(:,:,2),sigma);%高斯滤波
    D1(:,:,3)=imgaussfilt(absD1(:,:,3),sigma);%高斯滤波
    
    D2(:,:,1)=imgaussfilt(absD2(:,:,1),sigma);%高斯滤波
    D2(:,:,2)=imgaussfilt(absD2(:,:,2),sigma);%高斯滤波
    D2(:,:,3)=imgaussfilt(absD2(:,:,3),sigma);%高斯滤波
    
    ee=0.0001;
    
    Dx = mean(Dx, 3);
    Dy = mean(Dy, 3);
    D1 = mean(D1, 3);
    D2 = mean(D2, 3);
    Lx = mean(Lx, 3);
    Ly = mean(Ly, 3);
    L1 = mean(L1, 3);
    L2 = mean(L2, 3);
    
    RTVx=(Dx+ee)./(Lx+ee);
    RTVy=(Dy+ee)./(Ly+ee);
    
    RTVl=(D1+ee)./(L1+ee);
    RTVr=(D2+ee)./(L2+ee);


end