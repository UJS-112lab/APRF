function [h_input,v_input,l_input,r_input] = MyMapping(img,sigma1,alpha1)

% h_input = [diff(img,1,2), img(:,1,:) - img(:,end,:)];
% v_input = [diff(img,1,1); img(1,:,:) - img(end,:,:)];
h_input = [img(:,1,:) - img(:,end,:),diff(img,1,2)];
v_input = [img(1,:,:) - img(end,:,:);diff(img,1,1)];

Sx1=[img(:,end,:), img(:,1:end-1,:)];%ÓÒÒÆ1
S1=[Sx1(end,:,:); Sx1(1:end-1,:,:)];%ÏÂÒÆ1£¬×óÉÏ½Ç
    
Sx2=[img(:,2:end,:), img(:,1,:)];%×óÒÆ1
S2=[Sx2(end,:,:); Sx2(1:end-1,:,:)];%ÏÂÒÆ1£¬ÓÒÉÏ½Ç

l_input=img-S1;%Óë×óÉÏ½ÇÏñËØ²î
r_input=img-S2;%ÓëÓÒÉÏ½ÇÏñËØ²î

% [h_input,v_input]=myEdge(S1,3);

tt1 = (abs(h_input)>=sigma1);
tt2 = (abs(v_input)>=sigma1);
tt3 = (abs(l_input)>=sigma1);
tt4 = (abs(r_input)>=sigma1);

hh_input = h_input(tt1);
hv_input = v_input(tt2);
ll_input = l_input(tt3);
rr_input = r_input(tt4);

h_input = sign(h_input)*sigma1.*(abs(h_input)/sigma1).^alpha1;
v_input = sign(v_input)*sigma1.*(abs(v_input)/sigma1).^alpha1;
l_input = sign(l_input)*sigma1.*(abs(l_input)/sigma1).^alpha1;
r_input = sign(r_input)*sigma1.*(abs(r_input)/sigma1).^alpha1;

h_input(tt1) = hh_input;
v_input(tt2) = hv_input;
l_input(tt3) = ll_input;
r_input(tt4) = rr_input;

% h_input(:, end, :) = img(:,1,:) - img(:,end,:);
% v_input(end, :, :) = img(1,:,:) - img(end,:,:);

end