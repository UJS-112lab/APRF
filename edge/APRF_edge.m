function x = APRF_edge(S,lambda,alpha,q,omega,omegaq,sigma1,alpha1,rho)
%% Global constants and defaults

MAX_ITER =10;

%% Data preprocessing
S = gpuArray(double(S));% Transfer data to GPU


[h_input,v_input,l_input,r_input] = MyMapping(S,sigma1,alpha1);
[row, col, cha] = size(S);


Ss = S;

h_input1 = [Ss(:,1,:) - Ss(:,end,:),diff(Ss,1,2),];
v_input1 = [Ss(1,:,:) - Ss(end,:,:);diff(Ss,1,1) ];
[Ex,~]=max(abs(h_input1),[],3);
[Ey,~]=max(abs(v_input1),[],3);

wx = exp(-(Ex.^2)./(omega^2));
wy = exp(-(Ey.^2)./(omega^2));


% 
Sx1=[Ss(:,end,:), Ss(:,1:end-1,:)];%右移1
S1=[Sx1(end,:,:); Sx1(1:end-1,:,:)];%下移1，左上角



    
Sx2=[Ss(:,2:end,:), Ss(:,1,:)];%左移1
S2=[Sx2(end,:,:); Sx2(1:end-1,:,:)];%下移1，右上角



G1=Ss-S1;%与左上角像素差
G2=Ss-S2;%与右上角像素差
[E1,~]=max(abs(G1),[],3);
[E2,~]=max(abs(G2),[],3);

wl = exp(-(E1.^2)./(omega^2));%左上角权重
wr = exp(-(E2.^2)./(omega^2));%右上角权重

wxq = exp(-(Ex.^2)./(omegaq^2));
wyq = exp(-(Ey.^2)./(omegaq^2));
wlq = exp(-(E1.^2)./(omegaq^2));%左上角权重
wrq = exp(-(E2.^2)./(omegaq^2));%右上角权重



w1 = wx.* lambda/rho;
w2 = wy.* lambda/rho;
w3 = wl.* lambda/rho;
w4 = wr.* lambda/rho;

qx = q./wxq;
qy = q./wyq;
ql = q./wlq;
qr = q./wrq;


%% ADMM solver
x = gpuArray.zeros(row, col, cha); % Initialize x on the GPU
z1 = x; % Initialize z1 on the GPU
z2 = x; % Initialize z2 on the GPU
z3 = x; % Initialize z1 on the GPU
z4 = x; % Initialize z2 on the GPU
u1 = x; % Initialize u1 on the GPU
u2 = x; % Initialize u2 on the GPU
u3 = x; % Initialize u1 on the GPU
u4 = x; % Initialize u2 on the GPU

size2D = [row,col];

otfFx1 = psf2otf_Dx_GPU(size2D); % equal to otfFx = psf2otf(fx, sizeI2D) where fx = [-1, 1];
otfFy1 = psf2otf_Dy_GPU(size2D); % equal to otfFy = psf2otf(fy, sizeI2D) where fy = [-1; 1];
otfF1=gpuArray(psf2otf([1 0; 0 -1], size2D));
otfF2=gpuArray(psf2otf([0 1; -1 0], size2D));
Denormin = abs(otfFx1).^2  +abs(otfFy1).^2 + abs(otfF1).^2 + abs(otfF2).^2 ;
Denormin = repmat(Denormin, [1, 1, cha]);
Denormin = fft2(1)+rho/2*Denormin;
x = S ;
Normin = fft2(x);


for k = 1:MAX_ITER
   
   
    % x-update
    h = h_input + z1 - u1;
    v = v_input + z2 - u2;
    l = l_input + z3 - u3;
    r = r_input + z4 - u4;

    Norminh = [-diff(h,1,2),h(:,end,:) - h(:, 1,:)];
    Norminv = [-diff(v,1,1);v(end,:,:) - v(1, :,:)];
    
    l1=[l(:,2:end,:), l(:,1,:)];%左移1
    l1=[l1(2:end,:,:); l1(1,:,:)];%上移1，右下角
    
    r1=[r(:,end,:), r(:,1:end-1,:)];%右移1
    r1=[r1(2:end,:,:); r1(1,:,:)];%上移1，左下角
    
    Norminl=l-l1;
    Norminr=r-r1;

    Normin2 = Norminh + Norminv +Norminl +Norminr;
    %Normin3 = b_u;
    FS = (Normin+rho/2*fft2(Normin2))./Denormin;
    x = real(ifft2(FS)); 
    %Normin = FS;
    

    
    % z-update

    h1 = [x(:,1,:) - x(:,end,:),diff(x,1,2)];
    v1 = [x(1,:,:) - x(end,:,:);diff(x,1,1);];
    
    Ux1=[x(:,end,:), x(:,1:end-1,:)];%右移1
    U1=[Ux1(end,:,:); Ux1(1:end-1,:,:)];%下移1，左上角
    
    Ux2=[x(:,2:end,:), x(:,1,:)];%左移1
    U2=[Ux2(end,:,:); Ux2(1:end-1,:,:)];%下移1，右上角
    
    l1=x-U1;
    r1=x-U2; 
    
    Ax_hatx = alpha .* (h1 - h_input) + (1-alpha) .* z1 + u1;
    Ax_haty = alpha .* (v1 - v_input) + (1-alpha) .* z2 + u2;
    Ax_hatl = alpha .* (l1 - l_input) + (1-alpha) .* z3 + u3;
    Ax_hatr = alpha .* (r1 - r_input) + (1-alpha) .* z4 + u4;
    
    %solve with p-shrinkage
    z1 = max(abs(Ax_hatx) - w1.*(abs(Ax_hatx).^(qx-1)), 0) .* sign(Ax_hatx);
    z2 = max(abs(Ax_haty) - w2.*(abs(Ax_haty).^(qy-1)), 0) .* sign(Ax_haty);
    z3 = max(abs(Ax_hatl) - w3.*(abs(Ax_hatl).^(ql-1)), 0) .* sign(Ax_hatl);
    z4 = max(abs(Ax_hatr) - w4.*(abs(Ax_hatr).^(qr-1)), 0) .* sign(Ax_hatr);


    % u-update
    u1 = Ax_hatx - z1;
    u2 = Ax_haty - z2;
    u3 = Ax_hatl - z3;
    u4 = Ax_hatr - z4;

end
x = gather(x); % Transfer x back to the CPU
end
