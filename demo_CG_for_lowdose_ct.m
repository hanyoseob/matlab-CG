%%  REFERENCE
% https://en.wikipedia.org/wiki/Conjugate_gradient_method

%% COST FUNCTION
% x^* = argmin_x { 1/2 * || A(X) - Y ||_2^2 + lambda/2 * ( || D_x(X) ||_2^2 + || D_y(X) ||_2^2 ) }

%%
clear ;
close all;
home;

%% GPU Processing
% If there is GPU device on your board, 
% then isgpu is true. Otherwise, it is false.
bgpu    = false;
bfig    = true;

%%  SYSTEM SETTING
N       = 512;
VIEW    = 360;
THETA   = linspace(0, 180, VIEW + 1);   THETA(end) = [];

R       = @(x) radon(x, THETA);
RT      = @(y) iradon(y, THETA, 'none', N)/(pi/(2*length(THETA)));
RINV    = @(y) iradon(y, THETA, N);

%% DATA GENERATION
load('XCAT512.mat');
x       = imresize(double(XCAT512), [N, N]);
p       = R(x);
x_full  = RINV(p);

%% LOW-DOSE SINOGRAM GENERATION
i0     	= 5e4;
pn     	= exp(-p);
pn     	= i0.*pn;
pn     	= poissrnd(pn);
pn      = max(-log(max(pn,1)./i0),0);

y       = pn;
x_low   = RINV(y);

%% CONJUGATE GRADIENT METHOD (CG) INITIALIZATION
LAMBDA  = 5e2;
A       = @(x) (RT(R(x))  + LAMBDA*(Dxt(Dx(x)) + Dyt(Dy(x))));

b0      = RT(y);
x0      = zeros(N);

niter   = 2e1;

L2              = @(x) power(norm(x, 'fro'), 2);
COST.equation   = '1/2 * || A(X) - Y ||_2^2 + lambda/2 * ( || D_x(X) ||_2^2 + || D_y(X) ||_2^2 )';
COST.function	= @(x) 1/2 * L2(R(x) - y) + LAMBDA/2 * (L2(Dx(x)) + L2(Dy(x)));

%% RUN CONJUGATE GRADIENT METHOD (CG)
if bgpu
    b0 = gpuArray(b0);
    x0 = gpuArray(x0);
end

[x_cg, obj]	= CG(A, b0, x0, niter, COST, bfig);

%% CALCUATE QUANTIFICATION FACTOR 
x_low       = max(x_low, 0);
x_cg        = max(x_cg, 0);
nor         = max(x(:));

mse_x_low   = immse(x_low./nor, x./nor);
mse_x_cg    = immse(x_cg./nor, x./nor);

psnr_x_low 	= psnr(x_low./nor, x./nor);
psnr_x_cg 	= psnr(x_cg./nor, x./nor);

ssim_x_low  = ssim(x_low./nor, x./nor);
ssim_x_cg   = ssim(x_cg./nor, x./nor);

%% DISPLAY
wndImg  = [0, 0.03];

figure('name', 'Conjugate Gradient (CG) Method');
colormap(gray(256));

suptitle('Conjugate Gradient (CG) Method');
subplot(231);   imagesc(x,     	wndImg); 	axis image off;     title('ground truth');
subplot(232);   imagesc(x_full, wndImg);   	axis image off;     title(['full-dose_{FBP, view : ', num2str(VIEW) '}']);
subplot(234);   imagesc(x_low,  wndImg);   	axis image off;     title({['low-dose_{FBP, view : ', num2str(VIEW) '}'], ['MSE : ' num2str(mse_x_low, '%.4e')], ['PSNR : ' num2str(psnr_x_low, '%.4f')], ['SSIM : ' num2str(ssim_x_low, '%.4f')]});
subplot(235);   imagesc(x_cg,   wndImg);  	axis image off;     title({['recon_{CG}'], ['MSE : ' num2str(mse_x_cg, '%.4e')], ['PSNR : ' num2str(psnr_x_cg, '%.4f')], ['SSIM : ' num2str(ssim_x_cg, '%.4f')]});

subplot(2,3,[3,6]); semilogy(obj, '*-');    title(COST.equation);  xlabel('# of iteration');   ylabel('Objective'); 
                                            xlim([1, niter]);   grid on; grid minor;

