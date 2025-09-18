function [imageRe,STR_I,NUM_I,SSIM_map] = funRearrangeimage(x,SIM_M)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
x=x-min(min(x));
x=x./max(max(x));
imageRe=reshape(x,size(SIM_M,1),size(SIM_M,2));
R=funOtherIMquality(SIM_M,imageRe);
[SSIM_I,SSIM_map]=ssim(imageRe,SIM_M,"Exponents",[1,1,1]);
PSNR_I=R.PSNR;
NRMSE_I=R.NRMSE;
STR_I=strcat('NRMSE',num2str(NRMSE_I),',SSIM:',num2str(SSIM_I),',PSNR:',num2str(PSNR_I));
NUM_I=[NRMSE_I SSIM_I PSNR_I];
end

