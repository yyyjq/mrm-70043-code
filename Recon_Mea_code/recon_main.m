% Direct image reconstruction
% Optimize the delay
warning('off','all');
clear all;
close all;

recon_init;

%% Load 13-bottle phantom
load('phantom\phantom_13bottle_up_new.mat');
SIM_M = imresize(image_exp,[rs,rs]);
SIM_M = flipud(SIM_M);
SIM_m = reshape(SIM_M,numel(SIM_M),1);

figure(2);
h = pcolor(FOV_X,FOV_Y,SIM_M);
shading flat;
axis square;
% colormap gray;
title('Original phantom');
set(gcf,'color','w');
pbaspect([1 1 1]);

%delay_op = zeros(1,EP_angle_N);
%T2_op = ones(1,EP_angle_N)*10; % 0.1 ms

%% 3. Save Result

%[Signalt,Et] = CalSignal(SIM_m,para_list,EB0_long_all,EB1,zeros(1,EP_angle_N));
EB0_long_size = size(EB0_long_all);
EB0_long_size = EB0_long_size(1:2);
for i_an = 1:EP_angle_N
    EB0_long = reshape(EB0_long_all(:,:,i_an),EB0_long_size);
%     if phase_en
%         phase_factor = exp(1i*delay_op(i_an)*EP_dt*weight*delta_B_all(i_an)*const_gamma);
%     else
%         phase_factor = 1;
%     end
    for i_sample = 1:sampleLong
        if T2star(i_an) == 0
            T2factor = 1;
        else
            T2factor = exp(-(i_sample-1)*EP_dt/T2star(i_an));
        end
        EB0_long(:,i_sample) = EB0_long(:,i_sample)*T2factor;
    end
    EB0_long_all(:,:,i_an) = EB0_long;
end
EB0_long_size = size(EB0_long_all);
EB0_long_size = EB0_long_size(1:2);
Et = zeros(EB0_long_size(1),EP_sample_N,Coil_N,EP_angle_N);
Signalt = zeros(EP_sample_N,Coil_N,EP_angle_N);
for i_an = 1:EP_angle_N
    EB0_long = reshape(EB0_long_all(:,:,i_an),EB0_long_size);
    delayCorr = delayCorrArray(i_an);
    sampleStart = RxDelay/EP_dt+delayCorr;
    EB0 = EB0_long(:,sampleStart+1:sampleStart+EP_sample_N);

    for i_coil = 1:Coil_N
        Ea1 = EB0(:,:);
        Ea2 = repmat(EB1(:,i_coil),1,EP_sample_N);
        Ep1 = Ea1.*Ea2;
        % geneted ideal encoding matrix in time domain
        Et(:,:,i_coil,i_an) = Ep1;
        % Generate simulated signal
        s = (Ep1.')*SIM_m;
        Signalt(:,i_coil,i_an) = s;
    end
end


%% 5. Image Reconstruction (Optional)
for k = 34
    if reCon == 1
        % Only use the indices with similar signal
        % normal seq
    %          data_all = 1:4:144;
        % opt seq
         data_all = [130,133,127,136,124,139,121,69,118,72,75,115,78,112,81,109,84,141,106,87,103,90,131,135,128,100,125,138,93,122,97,71,119,95,74,116,113,110,108,105,77,102,89,99,85,91,80,83,107,98,101,104,111,92,96,94,67,114,88,86,82,79,117,76,143,120,73,123,65,70,126,129,1,68,63,134,132,66,137,140,64,2,142,62,144,3,61,4,60,59,58,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56
    ];
        data_left = data_all(1:0+k);
        data_trancate = setdiff(1:EP_angle_N,data_left); %input data_left
        sim_signal_idx_list = setdiff(1:EP_angle_N,[data_trancate]);
        num_angle = length(sim_signal_idx_list);
        At = reshape(Et(:,:,:,sim_signal_idx_list),rs*rs,EP_sample_N*num_angle*Coil_N);
        bt = reshape(Signalfilt(:,sim_signal_idx_list),EP_sample_N*num_angle*Coil_N,1); % Measurement data
    
        At = At.';
        
        Atcg = At'*At;
        btcg = At'*bt;
        
        for i=2
            CGtimes = i;
            [x2,fkg,xres2] = lsqr(Atcg,btcg,1e-10,CGtimes);
            f3=figure(300+i);
            [imageRecg,STR_Icg,NUM_Icg] = funRearrangeimage(abs(x2),SIM_M); %NUM_I=[NRMSE_I SSIM_I PSNR_I];
            pcolor(FOV_X,FOV_Y,imageRecg);
            shading flat;
            axis square;
            pbaspect([1 1 1])
            curr_fig = figure(300+i);
            filename = strcat('SpokeNum_',num2str(k),'.jpg');    
        end
    end
end