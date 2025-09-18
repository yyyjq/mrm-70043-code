% Parameter initialization & variable preparation for phase optimization
clear;
close all;

%% 1. field/FOV/sample/phantom setting

    % ===== FOV =====
    % center position [mm]
    xc = 30;
%     xc = xc_list(xc_idx);
    yc = 20; % a portion

    rs = 128;         %resolution  200 mm for human head length, so rs = 100 for resolution = 2mm

    l = 29; %length (mm) 
    w = 29; %width (mm)

    thre_k = 8; 


    step = l/rs;

    % ===== B1 =====
    Coil_N = 1;
    Coilsurf = 0; % 1: use Bio-savart law  0: uniform sensivity map
    %when Coilsurf = 1
    Coilautoposition = 1; % 1: generate the data by algrithm  0: basing on experiment data
    %wehn Coilautoposition = 0
    Coil_r = 80e-3; %m
    Coil_p(:,1)= [0;0;0;0;0;0;0;0]*1e-3; % x m
    Coil_p(:,2)= [0;0;0;0;0;0;0;0]*1e-3; % y m
    Coil_p(:,3)= [0;0;0;0;0;0;0;0]*1e-3; % z m
    Coil_Azi = [0 0 0 0 0 0 0 0]; % coil angle
    Coil_currentDir = [1 1 1 1 1 1 1 1];% 

    % ===== time convention =====
    time_con = 2;                   % 1: exp(-i omega t); 2: exp(i omega t)

    % ===== sample parameters =====
    EP_sample_N = 200/2;
    EP_angle_N = 144;
    EP_angle_Total = (2*pi + 3*pi/180 - pi/180*2.5); % middle 3
    EP_dt = (1e-6); % s
    RxDelay = 50*1e-6;%46*1e-6; % s
    T2 = 0;                         % 0 - no delay; 1 - custom delay mode
    T2_coeff = 1;                 % T2 = 0.1*T2_coeff ms
    key_T2 = 'T2_0p1ms_NoDelay';    % delay data key (only in T2 = 1 mode)
    EP_angle_start= 0;           
    % rotate direction of the magnet array
    a = 1;                          %a: 1 clockwise -1 counterclockwise
    EP_ratio = 1;
    EPsimulation = 0;               % 1: generate the data by algrithm  0: by people

    num_sample = EP_sample_N;

    % ===== Noise =====
    % simulated signal
    noise = 1;
    SNR_S = 20;                     % Noise level (dB)


    % ===== plot signal =====
    plot_k = 1;             % Whether to plot local k-space
    plot_map = 0;
    k_dim_x = 5;              % k-space dimension
    k_dim_y = 3; 

    % ===== reconstruction =====
    reCon = 1;              % Whether to do image reconstruction
    cgt=10;


    %% 2. Main parameter generator
    %% ===== FOV =====
    FOV_xs = -w/2+xc+step/2; % x starting point
    FOV_xe = w/2+xc-step/2; % x ending point
    FOV_ys = -l/2+yc+step/2; % y starting point
    FOV_ye = l/2+yc-step/2; % y ending point

    [FOV_X,FOV_Y] = meshgrid(FOV_xs:step:FOV_xe,FOV_ys:step:FOV_ye);
    FOV_x = reshape(FOV_X, numel(FOV_X), 1);
    FOV_y = reshape(FOV_Y, numel(FOV_Y), 1);

    %% ===== B0 map =====
    B0_filename = 'B0_measured_ploy55_Junqi_2022_11_2.mat';%'B0_meas_Halb.mat';%'B0_meas_Halb.mat';%
    load(B0_filename);
    B0_meas_Halb = B0_meas_Halb_poly55;

    B0_ix = reshape(B0_meas_Halb(:,:,1), numel(B0_meas_Halb(:,:,1)), 1);
    B0_iy = reshape(B0_meas_Halb(:,:,2), numel(B0_meas_Halb(:,:,2)), 1);
    B0_im = reshape(B0_meas_Halb(:,:,3).*1e-3, numel(B0_meas_Halb(:,:,3)), 1);

    B0fovi = griddata(B0_ix, B0_iy, B0_im, FOV_x, FOV_y);
    B0FOV = reshape(B0fovi, size(FOV_X,1), size(FOV_X,2));
    B0fovMax = max(B0_im(~isnan(B0_im)));
    B0fovMin = min(B0_im(~isnan(B0_im)));

    EP_f_center = mean(B0_im(~isnan(B0_im)))*42.577;
    B0fov = reshape(B0FOV, numel(FOV_X), 1);


    %% ===== Sample parameters =====
    if T2 == 0
        key_T2 = 'NoT2_NoDelay';
        delayCorrArray = zeros(1,EP_angle_N);
        T2star = zeros(EP_angle_N,1);
        selectDelay = 0;
    elseif T2 == 1
        T2star = ones(1,EP_angle_N)*1e-4*T2_coeff;
        delayCorrArray = zeros(1,EP_angle_N);
        selectDelay = 1;
    end

    E_fshift = 0;
    E_xcm = 0; % eccentric error
    E_ycm = 0; % eccentric error 
    E_tra = zeros(EP_angle_N,2);

    %% ===== Load phantom =====

    load('phantom_13bottle_up_new.mat');
    SIM_M = imresize(image_exp,[rs,rs]);
    SIM_M = flipud(SIM_M);
    SIM_m = reshape(SIM_M,numel(SIM_M),1);

    % SIM_M = zeros(rs);
    % % SIM_M(5,5) = 1;
    % SIM_M(ceil(rs/2),ceil(rs/2)) = 1;
    % SIM_m = reshape(SIM_M,numel(SIM_M),1);

    figure(2);
    h = pcolor(FOV_X,FOV_Y,SIM_M);
    shading flat;
    axis square;
%     colormap gray;
    title('Original phantom');
    set(gcf,'color','w');
    pbaspect([1 1 1]);

    % saveas(h,'phantom_C.jpg');

    %% ===== Build EB0_long as the base =====
    const_gamma = 42580000;

    fshift = E_fshift;
    xcm = E_xcm;
    ycm = E_ycm;

    idx_valid = ~isnan(B0_im);
    B0fit = fit([B0_ix(idx_valid), B0_iy(idx_valid)],B0_im(idx_valid),'poly55');

    % Add a constant B0 fit for main field
    EB0x = FOV_x-xcm;
    EB0y = FOV_y-ycm;

    addOnPoint = round(10*RxDelay/EP_dt);
    sampleLong = EP_sample_N+addOnPoint;
    EB0_long_all = zeros(numel(FOV_x),sampleLong,EP_angle_N);
    Gbx = zeros(size(FOV_X,1),size(FOV_X,2),EP_angle_N);
    Gby = zeros(size(FOV_X,1),size(FOV_X,2),EP_angle_N);
    for i_an = 1:EP_angle_N
        % ===== Produce the B0 maps =====
        R_angle = EP_angle_start+(i_an-1)*EP_angle_Total./(EP_angle_N);

        % rotate direction of the magnet array
        % (para)a = 1; %a: 1 clockwise -1 counterclockwise
        R_angle = a*R_angle;

        % calculate the postion
        EB0x = (FOV_x-xcm).*cos(R_angle) - (FOV_y-ycm).*sin(R_angle);
        EB0y = (FOV_x-xcm).*sin(R_angle) + (FOV_y-ycm).*cos(R_angle);

        EB0_org = B0fit.p00 + B0fit.p10.*EB0x + B0fit.p01.*EB0y + B0fit.p20.*EB0x.^2 + B0fit.p11.*EB0x.*EB0y + B0fit.p02.*EB0y.^2 + B0fit.p30.*EB0x.^3......
            + B0fit.p21.*EB0x.^2.*EB0y + B0fit.p12.*EB0x.*EB0y.^2 + B0fit.p03.*EB0y.^3 + B0fit.p40.*EB0x.^4 + B0fit.p31.*EB0x.^3.*EB0y......
            + B0fit.p22.*EB0x.^2.*EB0y.^2 + B0fit.p13.*EB0x.*EB0y.^3 + B0fit.p04.*EB0y.^4 + B0fit.p50.*EB0x.^5 + B0fit.p41.*EB0x.^4.*EB0y......
            + B0fit.p32.*EB0x.^3.*EB0y.^2 + B0fit.p23.*EB0x.^2.*EB0y.^3 + B0fit.p14.*EB0x.*EB0y.^4 + B0fit.p05.*EB0y.^5;


        % define the sensitivity map
            if Coilsurf == 1
                EB1 = EB1_x-1i*EB1_y;
            else
                EB1 = ones(numel(FOV_X),Coil_N).*1e-4;
%                 EB1 = EB1_z;
            end

        f=EB0_org*const_gamma-fshift;

        EB0_long = zeros(numel(f),sampleLong);
        for i_sample = 1:sampleLong
            if time_con ==1
                EB0_long(:,i_sample) = exp(-1i*2*pi*(f-EP_f_center)*(i_sample-1)*EP_dt);
            elseif time_con ==2
                EB0_long(:,i_sample) = exp(1i*2*pi*(f-EP_f_center)*(i_sample-1)*EP_dt);
            end        
        end
        EB0_long_all(:,:,i_an) = EB0_long;
        EB0M = reshape(EB0_org,size(FOV_X,1),size(FOV_X,2));       
        [Gbx(:,:,i_an),Gby(:,:,i_an)] = gradient(EB0M,step);

        if plot_map == 1
            figure(3);
            h=pcolor(FOV_X,FOV_Y,EB0M); 
            xlabel('x (mm)');
            ylabel('y (mm)');
            title(strcat('B0 at angle ',num2str(i_an),' of ',num2str(EP_angle_N),' angles'));
            shading flat; 
            axis square;
            colorbar;
            hold on;
            contour(FOV_X,FOV_Y,EB0M,20,'red');
            drawnow;
            clim([B0fovMin B0fovMax]);
        end
    end

    %% Plot local k-space
    idx = 1;
    randnum = randperm(144);
    NRMSE_SSIM_PSNR = zeros(3,144);
    for k = 144
        G_angle = [];
        G_abs = [];
%     data_left = randnum;
    data_left = [130,133,127,136,124,139,121,69,118,72,75,115,78,112,81,109,84,141,106,87,103,90,131,135,128,100,125,138,93,122,97,71,119,95,74,116,113,110,108,105,77,102,89,99,85,91,80,83,107,98,101,104,111,92,96,94,67,114,88,86,82,79,117,76,143,120,73,123,65,70,126,129,1,68,63,134,132,66,137,140,64,2,142,62,144,3,61,4,60,59,58,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57];
    data_left = data_left(1:k);
    data_trancate = setdiff(1:EP_angle_N,data_left); %input data_left
    sim_signal_idx_list = setdiff(1:EP_angle_N,[data_trancate]);
    num_angle = length(sim_signal_idx_list);

    if plot_k == 1

        GXK = Gbx*const_gamma*EP_dt;
        GYK = Gby*const_gamma*EP_dt;

        % plot local k space by useing GX,GY
        h2=figure(23);
        set(gcf, 'position', [0 0 500 500], 'color', 'w');
        [dM,dN] = size(FOV_X);
        dxstep = fix(dM/k_dim_x);dystep = fix(dN/k_dim_y);
        for ix = 1:k_dim_x
            for iy = 1:k_dim_y
                h = subplot(k_dim_x,k_dim_y,ix*k_dim_y+iy-k_dim_y);
                if k_dim_x < 10
                    set(h,'position',[0.11*iy 1-0.11*ix 0.1 0.1]);
                else
                    set(h,'position',[0.06*(iy-0.5) 1-0.06*(ix+0.1) 0.05 0.05]);
                end

                ix_data = k_dim_x + 1 - ix; % k-space revise
                for it = 1:length(sim_signal_idx_list)
                    it_real = sim_signal_idx_list(it);
                    mx = GXK(fix(dM/2)+(ix_data-ceil(k_dim_x/2))*dxstep,fix(dN/2)+(iy-ceil(k_dim_y/2))*dystep,it_real);
                    my = GYK(fix(dM/2)+(ix_data-ceil(k_dim_x/2))*dxstep,fix(dN/2)+(iy-ceil(k_dim_y/2))*dystep,it_real);

                    G_angle(it_real) = angle(mx + my*1i)/(2*pi)*360;
                    G_abs(it_real) = sqrt(mx.^2+my.^2);

                    n = (0:EP_sample_N/10:EP_sample_N); 
                    plot(n*mx,n*my,'.r','MarkerSize',2);
                    hold on;
                end
                k_bound = rs/(2*w)/7; % rs/w; %  1/del_w
                axis([-k_bound k_bound -k_bound k_bound]); pbaspect([1 1 1]);
                set(gca,'ytick',[]);
                set(gca,'xtick',[]);
            end
        end
    end

    % Prepare the simulated signal
    G_angle = wrapTo360(G_angle - G_angle(1)); % Initial is 360 or 0
    k_bound_real = k_bound;
    idx = idx +1;

    if plot_k
        Signalt = zeros(EP_sample_N,Coil_N,num_angle);
        EB0_long_size = size(EB0_long_all);
        EB0_long_size = EB0_long_size(1:2);
        Et = zeros(EB0_long_size(1),EP_sample_N,Coil_N,num_angle);
        for i_an = 1:length(sim_signal_idx_list)
            i_real = sim_signal_idx_list(i_an);
            EB0_long = reshape(EB0_long_all(:,:,i_real),EB0_long_size);
            delayCorr = delayCorrArray(i_real);
            sampleStart = RxDelay/EP_dt+delayCorr;
            EB0 = EB0_long(:,sampleStart+1:sampleStart+EP_sample_N);

            for i_coil = 1:Coil_N
                Ea1 = EB0(:,:);
                Ea2 = repmat(EB1(:,i_coil),1,EP_sample_N);
                Ep1 = Ea1.*Ea2;
                s = (Ep1.')*SIM_m;
                if noise == 0
                    % NO noise
                    sn = s;
                else
                    % add noise
                    sn = awgn(s,SNR_S,'measured');
                end
                % generate ideal signal in time domain
                Signalt(:,i_coil,i_real) = sn;
                % geneted ideal encoding matrix in time domain
                Et(:,:,i_coil,i_real) = Ep1;
            end
        end
    end

    % Use optimized encoding matrix & measured signal
        if reCon == 1
            num_angle = length(sim_signal_idx_list);
            Et_final = Et(:,1:num_sample,:,:);
            At = reshape(Et(:,:,:,sim_signal_idx_list),rs*rs,EP_sample_N*num_angle*Coil_N);
            bt = reshape(Signalt(:,sim_signal_idx_list),EP_sample_N*num_angle*Coil_N,1);
            At = At.';
    
            Atcg = At'*At;
            btcg = At'*bt;
    
            for i=cgt
                CGtimes = i;
                [x2,fkg,xres2] = lsqr(Atcg,btcg,1e-10,CGtimes);
                f3=figure(300+i);
                [imageRecg,STR_Icg,NUM_Icg] = funRearrangeimage(abs(x2),SIM_M); %NUM_I=[NRMSE_I SSIM_I PSNR_I];
                pcolor(FOV_X,FOV_Y,imageRecg);
                shading flat;
    %             colormap gray;
                pbaspect([1 1 1])
                set(gcf,'color','w');
            end
        end
    end
