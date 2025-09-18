% Parameter initialization & variable preparation for phase optimization
clear;
close all;

%% 1. field/FOV/sample/phantom setting

% ===== B0 map =====
B0_map = 1;                     % 1 - Halbach, 2 - IO-Ring, 3 - Constant main + rotated G

xc_list = 30;
xc_len = length(xc_list);
NRMSE_SSIM_PSNR_all = zeros(3,xc_len);

k_dim_x = 5;              % k-space dimension
k_dim_y = 3; 
thre_k = 8;                 %degree

Eval = zeros(length(xc_list),k_dim_x*k_dim_y);
Angle_all = zeros(2,k_dim_x*k_dim_y,144);
G_abs_all = zeros(length(xc_list),k_dim_x*k_dim_y);

for xc_idx = 1:xc_len

% ===== FOV =====
% center position [mm]
% xc = 30;
xc = xc_list(xc_idx);
yc = 19; 

rs = 128;         %resolution  200 mm for human head length, so rs = 100 for resolution = 2mm

l = 29; %78.45*2; %length (mm)
w = 29; %78.45*2; %width (mm)

step = l/rs;

% ===== time convention =====
time_con = 2;                   % 1: exp(-i omega t); 2: exp(i omega t)

% ===== sample parameters =====
EP_sample_N = 260;
EP_angle_N = 144;
EP_angle_Total = (2*pi + 3*pi/180 - pi/180*2.5); % middle 3
EP_dt = 0.5e-6; % s
RxDelay = 0;%46*1e-6; % s
T2 = 0;                         % 0 - no delay; 1 - custom delay mode
T2_coeff = 1;                 % T2 = 0.1*T2_coeff ms
key_T2 = 'T2_0p1ms_NoDelay';    % delay data key (only in T2 = 1 mode)
EP_f_center = 0;           %gj: 2.84MHz; yd1:2.84MHz; yd: 2.85MHz; small tube: 2.85; big tube: 2.84
EP_angle_start= 0;           % yidan: pi; gj: pi
% rotate direction of the magnet array
a = 1;                          %a: 1 clockwise -1 counterclockwise
EP_ratio = 1;
EPsimulation = 0;               % 1: generate the data by algrithm  0: by people

num_sample = EP_sample_N;

% ===== Noise =====
% simulated signal
noise = 0;
SNR_S = 15;                     % Noise level (dB)

if noise == 1
    key_noise = strcat('noisy_',num2str(SNR_S),'dB');
else
    key_noise = 'noiseless';
end


% ===== encoding matrix =====
B0_D = [0,1,0];

filekey = 'fullZshimlowG30PSF';
version = strcat('x0y0mm_size',num2str(floor(l)),'mm_rs',num2str(rs),'_',filekey);
range = strcat('x0y0mm_size',num2str(floor(l)),'mm_rs',num2str(rs),'_',filekey);

sigS = -1;

% ===== plot signal =====
plotSig = 1;            % Whether to plot the signal after optimization
phase_en = 0;           % Whether to add the phase correction
plot_k = 1;             % Whether to plot local k-space
plot_eval = 0;



% ===== reconstruction =====
reCon = 0;              % Whether to do image reconstruction
cgt=10;

delay = RxDelay*1e6;
key_angle = 'allAngles';

%% 2. Main parameter generator
%% ===== FOV =====
FOV_xs = -w/2+xc+step/2; % x starting point
FOV_xe = w/2+xc-step/2; % x ending point
FOV_ys = -l/2+yc+step/2; % y starting point
FOV_ye = l/2+yc-step/2; % y ending point

% FOV_X, FOV_Y are coordinate matrices
[FOV_X,FOV_Y] = meshgrid(FOV_xs:step:FOV_xe,FOV_ys:step:FOV_ye);
% FOV_x, FOV_y are coordinate arrays
FOV_x = reshape(FOV_X, numel(FOV_X), 1);
FOV_y = reshape(FOV_Y, numel(FOV_Y), 1);

% D is a scaling factor (do not know the function)
D=1;
[FOVD_X,FOVD_Y] = meshgrid(FOV_xs-0.5*(step-D*step):D*step:...
    FOV_xe+0.5*(step-D*step),FOV_ys-0.5*(step-D*step)...
    :D*step:FOV_ye+0.5*(step-D*step));
FOVD_x = reshape(FOVD_X, numel(FOVD_X), 1);
FOVD_y = reshape(FOVD_Y, numel(FOVD_Y), 1);

%% ===== B0 map =====
addpath ('C:\Users\snysb\Desktop\MATLAB\Magnetic field');
if B0_map == 1
    B0_filename = 'B0_measured_ploy55_Junqi_2022_11_2.mat';
    load(B0_filename);
    B0_meas_Halb = B0_meas_Halb_poly55;
    
    B0_ix = reshape(B0_meas_Halb(:,:,1), numel(B0_meas_Halb(:,:,1)), 1);
    B0_iy = reshape(B0_meas_Halb(:,:,2), numel(B0_meas_Halb(:,:,2)), 1);
    B0_im = reshape(B0_meas_Halb(:,:,3).*1e-3, numel(B0_meas_Halb(:,:,3)), 1);
%     B0_im = reshape(B0_meas_Halb(:,:,3).*5e-3, numel(B0_meas_Halb(:,:,3)), 1);
    
    B0fovi = griddata(B0_ix, B0_iy, B0_im, FOV_x, FOV_y);
    B0FOV = reshape(B0fovi, size(FOV_X,1), size(FOV_X,2));
    B0fovMax = max(B0_im(~isnan(B0_im)));
    B0fovMin = min(B0_im(~isnan(B0_im)));
elseif B0_map == 2
    B0_filename = 'Quadrupolar_main66p3_delta10.mat'; % 'Linear_main66p3_delta12p8.mat';%'IO_Bz_G35_mean100_large270_len480_base_slices_z60.mat';%'Linear_main111p5_delta8p2.mat';%'IO_Bz_G35_large260_custom_z20.mat';%'Linear_main111p4_delta5p91.mat';%'Halbach_main109p70_delta6p33.mat';%'IO_Bz_G25_pop100_gen100_Ng4_nbar10.mat';%'IO_Bz_latest_G25_new_Ng4_nbar10.mat';%'Linear_main109p70_delta6p33.mat';%'IO_Bz_A_lowG30_new.mat';%'Halbach_main111p33_delta6p52.mat';%'IO_Bz_fullZshim_0.mat';%'IO_field_main111p33_grad3.mat';%'IO_Bz_optimized_A_3D_60.mat';
    load(B0_filename);
    B0_im = reshape(Bz*1e-3,numel(Bz),1);
%     B0_im = reshape(Bz,numel(Bz),1);
%     B0_im = reshape((Bz)*1e-3,numel(Bz),1);
    B0_ix = reshape(x,numel(x),1);
    B0_iy = reshape(y,numel(y),1);
    
    B0_im(isnan(B0_im)) = 0;

    B0fovi = griddata(B0_ix, B0_iy, B0_im, FOV_x, FOV_y);
    B0FOV = reshape(B0fovi, size(FOV_X,1), size(FOV_X,2));
    B0_im(B0_im == 0) = NaN;
    B0fovi(B0fovi == 0) = NaN;
    B0FOV(B0FOV == 0) = NaN;
    B0fovMax = max(B0_im(~isnan(B0_im)));
    B0fovMin = min(B0_im(~isnan(B0_im)));
elseif B0_map == 3
    B0_filename = 'IO_Bz_subgrad.mat';%'IO_field_main204_grad11.mat';
    load(B0_filename);
    B0_im = reshape(Bz*1e-3,numel(Bz),1);
    B0_ix = reshape(x,numel(x),1);
    B0_iy = reshape(y,numel(y),1);
    
    % Load the unrotated main field
    load('IO_Bz_optimized_A_half_main.mat');
    B0_im_main = reshape(Bz*1e-3,numel(Bz),1);

    B0_im(isnan(B0_im)) = 0;

    B0fovi = griddata(B0_ix, B0_iy, B0_im, FOV_x, FOV_y);
    B0FOV = reshape(B0fovi, size(FOV_X,1), size(FOV_X,2));
    B0_im(B0_im == 0) = NaN;
    B0fovi(B0fovi == 0) = NaN;
    B0FOV(B0FOV == 0) = NaN;
    
    B0_sum = B0_im + B0_im_main;
    B0fovMax = max(B0_sum(~isnan(B0_sum)));
    B0fovMin = min(B0_sum(~isnan(B0_sum)));
end

B0fov = reshape(B0FOV, numel(FOV_X), 1);

%% === Coil senstivity ===
% choose the automatical coil position generation
    B1x = zeros(numel(FOV_X),1); 
    B1y = B1x; 
    B1z = ones(numel(FOV_X),1);

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
if B0_map == 3
    B0fit_main = fit([B0_ix(idx_valid), B0_iy(idx_valid)],B0_im_main(idx_valid),'poly55');
    EB0_main = B0fit_main.p00 + B0fit_main.p10.*EB0x + B0fit_main.p01.*EB0y + B0fit_main.p20.*EB0x.^2 + B0fit_main.p11.*EB0x.*EB0y + B0fit_main.p02.*EB0y.^2 + B0fit_main.p30.*EB0x.^3......
            + B0fit_main.p21.*EB0x.^2.*EB0y + B0fit_main.p12.*EB0x.*EB0y.^2 + B0fit_main.p03.*EB0y.^3 + B0fit_main.p40.*EB0x.^4 + B0fit_main.p31.*EB0x.^3.*EB0y......
            + B0fit_main.p22.*EB0x.^2.*EB0y.^2 + B0fit_main.p13.*EB0x.*EB0y.^3 + B0fit_main.p04.*EB0y.^4 + B0fit_main.p50.*EB0x.^5 + B0fit_main.p41.*EB0x.^4.*EB0y......
            + B0fit_main.p32.*EB0x.^3.*EB0y.^2 + B0fit_main.p23.*EB0x.^2.*EB0y.^3 + B0fit_main.p14.*EB0x.*EB0y.^4 + B0fit_main.p05.*EB0y.^5;
end

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

    % remove parallel component
    EB1_z = B1z;

    % define the sensitivity map
    EB1 = EB1_z;
    
    if B0_map == 3
        EB0_org = EB0_org + EB0_main;
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
    
    if plot_k == 1
        EB0M = reshape(EB0_org,size(FOV_X,1),size(FOV_X,2));
        
        [Gbx(:,:,i_an),Gby(:,:,i_an)] = gradient(EB0M,step);
        
%         figure(3);
%         h=pcolor(FOV_X,FOV_Y,EB0M); 
%         xlabel('x (mm)');
%         ylabel('y (mm)');
%         title(strcat('B0 at angle ',num2str(i_an),' of ',num2str(EP_angle_N),' angles'));
%         shading flat; 
%         axis square;
%         colorbar;
%         hold on;
%         contour(FOV_X,FOV_Y,EB0M,20,'red');
%         drawnow;
%         caxis([B0fovMin B0fovMax]);
%         
%         saveas(h,strcat('B0data_144angles\B0_',version,'_at_angle',num2str(i_an),'of',num2str(EP_angle_N),'angles.jpg'));
%         
%         Etitle = strcat('B0',version,'_',num2str(EP_angle_N),'angle',num2str(i_an),'.mat');
%         save(['B0data_144angles\',Etitle],'EB0_org');
        
        GbxGbyName = strcat('Gbx_Gby',range,'.mat');
        save(GbxGbyName,'Gbx','Gby');
    end
end

%% Plot local k-space

    index_angle = linspace(0,360,EP_angle_N);

    Gx = zeros(1,EP_angle_N);
    Gy = zeros(1,EP_angle_N);
    G_angle = zeros(1,EP_angle_N);
    G_abs = zeros(1,EP_angle_N);
    Eval_all = zeros(EP_angle_N,1,k_dim_y*k_dim_x);


mkdir('plot_local_k');
if plot_k == 1
    % define the k space boundary
    load(GbxGbyName);
    % transfer gradient to k space
    %     GXK = Gbx/step*const_gama*EP_dt*step;
    %     GYK = Gby/step*const_gama*EP_dt*step;
    G_subplot_angle = zeros(k_dim_x,k_dim_y,EP_angle_N);
    Gx_subplot_angle = zeros(k_dim_x,k_dim_y,EP_angle_N);
    Gy_subplot_angle = zeros(k_dim_x,k_dim_y,EP_angle_N);
    
    GXK = Gbx*const_gamma*EP_dt;
    GYK = Gby*const_gamma*EP_dt;
    
    % plot local k space by useing GX,GY
        h2=figure(23);
        set(gcf, 'position', [0 0 500 500], 'color', 'w');
        [dM,dN] = size(FOVD_X);
        dxstep = fix(dM/k_dim_x);dystep = fix(dN/k_dim_y);
        for ix = 1:k_dim_x
            for iy = 1:k_dim_y
                figure(23);
                h = subplot(k_dim_x,k_dim_y,ix*k_dim_y+iy-k_dim_y);
                if k_dim_x < 10
                    set(h,'position',[0.11*iy 1-0.11*ix 0.1 0.1]);
                else
                    set(h,'position',[0.06*(iy-0.5) 1-0.06*(ix+0.1) 0.05 0.05]);
                end

                mx = zeros(k_dim_x,k_dim_y,EP_angle_N);
                my = zeros(k_dim_x,k_dim_y,EP_angle_N);
                mx_error_max = zeros(k_dim_x,k_dim_y,EP_angle_N);
                my_error_max = zeros(k_dim_x,k_dim_y,EP_angle_N);
                m_error_temp = zeros(1,4);
                ix_data = k_dim_x + 1 - ix; % k-space revise
                for it = 1:EP_angle_N
                    it_real = it; 
                    center_loc_x = fix(dM/2)+( ix_data-ceil(k_dim_x/2))*dxstep;
                    center_loc_y = fix(dN/2)+(iy-ceil(k_dim_y/2))*dystep;
                    mx( ix_data,iy,it_real) = GXK(center_loc_x,center_loc_y,it_real);
                    my( ix_data,iy,it_real) = GYK(center_loc_x,center_loc_y,it_real);
                    n = (0:EP_sample_N/4:EP_sample_N);
%                     n = EP_sample_N;
                    plot(n*mx( ix_data,iy,it_real),n*my( ix_data,iy,it_real),'.r',MarkerSize=3);

                    k_bound = rs/(2*w)/7; 
                    axis([-k_bound k_bound -k_bound k_bound]); pbaspect([1 1 1]);
                    set(gca,'ytick',[]);
                    set(gca,'xtick',[]);
                    hold on;
                    % peripheral of local Gx
                    m_error_temp(1) = abs(GXK(center_loc_x,center_loc_y,it_real)-GXK(center_loc_x+ceil(dxstep/2)-1,center_loc_y+ceil(dystep/2)-1,it_real));
                    m_error_temp(2) = abs(GXK(center_loc_x,center_loc_y,it_real)-GXK(center_loc_x-ceil(dxstep/2)+1,center_loc_y+ceil(dystep/2)-1,it_real));
                    m_error_temp(3) = abs(GXK(center_loc_x,center_loc_y,it_real)-GXK(center_loc_x+ceil(dxstep/2)-1,center_loc_y-ceil(dystep/2)+1,it_real));
                    m_error_temp(4) = abs(GXK(center_loc_x,center_loc_y,it_real)-GXK(center_loc_x-ceil(dxstep/2)+1,center_loc_y-ceil(dystep/2)+1,it_real));
                    mx_error_max( ix_data,iy,it_real) = max(m_error_temp);
                    % peripheral of local Gx
                    m_error_temp(1) = abs(GYK(center_loc_y,center_loc_y,it_real)-GYK(center_loc_y+ceil(dystep/2)-1,center_loc_y+ceil(dystep/2)-1,it_real));
                    m_error_temp(2) = abs(GYK(center_loc_y,center_loc_y,it_real)-GYK(center_loc_y-ceil(dystep/2)+1,center_loc_y+ceil(dystep/2)-1,it_real));
                    m_error_temp(3) = abs(GYK(center_loc_y,center_loc_y,it_real)-GYK(center_loc_y+ceil(dystep/2)-1,center_loc_y-ceil(dystep/2)+1,it_real));
                    m_error_temp(4) = abs(GYK(center_loc_y,center_loc_y,it_real)-GYK(center_loc_y-ceil(dystep/2)+1,center_loc_y-ceil(dystep/2)+1,it));
                    my_error_max( ix_data,iy,it_real) = max(m_error_temp);
        
                    G_mag = sqrt(mx( ix_data,iy,it_real)^2+my( ix_data,iy,it_real)^2);                       
                    G_subplot_angle( ix_data,iy,it_real) = G_mag;
                    Gx_subplot_angle( ix_data,iy,it_real) = mx( ix_data,iy,it_real);
                    Gy_subplot_angle( ix_data,iy,it_real) = my( ix_data,iy,it_real);

                    Gx(it_real) = Gx_subplot_angle(ix_data,iy,it_real);
                    Gy(it_real) = Gy_subplot_angle(ix_data,iy,it_real);
                    G_angle = angle(Gx + Gy*1i)/(2*pi)*360;
                    G_angle = wrapTo360(G_angle - G_angle(1)); % Initial is 360 or 0
                    G_abs(it_real) = G_subplot_angle(ix_data,iy,it_real);
                end

                figure(101);
                subplot(k_dim_x,k_dim_y,ix*k_dim_y+iy-k_dim_y);
                plot(index_angle,Gx*1e3,".");
                axis([min(index_angle) max(index_angle) min(Gx_subplot_angle,[],'all')*1e3 max(Gx_subplot_angle,[],'all')*1e3])
                sgtitle('|G_x|')

                figure(102);
                subplot(k_dim_x,k_dim_y,ix*k_dim_y+iy-k_dim_y);
                plot(index_angle,Gy*1e3,".");
                axis([min(index_angle) max(index_angle) min(Gy_subplot_angle,[],'all')*1e3 max(Gy_subplot_angle,[],'all')*1e3])
                sgtitle('|G_y|')

                figure(100);
                subplot(k_dim_x,k_dim_y,ix*k_dim_y+iy-k_dim_y);
                yyaxis left
                plot(index_angle,G_angle,".");
                axis([min(index_angle) max(index_angle) 0 360])
                yyaxis right
                plot(index_angle,G_abs*1e3,".");
                axis([min(index_angle) max(index_angle) 0 max(G_subplot_angle,[],'all')*1e3])
                sgtitle('G phase and magnitude')
                k_bound_real = abs(k_bound./cos(G_angle/360*2*pi));
%                 Eval(xc_idx,k_dim_y*(ix-1)+iy) = Eval_k(G_abs,G_angle,thre_k,EP_angle_N);
                [Eval(xc_idx,k_dim_y*(ix-1)+iy),Eval_all(:,xc_idx,k_dim_y*(ix-1)+iy),Angle_all(:,k_dim_y*(ix-1)+iy,:)] = Eval_k2(G_abs*EP_sample_N,G_angle,thre_k,k_bound_real);
%                 title(num2str(Eval(xc_idx,k_dim_y*(ix-1)+iy)))
                G_abs_all(xc_idx,k_dim_y*(ix-1)+iy) = sum(G_abs);
                
                figure(110);
                subplot(k_dim_x,k_dim_y,ix*k_dim_y+iy-k_dim_y);
                plot(index_angle,Eval_all(:,xc_idx,k_dim_y*(ix-1)+iy),".");
                axis([min(index_angle) max(index_angle) 0 max(max(Eval_all))])
                sgtitle('Eval_all')



            end
        end
%         title(num2str(it))
%     saveas(h2,strcat('plot_local_k\',version,'_local_k.jpg'));
end

Eval_angle = sum(Eval_all,[2,3]);
figure(111)
plot(Eval_angle)
[~, sortedIndex] = sort(Eval_angle, 'descend');
sortedIndex

end
Angle_temp(:,:) = Angle_all(1,:,:);