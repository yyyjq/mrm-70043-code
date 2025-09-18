% Initialization script for image recon
% Parameter initialization & variable preparation for phase optimization
%% 1. field/FOV/sample/phantom setting

% ===== FOV =====
xc = 30; yc = 20;               % center location: 
rs = 128;                        % resolution
l = 29; w = 29;                 % FOV length (y)/width (x)

step = l/rs; %Delta w

% ===== B0 map =====
B0_map = 1;                     % 1 - Halbach, 2 - IO-Ring
truncate_B0 = 1;

% ===== B1 =====
Coil_N = 1;
Coilsurf = 0; % 1: use Bio-savart law  0: uniform sensivity map
%when Coilsurf = 1
Coilautoposition = 1; % 1: generate the data by algrithm  0: basing on experiment data
%wehn Coilautoposition = 0
Coil_r = 4.2e-3; %m
Coil_p(:,1)= [0;0;0;0;0;0;0;0]*1e-3; % x m
Coil_p(:,2)= [0;0;0;0;0;0;0;0]*1e-3; % y m
Coil_p(:,3)= [0;0;0;0;0;0;0;0]*1e-3; % z m
Coil_Azi = [0 0 0 0 0 0 0 0]; % coil angle
Coil_currentDir = [1 1 1 1 1 1 1 1];% 

% ===== time convention =====
time_con = 2;                   % 1: exp(-i omega t); 2: exp(i omega t)

% ===== Mea DATA=====
dataName = '70mm_copper_wobox\';

% ===== sample parameters =====
EP_sample_N = 260;
EP_angle_N = 144;
EP_angle_Total = (2*pi + 3*pi/180 - pi/180*2.5); % middle 3
EP_dt = 0.5e-6; % s
RxDelay = 50*1e-6; % s % without Q-damp: 50
T2 = 0;                         % 0 - no delay; 1 - custom delay mode
T2_coeff = 10;                 % T2 = 0.1*T2_coeff ms
key_T2 = 'T2_0p1ms_NoDelay';    % delay data key (only in T2 = 1 mode)
EP_f_center = 2.84475e6;           % 13 bottles: 2.85, 2.837
EP_angle_start= 0;           % yidan: pi; gj: pi; Tingou:pi/2
% rotate direction of the magnet array
a = 1;                          %a: 1 clockwise -1 counterclockwise
EP_ratio = 1;


sigS = -1;

% ===== plot signal =====
plot_k = 0;
k_dim = 7;
plotSig = 0;            % Whether to plot the signal after optimization
phase_en = 0;           % Whether to add the phase correction

% ===== reconstruction =====
reCon = 1;              % Whether to do image reconstruction

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

B0_filename = 'B0_measured_ploy55_Junqi_2022_11_2.mat';
load(B0_filename);
B0_meas_Halb = B0_meas_Halb_poly55;

B0_ix = reshape(B0_meas_Halb(:,:,1), numel(B0_meas_Halb(:,:,1)), 1);
B0_iy = reshape(B0_meas_Halb(:,:,2), numel(B0_meas_Halb(:,:,2)), 1);
B0_im = reshape(B0_meas_Halb(:,:,3).*1e-3, numel(B0_meas_Halb(:,:,3)), 1);

B0fovi = griddata(B0_ix, B0_iy, B0_im, FOV_x, FOV_y);
B0FOV = reshape(B0fovi, size(FOV_X,1), size(FOV_X,2));

B0fovMax = max(B0_im(~isnan(B0_im)));
B0fovMin = min(B0_im(~isnan(B0_im)));

B0fov = reshape(B0FOV, numel(FOV_X), 1);

%% === Coil senstivity ===
% choose the automatical coil position generation
if Coilsurf == 1
    if Coilautoposition == 1
        Coil_r = abs(w+1i*l)./2*sin(pi/Coil_N)*1e-3; %m
        r=abs(w+1i*l)./2; %mm
        deltaangle=2*pi./Coil_N;
        for i=1:Coil_N
            Angle = deltaangle*(i-1)-pi./2+2*pi./360*(2*rand(1)-1);
            Coil_p( i, 1) = r*cos(Angle)+xc;
            Coil_p( i, 2) = r*sin(Angle)+yc;
            Coil_p( i, 3) = 0;
            Coil_Azi( i ) = deltaangle*(i-1)-pi./2;
        end
        Coil_p = Coil_p*1e-3;  % normalize to m
        Coil_currentDir = ones(Coil_N,1);
    end
    B1x = zeros(numel(FOV_X),Coil_N); 
    B1y = B1x; 
    B1z = B1x;
    for n = 1:Coil_N
        [B1x(:,n), B1y(:,n), B1z(:,n), ~] = ...
            funOthercoilBiotSavart(FOV_X,FOV_Y,Coil_p(n,:),Coil_Azi(n),...
            0,Coil_r,36,Coil_currentDir(n),1,0,1);
    end
else
    % Calculated sensitivity map
    B1_filename = 'Bz_D60_coil.mat';
    load(B1_filename);
    
    B1_ix = reshape(x, numel(x), 1)+xc;
    B1_iy = reshape(y, numel(y), 1)+yc;
    B1_im = reshape(Bz.*1e-3, numel(Bz), 1);
%     
    B1z = griddata(B1_ix, B1_iy, B1_im, FOV_x, FOV_y);
    B1z = B1z./max(B1z,[],'all');
    B1x = zeros(numel(FOV_X),Coil_N);
    B1y = B1x;

end

%% ===== Sample parameters =====
if T2 == 0
    key_T2 = 'NoT2_NoDelay';
    delayCorrArray = zeros(1,EP_angle_N);
    T2star = zeros(EP_angle_N,1);
    selectDelay = 0;
elseif T2 == 1
    temp = load(strcat('correction_',key_T2,'.txt'));
    T2star = temp(:,1).*1e-4*T2_coeff;
    delayCorrArray = temp(:,2).*1e-4;
    selectDelay = 1;
end

E_fshift = 0;
E_xcm = -2; % eccentric error
E_ycm = 0; % eccentric error 
E_tra = zeros(EP_angle_N,2);

%% === load Signal & Sparse the Sampling ===
A_s_max = zeros(1,EP_angle_N);

for i_an = 1:1:EP_angle_N

        % load measured data
        load(strcat(dataName,num2str(i_an-1),'\data.csv')); % us, uV
        A_s = data(:,2)+sigS.*1i.*data(:,3);
        t_exp = data(:,1);
        A_s_max(i_an) = max(abs(A_s));
        if i_an == 1
            SignalexpA = zeros(numel(A_s),EP_angle_N);
        end
        SignalexpA(:,i_an)=A_s;
        
end

start_threshold = 0;

Signalexpi = SignalexpA(1+start_threshold:EP_sample_N+start_threshold,:);
clear Signalexp;
Signalexpt = reshape(Signalexpi,EP_sample_N,Coil_N,EP_angle_N); % size: 128 x 1 x 144, 144 angles

%% ===== Build EB0_long as the base =====
const_gamma = 42580000;

fshift = E_fshift;
xcm = E_xcm;
ycm = E_ycm;

idx_valid = ~isnan(B0_im);
B0fit = fit([B0_ix(idx_valid), B0_iy(idx_valid)],B0_im(idx_valid),'poly55');


addOnPoint = round(10*RxDelay/EP_dt);
sampleLong = EP_sample_N+addOnPoint;
EB0_long_all = zeros(numel(FOV_x),sampleLong,EP_angle_N);
Gbx = zeros(size(FOV_X,1),size(FOV_X,2),EP_angle_N);
Gby = zeros(size(FOV_X,1),size(FOV_X,2),EP_angle_N);
G_indicator = zeros(1,EP_angle_N);
% Inhomogeneity matrix
boundary_f_max_min = zeros(2,EP_angle_N);
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
    
    % Inhomogeneity
    boundary_f_max_min(1,i_an) = max(f-EP_f_center);
    boundary_f_max_min(2,i_an) = min(f-EP_f_center);
    
    EB0M = reshape(EB0_org,size(FOV_X,1),size(FOV_X,2));
    [Gbx(:,:,i_an),Gby(:,:,i_an)] = gradient(EB0M,step); 
end

