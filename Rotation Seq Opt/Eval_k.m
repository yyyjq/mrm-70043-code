%%
% Calculate the outline area of k space

% Input: Magnitude of gradient (G_abs); Phase of gradient (G_angle); The k
% bound (k_bound); Angle interval (thre); 

% Output: Eval: The value of evaluation; Angle: The truncated angle
% distribution; Eval_all: The evaluation value for every spokes; 

function [Eval,Eval_all,Angle] = Eval_k2(G_abs,G_angle,thre,k_bound)

thre = thre/2;

[~,Index] = sort(G_abs,'ascend'); % sort the data as ascend
Effect_area = zeros(1,length(G_abs));
Angle_detail = zeros(2,length(G_abs));

for i = 1:length(G_angle)
    Index_current = Index(i);
    Angle_min = G_angle(Index_current)-thre;
    Angle_max = G_angle(Index_current)+thre;
%     if Angle_min <= 0
%         Angle_min = Angle_min+360;
%     else if Angle_min >= 360
%             Angle_min = Angle_min - 360;
%     end
%     end
%     if Angle_max <= 0
%         Angle_max = Angle_max+360;
%     else if Angle_max >= 360
%         Angle_max = Angle_max - 360;
%     end
%     end

    if i == length(G_angle)
    else
    for j = i + 1:length(G_angle)
        Index_temp = Index(j);
        if (Angle_min > G_angle(Index_temp)+thre)||(G_angle(Index_temp)-thre > Angle_max)    % three situations
        elseif (Angle_min < G_angle(Index_temp)+thre)&&(Angle_max > G_angle(Index_temp)+thre)
                Angle_min = G_angle(Index_temp) + thre;
        elseif (Angle_max > G_angle(Index_temp)-thre)&&(Angle_min < G_angle(Index_temp)-thre)
                Angle_max = G_angle(Index_temp) - thre;
        else
                Angle_max = Angle_min;
        end
    end
    end
            Effect_area(1,Index_current) = (Angle_max - Angle_min)/360*2*pi;
            Angle_detail(1,Index_current) = Angle_max;
            Angle_detail(2,Index_current) = Angle_min;
end
for i = 1:length(G_angle)
    if G_abs(i) > k_bound
        G_abs(i) = k_bound;
    end
end
Angle = Angle_detail;
Eval_all = (1/2*Effect_area.*G_abs.^2)./(2*thre./360*2*pi/2*(k_bound).^2);
% Eval_all = (Effect_area)./(thre/360*2*pi);
% Eval_all = (Effect_area.*G_abs)./(2*thre./360*2*pi*(k_bound));
Eval = sum(Eval_all,'all')/length(G_abs);
end

%%
% FOV is dependent to the angle, because of the rectangle shape of image