%%
% Calculate the outline area of k space

% Input: Magnitude of gradient (G_abs); Phase of gradient (G_angle); The k
% bound (k_bound); Angle interval (thre); 

% Output: Eval: The value of evaluation; Angle: The truncated angle
% distribution; Eval_all: The evaluation value for every spokes; 

function [Eval,Eval_all,Angle] = Eval_k3(G_abs,G_angle,thre,k_bound)

thre = thre/2;

[~,Index] = sort(G_abs,'ascend'); % sort the data as ascend
Effect_area = zeros(1,length(G_abs));
Angle_detail = zeros(4,length(G_abs));

for i = 1:length(G_angle)
    Index_current = Index(i);
    Angle_min = G_angle(Index_current)-thre;
    Angle_max = G_angle(Index_current)+thre;

    Angle_test = CheckAngleInterval(Angle_min,Angle_max); 

    if i == length(G_angle)
    else
    for j = i + 1:length(G_angle)
        Index_temp = Index(j);
        Angle_temp = CheckAngleInterval(G_angle(Index_temp)-thre, G_angle(Index_temp)+thre);
        [hasIntersection, intersection] = checkIntervalIntersection(Angle_test, Angle_temp);
        Angle_inter = Angle_test;
        if hasIntersection(1,1) || hasIntersection(1,2)
            Angle_inter(1,:) = removeIntersection(Angle_test(1,:),intersection(1,:));
            if hasIntersection(2,1) || hasIntersection(2,2)
                Angle_inter(2,:) = removeIntersection(Angle_test(2,:),intersection(2,:));
            end
        else
            Angle_inter(2,:) = removeIntersection(Angle_test(2,:),intersection(1,:));
        end
        Angle_test = Angle_inter;
    end
    end
            Effect_area(1,Index_current) = (sum(Angle_test(:,2) - Angle_test(:,1)))/360*2*pi;
            Angle_detail(1:2,Index_current) = Angle_test(1,:);
            Angle_detail(3:4,Index_current) = Angle_test(2,:);
            
end
for i = 1:length(G_angle)
    if G_abs(i) > k_bound
        G_abs(i) = k_bound;
    end
end

%% Angle detail
Angle = Angle_detail;

%% Area calculation
Eval_all = (1/2*Effect_area.*G_abs.^2)./(2*thre./360*2*pi/2*(k_bound).^2);
% Eval_all = (Effect_area)./(thre/360*2*pi);
% Eval_all = (Effect_area.*G_abs)./(2*thre./360*2*pi*(k_bound));
Eval = sum(Eval_all,'all')/length(G_abs);
end

%%
% FOV is dependent to the angle, because of the rectangle shape of image

%%  CheckAngleInterval
function AngleInter = CheckAngleInterval(Angle_min, Angle_max)

AngleInter = [Angle_min, Angle_max;
    0,0];

if Angle_min <= 0
    AngleInter = [0, Angle_max;
        Angle_min + 360, 360];
end
if Angle_max >= 360
    AngleInter = [0, Angle_max - 360;
        Angle_min, 360];
end
end

%% RemoveIntersection
function remaining = removeIntersection(interval1, interval2)
remaining = [0,0];
  %% interval1 and interval2 are two intervals represented as [start, end]
    % Ensure the input intervals are valid, i.e., start <= end
    if interval1(1) > interval1(2) || interval2(1) > interval2(2)
        error('Invalid intervals, make sure start is less than or equal to end.');
    end

    % Calculate the intersection part
    start = max(interval1(1), interval2(1));
    finish = min(interval1(2), interval2(2));

    % Check if there's an intersection
    if start > finish
        % No intersection, intervals 1 and 2 remain the same
        remaining1 = interval1;
        remaining2 = interval2;
    else
        % There is an intersection, remove the intersection part
        if interval1(1) < interval2(1)
            remaining1 = [interval1(1), start];
            remaining2 = [finish, interval2(2)];
        else
            remaining2 = [interval2(1), start];
            remaining1 = [finish, interval1(2)];
        end
    end
    if remaining1(1) ~= remaining1(2)
            remaining = remaining1;
    else if remaining2(1) ~= remaining2(2)
            remaining = remaining2;
    else
    end
    end
end

%% checkIntervalIntersection
function [hasIntersection, intersection] = checkIntervalIntersection(Input1, Input2)
    % interval1 and interval2 represent two intervals, each in the form [start, end]
    hasIntersection = zeros(2,2);
    k = 1;
    intersection = zeros(2,2);
    for i = 1:2
        for j = 1:2
            interval1 = Input1(i,:); interval2 = Input2(j,:); 
            % Ensure that the input intervals are valid, i.e., start <= end
            if interval1(1) >= interval1(2) || interval2(1) >= interval2(2)
                hasIntersection(i,j) = false; 
            else    
            % Check for intersection
            if interval1(2) < interval2(1) || interval2(2) < interval1(1)
                hasIntersection(i,j) = false;
            else
                hasIntersection(i,j) = true;
        
                % Calculate the intersection part
                start = max(interval1(1), interval2(1));
                finish = min(interval1(2), interval2(2));
                if start ~= finish     
                    intersection(k,:) = [start, finish];
                    k = k+1;
                end
            end
            end
        end
    end
end
