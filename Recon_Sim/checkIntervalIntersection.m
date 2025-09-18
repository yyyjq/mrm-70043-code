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
