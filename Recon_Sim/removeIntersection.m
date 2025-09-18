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
