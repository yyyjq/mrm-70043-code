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