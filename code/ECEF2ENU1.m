function [E, N, U] = ECEF2ENU1(Xe, Ye, Ze, long0, lat0)
% ECEFϵתENUϵ
    d2r = pi / 180;
    long = long0 * d2r;
    lat = lat0 * d2r;

%     Cen = [-sin(long), cos(long), 0;
%         -sin(lat)*cos(long), -sin(lat)*sin(long), cos(lat);
%         cos(lat)*cos(long), cos(lat)*sin(long), sin(lat)];
    E = -Xe.*sin(long) + Ye.*cos(long);
    N = -Xe.*sin(lat).*cos(long) - Ye.*sin(lat).*sin(long) + Ze.*cos(lat);
    U = Xe.*cos(lat).*cos(long) + Ye.*cos(lat).*sin(long) + Ze.*sin(lat);
end
