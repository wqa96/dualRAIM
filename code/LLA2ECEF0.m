function [X, Y, Z] = LLA2ECEF0(long, lat, alt)
% 经纬高(long, lat, alt)转ECEF坐标(X, Y, Z)
    Re = 6378137.0;             %椭球长半轴
%     f = 1/298.257223563;        %扁率
    e = 0.081819190842622;      %偏心率
    e2 = e * e;
    d2r = pi/180;
    
    L = long * d2r;
	B = lat * d2r;
	H = alt;

	RN = Re ./ sqrt(1 - e2 * sin(B) .* sin(B));
    RNH = RN + H;
	X = RNH .* cos(B) .* cos(L);
	Y = RNH .* cos(B) .* sin(L);
	Z = (RN * (1 - e2) + H) .* sin(B);
    
%     xyz = lla2ecef([lat, long, alt], 'WGS84');
%     X = xyz(:, 1);
%     Y = xyz(:, 2);
%     Z = xyz(:, 3);
end