function [X, Y, Z] = LLA2ECEF0(long, lat, alt)
% ��γ��(long, lat, alt)תECEF����(X, Y, Z)
    Re = 6378137.0;             %���򳤰���
%     f = 1/298.257223563;        %����
    e = 0.081819190842622;      %ƫ����
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