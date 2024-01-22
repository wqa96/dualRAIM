function [pos, r, perr, G, H, cnt] = pvtpos0(pr, sat, pos0, th, maxcnt)
%迭代最小二乘法求解位置
    c = 299792458;              % 光速（m/s）
    cnt = 0;
    pos = pos0;
    while (1)
        cnt = cnt + 1;
        earth_rot = 7.292115e-5 * (sat(:, 1) * pos(2) - sat(:, 2) * pos(1)) / c;
        r = sqrt(sum((sat - pos(1:3)).^2, 2));
        G = [-(sat - pos(1:3))./r, ...
            ones(size(sat, 1), 1)];
        H = inv(G'*G);
        perr = pr - r - pos(4) - earth_rot;
        if rcond(H) < 1e-15
            warning('矩阵奇异');
            break;
        end
        dp = H * G'* perr;
        pos = pos + dp';
        if(sum(abs(dp)) < th) || (cnt >= maxcnt)
            break;
        end
    end
end