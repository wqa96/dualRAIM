function [r] = mynorm(m, p, dim)
% m: ���ݾ���
% p: ����
% dim: ά��
    if (p == Inf)
        r = max(abs(m), [], dim);
    elseif (p == -Inf)
        r = min(abs(m), [], dim);
    elseif (p > 0)
        r = sum(abs(m).^p, dim).^(1/p);
    else
        error('p must be a positive real value or Inf or -Inf.');
    end
end