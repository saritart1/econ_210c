function [residual, T_order, T] = static_resid(y, x, params, T_order, T)
if nargin < 5
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = nkm_3eq.sparse.static_resid_tt(y, x, params, T_order, T);
residual = NaN(5, 1);
    residual(1) = (y(1)) - (params(3)*params(9)*y(2)+y(1)*params(4));
    residual(2) = (y(2)) - (y(2)-(1-params(5))/params(1)*(y(3)-y(1)-y(4)));
    residual(3) = (y(4)) - ((-(params(5)*params(1)*params(2)*(params(6)-1)*y(5)))/(params(1)+(1-params(5))*params(2)));
    residual(4) = (y(3)) - (y(1)*params(7)+y(2)*params(10));
    residual(5) = (y(5)) - (params(6)*y(5)+x(1));
end
