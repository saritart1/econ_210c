function [residual, T_order, T] = dynamic_resid(y, x, params, steady_state, T_order, T)
if nargin < 6
    T_order = -1;
    T = NaN(0, 1);
end
[T_order, T] = nkm_3eq.sparse.dynamic_resid_tt(y, x, params, steady_state, T_order, T);
residual = NaN(5, 1);
    residual(1) = (y(6)) - (params(3)*params(9)*y(7)+params(4)*y(11));
    residual(2) = (y(7)) - (y(12)-(1-params(5))/params(1)*(y(8)-y(11)-y(9)));
    residual(3) = (y(9)) - ((-(params(5)*params(1)*params(2)*(params(6)-1)*y(10)))/(params(1)+(1-params(5))*params(2)));
    residual(4) = (y(8)) - (y(6)*params(7)+y(7)*params(10));
    residual(5) = (y(10)) - (params(6)*y(5)+x(1));
end
