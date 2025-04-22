function [lhs, rhs] = static_resid(y, x, params)
T = NaN(0, 1);
lhs = NaN(5, 1);
rhs = NaN(5, 1);
lhs(1) = y(1);
rhs(1) = params(3)*params(9)*y(2)+y(1)*params(4);
lhs(2) = y(2);
rhs(2) = y(2)-(1-params(5))/params(1)*(y(3)-y(1)-y(4));
lhs(3) = y(4);
rhs(3) = (-(params(5)*params(1)*params(2)*(params(6)-1)*y(5)))/(params(1)+(1-params(5))*params(2));
lhs(4) = y(3);
rhs(4) = y(1)*params(7)+y(2)*params(10);
lhs(5) = y(5);
rhs(5) = params(6)*y(5)+x(1);
end
