function [y, T] = dynamic_1(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(10)=params(6)*y(5)+x(1);
  y(9)=(-(params(5)*params(1)*params(2)*(params(6)-1)*y(10)))/(params(1)+(1-params(5))*params(2));
end
