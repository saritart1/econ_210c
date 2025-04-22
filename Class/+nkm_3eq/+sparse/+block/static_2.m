function [y, T] = static_2(y, x, params, sparse_rowval, sparse_colval, sparse_colptr, T)
  y(4)=(-(params(5)*params(1)*params(2)*(params(6)-1)*y(5)))/(params(1)+(1-params(5))*params(2));
end
