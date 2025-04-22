function [y, T, residual, g1] = dynamic_2(y, x, params, steady_state, sparse_rowval, sparse_colval, sparse_colptr, T)
residual=NaN(3, 1);
  residual(1)=(y(7))-(y(12)-(1-params(5))/params(1)*(y(8)-y(11)-y(9)));
  residual(2)=(y(8))-(y(6)*params(7)+y(7)*params(10));
  residual(3)=(y(6))-(params(3)*params(9)*y(7)+params(4)*y(11));
if nargout > 3
    g1_v = NaN(7, 1);
g1_v(1)=(1-params(5))/params(1);
g1_v(2)=1;
g1_v(3)=1;
g1_v(4)=(-params(10));
g1_v(5)=(-(params(3)*params(9)));
g1_v(6)=(-params(7));
g1_v(7)=1;
    if ~isoctave && matlab_ver_less_than('9.8')
        sparse_rowval = double(sparse_rowval);
        sparse_colval = double(sparse_colval);
    end
    g1 = sparse(sparse_rowval, sparse_colval, g1_v, 3, 3);
end
end
