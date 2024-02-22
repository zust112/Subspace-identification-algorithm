function N = orthcomp(A)
%�������������
%   A�Ǵ������
%   N�Ǿ���A��������

A=A';
[R, pivcol] = rref(A, sqrt(eps)); %��A���ɽ���״,epsΪrank�����и����Ĺ���
[m, n] = size(A);
r = length(pivcol);
freecol = 1:n;
freecol(pivcol) = [];
N = zeros(n, n-r);
N(freecol, : ) = eye(n-r);
N(pivcol,  : ) = -R(1:r, freecol);