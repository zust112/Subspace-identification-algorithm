%%%%%%%%%%计算矩阵的正交补
function Ac=orthc(A)
[mm,nn]=size(A);
if (mm < nn)
    A=A';
    [mm,nn]=size(A);
end
[QQ,~]=qr(A);
Ac=QQ(:,nn+1:mm);
end