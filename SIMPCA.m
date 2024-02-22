%%%%%%%%%%参考文献：A new subspace identification approach based on principal component analysis
function [A, B, C, D,ss]=SIMPCA(y,u,f,p,n)
%   y是输出数据
%   u是输入数据
%   f是未来Hankel矩阵行维度
%   p是过去Hankel矩阵行维度
%   n是系统阶数
%   A,B,C,D是状态空间矩阵


%   判断输入参数
if nargin < 5
    n = [];
end

%   判断输入参数
[m, num_u] = size(u);
if(m > num_u)
    u = u';
    [m, num_u] = size(u);
end

[l, num_y] = size(y);
if(l > num_y)
    y = y';
    [l, num_y] = size(y);
end

if(num_u ~= num_y)
    error('输入输出向量长度必须相等');
end

num = num_u;


%   计算Hankel矩阵
for i=1:num-f-p+1
    for j=1:p+f
        Y((j-1)*l+1:j*l,i)=y(:,j+i-1);
    end
end

for i=1:num-f-p+1
    for j=1:p+f
        U((j-1)*m+1:j*m,i)=u(:,j+i-1);
    end
end

Yp = Y(1:p*l, :);
Yf = Y(p*l + 1:(p+f) * l, :);
Up = U(1:p*m, :);
Uf = U(p*m + 1:(p+f) * m, :);
Zp = [Yp; Up];  %注意与基本的SIM不同（使用了辅助变量）
Zf = [Yf; Uf];  

%利用SVD进行PCA
 [P, S, V] = svd((Zf * (Zp)') / num);%%式25


%   阶次辨识
if (isempty(n))
    ss=diag(S);
    figure
    bar(ss);
    n=input('阶次为？');
end
ss=diag(S);

Pr = P(:, (m * f)+n+1:end); %提取残余项的载荷矩阵
Pry = Pr(1:l * f, :);%%扩展可观矩阵的正交补
Pru = Pr(l * f + 1:end, :);
% gammaf_per = Pry'; %式28中，M取单位阵
gammaf = orthcomp(Pry); %计算正交补
%  gammaf = null(gammaf_per); %计算正交补的另一种方法，数值结果不同，最终A,B,C,D阵也不同，但系统特性基本一样

%   计算A和C
C = gammaf(1:l, :);
A = pinv(gammaf(1:l * (f - 1), :)) * gammaf(l + 1:l * f, :);

%   准备计算B和D，参考式32-36
phi = -Pry';%式32
Phi = Pru';%式33
Lhs = zeros((l * f - n) * f, l * f);
% Rhs = zeros((l * m - n) * i, m);
for k = 1 : f
    Lhs((k - 1) * (l * f - n) + 1:k * (l * f - n), 1:(f - k + 1) * l) = phi(:, (k - 1) * l + 1:f * l);
    Rhs((k - 1) * (l * f - n) + 1:k * (l * f - n), :) = Phi(:,(k - 1) * m + 1:m * k);
end

%   求解最小二乘
Hf = Lhs \ Rhs;
sol = pinv([eye(l) zeros(l, n); zeros(l * (f - 1), l) gammaf(1:l * (f - 1), :)]) * Hf; %式（40）

%   得到B和D矩阵
D = sol(1:l, :);
B = sol(l + 1:l+n, :);

end