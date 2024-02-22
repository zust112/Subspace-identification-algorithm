%%%%%%%%%%�ο����ף�A new subspace identification approach based on principal component analysis
function [A, B, C, D,ss]=SIMPCA(y,u,f,p,n)
%   y���������
%   u����������
%   f��δ��Hankel������ά��
%   p�ǹ�ȥHankel������ά��
%   n��ϵͳ����
%   A,B,C,D��״̬�ռ����


%   �ж��������
if nargin < 5
    n = [];
end

%   �ж��������
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
    error('��������������ȱ������');
end

num = num_u;


%   ����Hankel����
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
Zp = [Yp; Up];  %ע���������SIM��ͬ��ʹ���˸���������
Zf = [Yf; Uf];  

%����SVD����PCA
 [P, S, V] = svd((Zf * (Zp)') / num);%%ʽ25


%   �״α�ʶ
if (isempty(n))
    ss=diag(S);
    figure
    bar(ss);
    n=input('�״�Ϊ��');
end
ss=diag(S);

Pr = P(:, (m * f)+n+1:end); %��ȡ��������غɾ���
Pry = Pr(1:l * f, :);%%��չ�ɹ۾����������
Pru = Pr(l * f + 1:end, :);
% gammaf_per = Pry'; %ʽ28�У�Mȡ��λ��
gammaf = orthcomp(Pry); %����������
%  gammaf = null(gammaf_per); %��������������һ�ַ�������ֵ�����ͬ������A,B,C,D��Ҳ��ͬ����ϵͳ���Ի���һ��

%   ����A��C
C = gammaf(1:l, :);
A = pinv(gammaf(1:l * (f - 1), :)) * gammaf(l + 1:l * f, :);

%   ׼������B��D���ο�ʽ32-36
phi = -Pry';%ʽ32
Phi = Pru';%ʽ33
Lhs = zeros((l * f - n) * f, l * f);
% Rhs = zeros((l * m - n) * i, m);
for k = 1 : f
    Lhs((k - 1) * (l * f - n) + 1:k * (l * f - n), 1:(f - k + 1) * l) = phi(:, (k - 1) * l + 1:f * l);
    Rhs((k - 1) * (l * f - n) + 1:k * (l * f - n), :) = Phi(:,(k - 1) * m + 1:m * k);
end

%   �����С����
Hf = Lhs \ Rhs;
sol = pinv([eye(l) zeros(l, n); zeros(l * (f - 1), l) gammaf(1:l * (f - 1), :)]) * Hf; %ʽ��40��

%   �õ�B��D����
D = sol(1:l, :);
B = sol(l + 1:l+n, :);

end