function [A,B,C,D]=SIMPCA_CN(y,u,f,p,n)
%*******************************************
%   入参出参说明
%   y是输出数据序列
%   u是输入数据序列
%   f是未来Hankel矩阵行维度
%   p是过去Hankel矩阵行维度
%   n是系统阶数（可为空）
%   A,B,C,D是辨识得到的状态空间矩阵
%*******************************************

%   判断输入参数并转为行向量
if nargin < 5
    n = [];
end
%   nu,ny为输入输出变量数、sm为样本数
[ny,sm]=size(y);
if sm<ny
    y=y';
    [ny,~]=size(y);
end
[nu,sm]=size(u);
if sm<nu
    u=u';
    [nu,~]=size(u);
end 

%********************************************
%   构建Hankel矩阵
U=Hankel(p,f,u);
Y=Hankel(p,f,y);

%********************************************
%   数据矩阵划分
Yp=Y(1:p*ny,:);
Yf=Y(p*ny+1:(p+f)*ny,:);
Up=U(1:p*nu,:);
Uf=U(p*nu+1:(p+f)*nu,:);

%********************************************
%   噪声估计(式19-20)
ZQR = [Yp;Up;Uf;Yf];
[Q,L]= qr(ZQR',0);%%将大矩阵转置使其可用QR中的经济分解压缩矩阵大小，提高计算速度
Q=Q';
L=L';
L44 = L((nu+ny)*p+nu*f+1:end,(nu+ny)*p+nu*f+1:end);
Q4 = Q((nu+ny)*p+nu*f+1:end,:);
Ef=L44*Q4;
Yfb=Yf-Ef;

%********************************************
%   PCA（数值上使用svd计算）
[P,S,~]=svd([Yfb;Uf]);%%式9

%********************************************
%   阶次辨识
if (isempty(n))
    ss=diag(S);
    figure
    bar(ss);
    n=input('阶次为？');
end
% ss=diag(S);

%********************************************
%   提取残余项载荷矩阵
Pr=P(:,(nu*f)+n+1:end);%%式10

%********************************************
%   继续拆分(式11)
Pry=Pr(1:ny*f,:);
Pru=Pr(ny*f+1:end,:);

% Pre=Pr(ny*(f+p)+1:end,:);

%********************************************
%******************计算正交补**************%
%   方法一
% Fi=null(Pry');%%Pry'必须为满秩矩阵，否则运行时可能有bug

%   方法二
% Fi = orthcomp(Pry); %计算正交补

%   方法三
[mm,nn]=size(Pry);
if (mm < nn)
    Pry=Pry';
    [mm,nn]=size(Pry);
end
[QQ,~]=qr(Pry);
Fi=QQ(:,nn+1:mm);

Fi=Fi(:,1:n);

%********************************************
%   提取A和C (式21)
C=Fi(1:ny,:);
Fi1=Fi(1:ny*(f-1),:); %%去掉后ny行
A=pinv(Fi1)*Fi(ny+1:ny*f,:);

%********************************************
%   原方法提取B和D
% Fai=-(Pry');
% Psai=Pru';
% 
% FFai=zeros((ny*f-n)*f,(ny*f));
% for i=1:f
%   FFai((ny*f-n)*(i-1)+1:(ny*f-n)*i,1:(f-i+1)*ny)=Fai(:,(i-1)*ny+1:f*ny);
%   PPsai((i-1)*(ny*f-n)+1:(ny*f-n)*i,:)=Psai(:,(i-1)*nu+1:i*nu);
% end
% 
% Hf1=FFai\PPsai;
% DB=pinv([eye(ny)  zeros(ny,n); zeros(ny*(f-1),ny) Fi(1:ny*(f-1),:)])*Hf1;
% D=DB(1:ny,:);
% B=DB(ny+1:ny+n,:);
%********************************************
%   改进方法

%   截短为方阵
L=L(:,1:(f+p)*(ny+nu));

%   从A和C中重新计算Fi和Fi1 (式22)
Fi=C;
for k=2:f
   Fi((k-1)*ny+1:k*ny,:)=Fi((k-2)*ny+1:(k-1)*ny,:)*A;  
end
Fi1=Fi(1:ny*(f-1),:);%%去掉后L行的扩展可观
Fi_inv=pinv(Fi);%%扩展可观伪逆
Fi1_inv=pinv(Fi1);%%去掉后L行的扩展可观伪逆

%   用新的Fi和Fi1重新确定线性方程组 (式24)
Zi=L((nu+ny)*p+nu*f+1:end,1:(nu+ny)*p+nu*f);%%Zi=[R41 R42 R43]
Zi1=L((nu+ny)*p+nu*f+1+ny:end,1:(nu+ny)*p+nu*f+ny);%%Z(i+l)
Rhs=[Fi_inv*Zi          zeros(n,ny);
     L((nu+ny)*p+1:(nu+ny)*p+nu*f,1:(nu+ny)*p+nu*f+ny)];
 
Lhs=[Fi1_inv*Zi1;
    L((nu+ny)*p+nu*f+1:(nu+ny)*p+nu*f+ny,1:(nu+ny)*p+nu*f+ny)];

%   计算T和Qf
T=Lhs-[A;C]*Rhs(1:n,:);%%式23
T=T(:,1:(nu+ny)*p);
Qf=L((nu+ny)*p+1:(nu+ny)*p+nu*f,1:(nu+ny)*p);%%未来输入

%   建立包含B和D线性方程组
L1=A*Fi_inv;
L2=C*Fi_inv;
M=[zeros(n,ny),Fi1_inv];
X=[eye(ny),zeros(ny,n);zeros(ny*(f-1),ny),Fi1];%%式26

%   计算N和kron内积
tot=0;
for k=1:f
   N=[M(:,(k-1)*ny+1:ny*f)-L1(:,(k-1)*ny+1:ny*f)  zeros(n,(k-1)*ny);
      -L2(:,(k-1)*ny+1:ny*f)                      zeros(ny,(k-1)*ny)];
  if k == 1
      N(n+1:n+ny,1:ny) = eye(ny) + N(n+1:n+ny,1:ny);
  end
  N=N*X;%%式27
  tot = tot + kron(Qf((k-1)*nu+1:k*nu,:)',N);    
end

%   解最小二乘 (式28)
T=T(:);
sol=tot\T;

%    提取B和D
sol2=reshape(sol,(n+ny),nu);
D=sol2(1:ny,:);
B=sol2(ny+1:ny+n,:);

end
