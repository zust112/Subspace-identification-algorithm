%%%%%%%%%%%参考文献：Orthogonal projection based subspace identification against colored noise
function [A,B,C,D]=TORT_SIM(y,u,f,p,n)
%*******************************************
%   入参出参说明
%   y是输出数据
%   u是输入数据
%   f是未来Hankel矩阵行维度
%   p是过去Hankel矩阵行维度
%   n是系统阶数
%   A,B,C,D是状态空间矩阵



%**********判断输入参数*********%
if nargin < 5
    n = [];
end

%**********检查参数并转为行向量********%
[ny,numy]=size(y);
if numy < ny
    y=y';
    [ny,numy]=size(y);
end
[nu,numu]=size(u);
if numu < nu
    u=u';
    [nu,numu]=size(u);
end 
if numy~=numu
    error('输入输出样本数必须相同');
end
N = numy;%%总样本数

%构建Hankel矩阵
U=Hankel(p,f,u);
Y=Hankel(p,f,y);

%划分过去与未来
Yp=Y(1:p*ny,:);
Yf=Y(p*ny+1:end,:);
Up=U(1:p*nu,:);
Uf=U(p*nu+1:end,:);

%%%%投影的数值算法可通过QR分解来实现
%%%%第一次投影，将Yf正交投影到Uf的正交补上%%%%式8
[Q0,L0]=qr([Uf;Yf]',0);
Q0=Q0';
L0=L0';
Ypro1=L0(f*nu+1:end,f*nu+1:end)*Q0(f*nu+1:end,:);

%%%%第二次投影，消除噪声(投影到Up行空间上)%%%%式15
[Q1,L1]=qr([Up;Ypro1]',0);
Q1=Q1';
L1=L1';
% L11=L1(1:p*nu,1:p*nu);
L21=L1(p*nu+1:end,1:p*nu);
Ypro2=L21*Q1(1:p*nu,:);

%%%%若输入是持续性激励，阶次=max(p,f)，则通过上述两次投影后噪声项为0
%%%%SVD分解式18
[P,S,V]=svd(Ypro2);

%%%%阶次辨识
if isempty(n)
    ss=diag(S);
    figure
    bar(ss);
    n=input('阶次为？');
end

if (f*ny-n) < 1
    error('阶数太大或用户索引f太小');
end
%%%%估计扩展可观矩阵
P1=P(:,1:n);
Fi=P1;

%%%%估计系统矩阵A和C
J1=[eye((f-1)*ny),zeros((f-1)*ny,ny)];
J2=[zeros((f-1)*ny,ny),eye((f-1)*ny)];
J3=[eye(ny),zeros(ny,(f-1)*ny)];
A=pinv(J1*Fi)*J2*Fi;
C=J3*Fi;

%%%%估计系统矩阵B和D
Fi_orc=orthc(Fi);%%计算扩展可观矩阵的正交补

%对Yf右乘Uf的伪逆，左乘Fi的正交补%%%%式22
M=Fi_orc'*Yf*pinv(Uf);%%式23维度为(fny-n)*fnu
L=Fi_orc';%%式24维度为(fny-n)*fny
%若输入是持续性激励，阶次为f，则通过上述噪声项为0
%构建式25
X=[eye(ny) zeros(ny,n);
   zeros((f-1)*ny,ny) J1*Fi];
Lhs=zeros((f*ny-n)*f,f*ny);
Mhs=zeros((f*ny-n)*f,nu);

for k=1:f
    Mhs( (k-1)*(f*ny-n)+1:k*(f*ny-n),: )=M( :,(k-1)*nu+1:k*nu );
    Lhs( (k-1)*(f*ny-n)+1:k*(f*ny-n),1:ny*(f-k+1) )=L( :,(k-1)*ny+1:end );
end

sol=pinv(Lhs*X)*Mhs;%%原式25中计算位置出现错误
D=sol(1:ny,:);
B=sol(ny+1:end,:);

end