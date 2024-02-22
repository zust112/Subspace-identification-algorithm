%%%%%%%%%%%�ο����ף�Orthogonal projection based subspace identification against colored noise
function [A,B,C,D]=TORT_SIM(y,u,f,p,n)
%*******************************************
%   ��γ���˵��
%   y���������
%   u����������
%   f��δ��Hankel������ά��
%   p�ǹ�ȥHankel������ά��
%   n��ϵͳ����
%   A,B,C,D��״̬�ռ����



%**********�ж��������*********%
if nargin < 5
    n = [];
end

%**********��������תΪ������********%
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
    error('�������������������ͬ');
end
N = numy;%%��������

%����Hankel����
U=Hankel(p,f,u);
Y=Hankel(p,f,y);

%���ֹ�ȥ��δ��
Yp=Y(1:p*ny,:);
Yf=Y(p*ny+1:end,:);
Up=U(1:p*nu,:);
Uf=U(p*nu+1:end,:);

%%%%ͶӰ����ֵ�㷨��ͨ��QR�ֽ���ʵ��
%%%%��һ��ͶӰ����Yf����ͶӰ��Uf����������%%%%ʽ8
[Q0,L0]=qr([Uf;Yf]',0);
Q0=Q0';
L0=L0';
Ypro1=L0(f*nu+1:end,f*nu+1:end)*Q0(f*nu+1:end,:);

%%%%�ڶ���ͶӰ����������(ͶӰ��Up�пռ���)%%%%ʽ15
[Q1,L1]=qr([Up;Ypro1]',0);
Q1=Q1';
L1=L1';
% L11=L1(1:p*nu,1:p*nu);
L21=L1(p*nu+1:end,1:p*nu);
Ypro2=L21*Q1(1:p*nu,:);

%%%%�������ǳ����Լ������״�=max(p,f)����ͨ����������ͶӰ��������Ϊ0
%%%%SVD�ֽ�ʽ18
[P,S,V]=svd(Ypro2);

%%%%�״α�ʶ
if isempty(n)
    ss=diag(S);
    figure
    bar(ss);
    n=input('�״�Ϊ��');
end

if (f*ny-n) < 1
    error('����̫����û�����f̫С');
end
%%%%������չ�ɹ۾���
P1=P(:,1:n);
Fi=P1;

%%%%����ϵͳ����A��C
J1=[eye((f-1)*ny),zeros((f-1)*ny,ny)];
J2=[zeros((f-1)*ny,ny),eye((f-1)*ny)];
J3=[eye(ny),zeros(ny,(f-1)*ny)];
A=pinv(J1*Fi)*J2*Fi;
C=J3*Fi;

%%%%����ϵͳ����B��D
Fi_orc=orthc(Fi);%%������չ�ɹ۾����������

%��Yf�ҳ�Uf��α�棬���Fi��������%%%%ʽ22
M=Fi_orc'*Yf*pinv(Uf);%%ʽ23ά��Ϊ(fny-n)*fnu
L=Fi_orc';%%ʽ24ά��Ϊ(fny-n)*fny
%�������ǳ����Լ������״�Ϊf����ͨ������������Ϊ0
%����ʽ25
X=[eye(ny) zeros(ny,n);
   zeros((f-1)*ny,ny) J1*Fi];
Lhs=zeros((f*ny-n)*f,f*ny);
Mhs=zeros((f*ny-n)*f,nu);

for k=1:f
    Mhs( (k-1)*(f*ny-n)+1:k*(f*ny-n),: )=M( :,(k-1)*nu+1:k*nu );
    Lhs( (k-1)*(f*ny-n)+1:k*(f*ny-n),1:ny*(f-k+1) )=L( :,(k-1)*ny+1:end );
end

sol=pinv(Lhs*X)*Mhs;%%ԭʽ25�м���λ�ó��ִ���
D=sol(1:ny,:);
B=sol(ny+1:end,:);

end