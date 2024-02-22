function [A,B,C,D]=SIMPCA_CN(y,u,f,p,n)
%*******************************************
%   ��γ���˵��
%   y�������������
%   u��������������
%   f��δ��Hankel������ά��
%   p�ǹ�ȥHankel������ά��
%   n��ϵͳ��������Ϊ�գ�
%   A,B,C,D�Ǳ�ʶ�õ���״̬�ռ����
%*******************************************

%   �ж����������תΪ������
if nargin < 5
    n = [];
end
%   nu,nyΪ���������������smΪ������
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
%   ����Hankel����
U=Hankel(p,f,u);
Y=Hankel(p,f,y);

%********************************************
%   ���ݾ��󻮷�
Yp=Y(1:p*ny,:);
Yf=Y(p*ny+1:(p+f)*ny,:);
Up=U(1:p*nu,:);
Uf=U(p*nu+1:(p+f)*nu,:);

%********************************************
%   ��������(ʽ19-20)
ZQR = [Yp;Up;Uf;Yf];
[Q,L]= qr(ZQR',0);%%�������ת��ʹ�����QR�еľ��÷ֽ�ѹ�������С����߼����ٶ�
Q=Q';
L=L';
L44 = L((nu+ny)*p+nu*f+1:end,(nu+ny)*p+nu*f+1:end);
Q4 = Q((nu+ny)*p+nu*f+1:end,:);
Ef=L44*Q4;
Yfb=Yf-Ef;

%********************************************
%   PCA����ֵ��ʹ��svd���㣩
[P,S,~]=svd([Yfb;Uf]);%%ʽ9

%********************************************
%   �״α�ʶ
if (isempty(n))
    ss=diag(S);
    figure
    bar(ss);
    n=input('�״�Ϊ��');
end
% ss=diag(S);

%********************************************
%   ��ȡ�������غɾ���
Pr=P(:,(nu*f)+n+1:end);%%ʽ10

%********************************************
%   �������(ʽ11)
Pry=Pr(1:ny*f,:);
Pru=Pr(ny*f+1:end,:);

% Pre=Pr(ny*(f+p)+1:end,:);

%********************************************
%******************����������**************%
%   ����һ
% Fi=null(Pry');%%Pry'����Ϊ���Ⱦ��󣬷�������ʱ������bug

%   ������
% Fi = orthcomp(Pry); %����������

%   ������
[mm,nn]=size(Pry);
if (mm < nn)
    Pry=Pry';
    [mm,nn]=size(Pry);
end
[QQ,~]=qr(Pry);
Fi=QQ(:,nn+1:mm);

Fi=Fi(:,1:n);

%********************************************
%   ��ȡA��C (ʽ21)
C=Fi(1:ny,:);
Fi1=Fi(1:ny*(f-1),:); %%ȥ����ny��
A=pinv(Fi1)*Fi(ny+1:ny*f,:);

%********************************************
%   ԭ������ȡB��D
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
%   �Ľ�����

%   �ض�Ϊ����
L=L(:,1:(f+p)*(ny+nu));

%   ��A��C�����¼���Fi��Fi1 (ʽ22)
Fi=C;
for k=2:f
   Fi((k-1)*ny+1:k*ny,:)=Fi((k-2)*ny+1:(k-1)*ny,:)*A;  
end
Fi1=Fi(1:ny*(f-1),:);%%ȥ����L�е���չ�ɹ�
Fi_inv=pinv(Fi);%%��չ�ɹ�α��
Fi1_inv=pinv(Fi1);%%ȥ����L�е���չ�ɹ�α��

%   ���µ�Fi��Fi1����ȷ�����Է����� (ʽ24)
Zi=L((nu+ny)*p+nu*f+1:end,1:(nu+ny)*p+nu*f);%%Zi=[R41 R42 R43]
Zi1=L((nu+ny)*p+nu*f+1+ny:end,1:(nu+ny)*p+nu*f+ny);%%Z(i+l)
Rhs=[Fi_inv*Zi          zeros(n,ny);
     L((nu+ny)*p+1:(nu+ny)*p+nu*f,1:(nu+ny)*p+nu*f+ny)];
 
Lhs=[Fi1_inv*Zi1;
    L((nu+ny)*p+nu*f+1:(nu+ny)*p+nu*f+ny,1:(nu+ny)*p+nu*f+ny)];

%   ����T��Qf
T=Lhs-[A;C]*Rhs(1:n,:);%%ʽ23
T=T(:,1:(nu+ny)*p);
Qf=L((nu+ny)*p+1:(nu+ny)*p+nu*f,1:(nu+ny)*p);%%δ������

%   ��������B��D���Է�����
L1=A*Fi_inv;
L2=C*Fi_inv;
M=[zeros(n,ny),Fi1_inv];
X=[eye(ny),zeros(ny,n);zeros(ny*(f-1),ny),Fi1];%%ʽ26

%   ����N��kron�ڻ�
tot=0;
for k=1:f
   N=[M(:,(k-1)*ny+1:ny*f)-L1(:,(k-1)*ny+1:ny*f)  zeros(n,(k-1)*ny);
      -L2(:,(k-1)*ny+1:ny*f)                      zeros(ny,(k-1)*ny)];
  if k == 1
      N(n+1:n+ny,1:ny) = eye(ny) + N(n+1:n+ny,1:ny);
  end
  N=N*X;%%ʽ27
  tot = tot + kron(Qf((k-1)*nu+1:k*nu,:)',N);    
end

%   ����С���� (ʽ28)
T=T(:);
sol=tot\T;

%    ��ȡB��D
sol2=reshape(sol,(n+ny),nu);
D=sol2(1:ny,:);
B=sol2(ny+1:ny+n,:);

end
