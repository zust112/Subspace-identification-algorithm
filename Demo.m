%   ��ֵʾ��(pole)
clear
clc
close all
addpath(genpath('.'));

% SIMOϵͳ
a = [0.67 0.67 0 0;-0.67 0.67 0 0;0 0 -0.67 -0.67;0 0 0.67 -0.67];
b = [0.7123;1.9698;1.3171;-1.1752];
c = [-0.2324 1.0751 -0.5225 0.1130;1.4027 1.7543 -1.2159 2.1982];
d = [0.7139;-0.9632];
%%���漫��
GG_CN=zeros(4,100);
GG_ORT=zeros(4,100);
GG_PCA=zeros(4,100);

% % ���ؿ����������
for g=1:100
u=zeros(1000,1);
for k=1:1000
    MM=0;
     for i=1:10
       MM=MM+sin(0.3898*pi*i*k);%%�����ź����뼤��
    end
    u(k)=MM*rand;
end
u0=u(1);

 % % ��ɫ��������Ƶ��
f = 50;          % ��Ƶ�ź�Ƶ��Hz
T_pulse = 1000;    % ������
T = 1;% ��������
n_number = T_pulse * T;
k = 1: T :T_pulse;
s = cos(2.713*pi*f*k);

u=u+1*s';%%����ɫ����


y = dlsim(a, b, c, d, u);
y0=y(1,:);
y = y + [1*s',1*s'];%%%����ɫ����


% % % % % % % % % % % 
% % �׻��˲�
%�׻��˲�u
e1=u;
wu0=u0;%%�׻���ֵ
e0=u(1);%%��ɫ��ֵ
wu=zeros(1000,1);
wu(1)=0.5*wu0+e1(1)-e0;
for k=2:1000
    wu(k)=0.5*wu(k-1)+e1(k)-e1(k-1);
end
uCN=wu;%�׻�������

% %�׻��˲�y
e2=y;
wy11=y0(1);%%�׻���ֵ
wy22=y0(2);
e11=y(1,1);%%��ɫ��ֵ
e22=y(1,2);
wy=zeros(1000,2);
wy(1,1)=0.5*wy11+e2(1,1)-e11;
wy(1,2)=0.5*wy22+e2(1,2)-e22;
for k=2:1000
    wy(k,1)=0.5*wy(k-1,1)+e2(k,1)-e2(k-1,1);
    wy(k,2)=0.5*wy(k-1,2)+e2(k,2)-e2(k-1,2);
end
yCN=wy;%�׻�������

% % % % % % %��ʶ
% % 
[A_CN, B_CN, C_CN, D_CN] = SIMPCA_CN(yCN,uCN,5,5,4);
[A_PCA, B_PCA, C_PCA, D_PCA] = SIMPCA(y,u,5,5,4);
[A_ORT, B_ORT, C_ORT, D_ORT] = TORT_SIM(y,u,5,5,4);



% % % % % % % % ����ϵͳ
sys_pca = ss(A_PCA, B_PCA, C_PCA, D_PCA);
sys_CN = ss(A_CN, B_CN, C_CN, D_CN);
sys_ORT = ss(A_ORT, B_ORT, C_ORT, D_ORT);


% % % % % % % % % % % % % % % % % % % % % 
[p2,z2]=pzmap(sys_CN);
[p3,z3]=pzmap(sys_pca);
[p4,z4]=pzmap(sys_ORT);
GG_CN(:,g)=p2;
GG_PCA(:,g)=p3;
GG_ORT(:,g)=p4;

end

sys_real = ss(a,b,c,d);
[p1,z1]=pzmap(sys_real);
%%%2ORT-SIM����ͼ
figure(1)
set(figure(1),'position',[549 146.6 560 512.8]);
[~,hp1]=zplane(z1,p1); 
set(hp1,'markersize',20,'color','r','linewidth',1.5,'marker','+')
hold on
[~,hp2]=zplane([],GG_ORT);
set(hp2,'color','b','markersize',10,'linewidth',1);
xlabel('ʵ��')
ylabel('����')
axis([-1 1 -1 1])
set(gca,'FontSize',13,'position',[0.121 0.099 0.811 0.886]);
yticks([-1 -0.5 0 0.5 1])
% title('2ORT-SIM')


%%%SIMPCA����ͼ
figure(2)
set(figure(2),'position',[549 146.6 560 512.8]);
[~,hp1]=zplane(z1,p1); 
set(hp1,'markersize',20,'color','r','linewidth',1.5,'marker','+')
hold on
[~,hp3]=zplane([],GG_PCA);
set(hp3,'color','b','markersize',10,'linewidth',1);
axis([-1 1 -1 1])
set(gca,'FontSize',13,'position',[0.121 0.099 0.811 0.886]);
yticks([-1 -0.5 0 0.5 1])
% title('SIMPCA')
xlabel('ʵ��')
ylabel('����')

%%%SIMPCA-CN����ͼ
figure(3)
set(figure(3),'position',[549 146.6 560 512.8]);
[~,hp1]=zplane(z1,p1); 
set(hp1,'markersize',20,'color','r','linewidth',1.5,'marker','+')
hold on
[~,hp4]=zplane([],GG_CN);
set(hp4,'color','b','markersize',10,'linewidth',1);
axis([-1 1 -1 1])
set(gca,'FontSize',13,'position',[0.121 0.099 0.811 0.886]);
yticks([-1 -0.5 0 0.5 1])
% title('SIMPCA-CN')
xlabel('ʵ��')
ylabel('����')

