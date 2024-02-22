%%%%%%%%%% Hankel����Ĺ��� %%%%%%%%%%
%%%�����û��Զ���Ĺ�ȥ��δ��ʱ�̳���p,f�Լ�������N
%%%���������yΪ��
function [H]=Hankel(p,f,y)

[ny,N]=size(y);
H=zeros((p+f)*ny,N-f-p+1);
for i=1:N-f-p+1
    for j=1:p+f
        H((j-1)*ny+1:j*ny,i)=y(:,j+i-1);
    end
end

end
