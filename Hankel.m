%%%%%%%%%% Hankel矩阵的构建 %%%%%%%%%%
%%%基于用户自定义的过去、未来时刻长度p,f以及样本数N
%%%以输出序列y为例
function [H]=Hankel(p,f,y)

[ny,N]=size(y);
H=zeros((p+f)*ny,N-f-p+1);
for i=1:N-f-p+1
    for j=1:p+f
        H((j-1)*ny+1:j*ny,i)=y(:,j+i-1);
    end
end

end
