function [v]=GS(A)

v(:,1)=A(:,1)/norm(A(:,1));%��һ��
[Ahang,Alie]=size(A); %������к���
for k=2:Alie %����j����������
    res(:,1)=A(:,k);
    for i=1:k-1%��ȥ�������������������ϵ�ͶӰ
        res=res-v(:,i)'*res*v(:,i);
    end
    v(:,k)=res/norm(res);%��һ��
end
end