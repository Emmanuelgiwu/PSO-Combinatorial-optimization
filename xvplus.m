function x1=xvplus(x,v,N)
%λ�ü��ٶȵõ��µ�λ��,ע��������Ҫ�޸�50=N����Ⱥ��ģ����
j=1;
x;v;
for i=1:N
    while j<size(v(1,:),2)   %����v�������е��У����н�������
        m=find(x(i,:)==v(i,j));
        q=find(x(i,:)==v(i,j+1));
        k=x(i,m);
        l=x(i,q);
        x(i,m)=l;
        x(i,q)=k;
        j=j+1;
    end
end
x1=x;
end