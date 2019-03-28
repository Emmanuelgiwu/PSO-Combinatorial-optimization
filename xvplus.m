function x1=xvplus(x,v,N)
%位置加速度得到新的位置,注意这里需要修改50=N，种群规模个数
j=1;
x;v;
for i=1:N
    while j<size(v(1,:),2)   %对于v遍历所有的列，进行交换算子
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