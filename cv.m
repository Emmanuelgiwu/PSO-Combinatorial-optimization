function v1=cv(c,v)
%�������ٶȵĳ˻�
length=round(c*size(v(1,:),2));
v1=v(:,1:length);
end