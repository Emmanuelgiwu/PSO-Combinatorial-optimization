function v1=cv(c,v)
%常数和速度的乘积
length=round(c*size(v(1,:),2));
v1=v(:,1:length);
end