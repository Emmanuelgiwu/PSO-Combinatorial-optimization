function fitval=fitness(xx)
global a;
 a=xx;

 % define parameters first
% Mutation function for constrained minimization
% options = optimoptions('ga','MutationFcn',@mutationadaptfeasible);
%  x = ga(@(x)func_obj(x),2,[],[],[],[],[1;2],[50;50;],@(x)func_cons(x),options)
[x,fval,exitflag,output]=fmincon(@func_obj,[2 4],[],[],[],[],zeros(2,1),[],@func_cons);
if  exitflag<0
    fitval=[0,0];
else
 fitval=x
end

function f=func_obj(t)
f=-t(2);
end
function [c, ceq]=func_cons(t)
k=5;n=5;theta=30;
% x=[3,2,1];
RR=0.2; %lowest reliability required
sigma=zeros(1,n);
mu=zeros(1,n);
for i =1:1:n
    mu(i)=0.5.*i;
    sigma(i)=0.5*i;
end
Re=func_Relability(k,n,t,theta,mu,sigma);
c(1)=t(1)-t(2)+2;
c(2)=RR-Re(2);
c(3)=RR-Re(1);
c(3)=t(2)-100;
ceq=[];
 end
function [Re]=func_Relability(k,n,t,theta,mu,sigma)
%WD是即时的退化矩阵
%返回值是t1.t2时刻的两次系统可靠性值R=[R1.R2]
%Y=（Y1,Y2）是系统的t1,t2的退化水平

Y=zeros(1,2);
Re=zeros(1,2);
[WD,NN]=func_degradation( n,t,mu,sigma);
[a1,index1]=sort(WD(:,NN(1))); %对于k n系统排序，取出退化水平处于第k个的部件退化值
b1=index1(k);%b 是退化水平排k的位置
Y(1)=a1(b1);%系统退化水平等于b位置的退化水平

%判定X是否已经超过theta，若超过，R=0
if theta>Y(1)
   Re(1)=func_IG((theta-Y(1))./mu(b1),((theta-Y(1))./sigma(b1)).^2,t(1))
else
   Re(1)=0;
end

[a2,index2]=sort(WD(:,NN(2))); %对于k n系统排序，取出退化水平处于第k个的部件退化值
b2=index2(k);
Y(2)=a2(b2);

if theta>Y(2)
  Re(2)=func_IG((theta-Y(2))./mu(b2),((theta-Y(2))./sigma(b2)).^2,t(2))
else
   Re(2)=0;
end



end
function [ W,NN] = func_degradation( n,t,mu,sigma)
N1=1000;
randn('state',1000)          % set the state of randn
 
  N11=ceil(N1/t(1));
   dt1=1/N1;
W1=zeros(n,N11);
for i =1:1:n
  dW = normrnd(mu(i).*dt1,sigma(i).*sqrt(dt1),1,N11);   % increments
  W1(i,:) = cumsum(dW);            % cumulative sum

%plot([0:dt:T],[0,W(1,:)],'r-') ;  % plot W against t
end

N2=5000;
 randn('state',1000)          % set the state of randn
 
  if (t(2)-t(1))/N2 <=0
      dt2=1;N22=1; 
  else
   
   N22=ceil(N2/(t(2)-t(1)));
   dt2=1/N22;
  end
  
W2=zeros(n,N22);
% EW=zeros(1,3);  %这里的EW储存了W1最后的退化水平
% 
% for i=1:1:n
%  b=a(i);
% EW(i)=W1(b,N11)
% end

for i =1:1:n
  dW = normrnd(mu(i).*dt2,sigma(i).*sqrt(dt2),1,N22);   % increments  
  W2(i,:) = cumsum(dW);% cumulative sum 
  W2(i,:)=W2(i,:)+W1(a(i),N11);
end


W=zeros(n,N11+N22);
% for i=1:1:n
%    W(i,:)=[W1(i,:),W2(i,:)] ;
% end
W=[W1,W2];
NN=[N11,N11+N22];
%  plot([0:dt1:tao1],[0,W1(1,:)],'r-')
%  hold on;  % plot W against t
%   plot([tao1:dt2:tao2],[0,W2(2,:)],'g-')
end
function[y]=func_IG(a1,a2,m)
pd = makedist('InverseGaussian','mu',a1,'lambda',a2);
y = 1-cdf(pd,m);


end

end
