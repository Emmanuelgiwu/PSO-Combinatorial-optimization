clc,clear;
n=5;tic
%% 初始化种群
N = 30;                         % 初始种群个数
d = n;                          % 空间维数,即部件的个数n=3
ger = 10;                      % 最大迭代次数     
% limit = [0, 20];                % 设置位置参数限制
%vlimit = [n,nchoosek(n,2)];               % 设置速度限制,这里可以改，当n变大时，应调整此值
w = 0.2;                        % 惯性权重
c1 = 1;                       % 自我学习因子
c2 = 1;                       % 群体学习因子 
%xm=zeros(1,n)               %最优x
x=zeros(N,d);
perm=perms([1,2,3,4,5]);
for i=1:N
x(i,:) = perm(randi(factorial(d)),:);%初始种群的位置

% v(i,:)={randi(d,1,2),randi(d,1,2),randi(d,1,2)}; % 初始种群的速度  %这里需要随n修改

end
v=randi(n,N,2*n);           %这里的v速度是可变长的
xm = x;                          % 每个个体的历史最佳位置
ym = perm(randi(d),:);                % 种群的历史最佳位置,随机生成初始。
fxm = zeros(1,N);               % 每个个体的历史最佳适应度
ft1=zeros(1,N);              %每个个体的t1值
fx=zeros(1,N);
fym = -inf;                      % 种群历史最佳适应度
fmt1=-inf;                      %最优的t1的取值
%hold on
%plot(xm, f(xm), 'ro');title('初始状态图');
%figure(2)
%% 群体更新
iter = 1;
record = zeros(ger, 1);          % 记录器
while iter <= ger
    for k=1:N
        ff=fitness(x(k,:));
        fx(k)=ff(2);
        ft1(k)=ff(1);
    end
     %fx = fitness(x) ; % 个体当前适应度   
     for i = 1:N      
        if fxm(i) < fx(i)
            fxm(i) = fx(i);     % 更新个体历史最佳适应度
            ft1(i)=ft1(i);
            xm(i,:) = x(i,:);   % 更新个体历史最佳位置
        end 
     end
if fym < max(fxm)
        [fym, nmax] = max(fxm);   % 更新群体历史最佳适应度
        [mmm,index]=max(fxm);  %找到最大值的下标
        fmt1=ft1(index);
        ym = xm(nmax, :);      % 更新群体历史最佳位置
end
 
%这里需要重新定义+和-算子
%这里需要对w权重定义调整函数

    v =[cv(w,v),  cv(rand,xminus(xm,x)), cv(rand, xminus(repmat(ym,N,1),x))]  ;% 速度更新,repmat新建矩阵
    % 边界速度处理
%     v(v > vlimit(2)) = vlimit(2);
%     v(v < vlimit(1)) = vlimit(1);
    x = xvplus(x,v,N);% 位置更新
%     % 边界位置处理
%     x(x > limit(2)) = limit(2);
%     x(x < limit(1)) = limit(1);
    record(iter) = fym;%最大值记录
%     x0 = 0 : 0.01 : 20;
%     plot(x0, f(x0), 'b-', x, f(x), 'ro');title('状态位置变化')
%     pause(0.1)
    iter = iter+1;
end
% figure(3);plot(record);title('收敛过程')
% x0 = 0 : 0.01 : 20;
% figure(4);plot(x0, f(x0), 'b-', x, f(x), 'ro');title('最终状态位置')
% disp(['最大值：',num2str(fym)]);
% disp(['变量取值：',num2str(ym)]);
toc