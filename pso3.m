clc,clear;
n=5;tic
%% ��ʼ����Ⱥ
N = 30;                         % ��ʼ��Ⱥ����
d = n;                          % �ռ�ά��,�������ĸ���n=3
ger = 10;                      % ����������     
% limit = [0, 20];                % ����λ�ò�������
%vlimit = [n,nchoosek(n,2)];               % �����ٶ�����,������Ըģ���n���ʱ��Ӧ������ֵ
w = 0.2;                        % ����Ȩ��
c1 = 1;                       % ����ѧϰ����
c2 = 1;                       % Ⱥ��ѧϰ���� 
%xm=zeros(1,n)               %����x
x=zeros(N,d);
perm=perms([1,2,3,4,5]);
for i=1:N
x(i,:) = perm(randi(factorial(d)),:);%��ʼ��Ⱥ��λ��

% v(i,:)={randi(d,1,2),randi(d,1,2),randi(d,1,2)}; % ��ʼ��Ⱥ���ٶ�  %������Ҫ��n�޸�

end
v=randi(n,N,2*n);           %�����v�ٶ��ǿɱ䳤��
xm = x;                          % ÿ���������ʷ���λ��
ym = perm(randi(d),:);                % ��Ⱥ����ʷ���λ��,������ɳ�ʼ��
fxm = zeros(1,N);               % ÿ���������ʷ�����Ӧ��
ft1=zeros(1,N);              %ÿ�������t1ֵ
fx=zeros(1,N);
fym = -inf;                      % ��Ⱥ��ʷ�����Ӧ��
fmt1=-inf;                      %���ŵ�t1��ȡֵ
%hold on
%plot(xm, f(xm), 'ro');title('��ʼ״̬ͼ');
%figure(2)
%% Ⱥ�����
iter = 1;
record = zeros(ger, 1);          % ��¼��
while iter <= ger
    for k=1:N
        ff=fitness(x(k,:));
        fx(k)=ff(2);
        ft1(k)=ff(1);
    end
     %fx = fitness(x) ; % ���嵱ǰ��Ӧ��   
     for i = 1:N      
        if fxm(i) < fx(i)
            fxm(i) = fx(i);     % ���¸�����ʷ�����Ӧ��
            ft1(i)=ft1(i);
            xm(i,:) = x(i,:);   % ���¸�����ʷ���λ��
        end 
     end
if fym < max(fxm)
        [fym, nmax] = max(fxm);   % ����Ⱥ����ʷ�����Ӧ��
        [mmm,index]=max(fxm);  %�ҵ����ֵ���±�
        fmt1=ft1(index);
        ym = xm(nmax, :);      % ����Ⱥ����ʷ���λ��
end
 
%������Ҫ���¶���+��-����
%������Ҫ��wȨ�ض����������

    v =[cv(w,v),  cv(rand,xminus(xm,x)), cv(rand, xminus(repmat(ym,N,1),x))]  ;% �ٶȸ���,repmat�½�����
    % �߽��ٶȴ���
%     v(v > vlimit(2)) = vlimit(2);
%     v(v < vlimit(1)) = vlimit(1);
    x = xvplus(x,v,N);% λ�ø���
%     % �߽�λ�ô���
%     x(x > limit(2)) = limit(2);
%     x(x < limit(1)) = limit(1);
    record(iter) = fym;%���ֵ��¼
%     x0 = 0 : 0.01 : 20;
%     plot(x0, f(x0), 'b-', x, f(x), 'ro');title('״̬λ�ñ仯')
%     pause(0.1)
    iter = iter+1;
end
% figure(3);plot(record);title('��������')
% x0 = 0 : 0.01 : 20;
% figure(4);plot(x0, f(x0), 'b-', x, f(x), 'ro');title('����״̬λ��')
% disp(['���ֵ��',num2str(fym)]);
% disp(['����ȡֵ��',num2str(ym)]);
toc