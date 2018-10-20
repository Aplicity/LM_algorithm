clc;clear
% ����f(x) = x^2 ����Сֵʱx��ȡֵ

fx = @(x) (x^2)
Jf = @(x) (2*x)

x = 10;  % ��ֵx0
iter = 1000; % ��������
error = 1e-12;  % �����
u = 1;   % ��ʼlamda

record_x = [x];     % x�ĵ�����¼
record_lamda = [u]; % lamda�ĵ�����¼
record_grad = [0];  % �����ĵ�����¼

for i = 1:iter
    H = Jf(x)' * Jf(x) ; % Hessian ����
    G = H + u*eye(size(H));  % LM�㷨����Hessian���� ,����IΪdiag(1,1...)
    last = x;   % ��¼�ϴε�����xֵ
    while det(G) == 0  % ������Hessian���������
        u = u*4 ;       % ����uֵ��ֱ������Hessian��������
    end
    
    G = H + u*eye(size(H));   %����Hessian����
    r = (G)^(-1)* Jf(x)' * fx(x);   % ��������������
    x = x - r   % ����xֵ
    
    record_x(i+1) = x;
    record_grad(i) = r;
    
    if 0 < r < 0.25   % �����������̫С��ʱ������һ��������ʱ���ʵ��Ӵ��������
        u = u * 4;
        
    elseif 0.25 <= r <= 0.75 % �������������ʣ����ֲ���
        u = u;
    
    elseif r > 0.75 % �����������̫���ʱ������һ��������ʱ���ʵ���С�������룬����Խ�����ŵ�
        u = u / 2;
 
    end
    
    record_lamda(i+1) = u ;
    
    if r < 0    %��rk��0 ��˵������ֵ���������������½������Ʊ仯�ˣ������Ż���Ŀ���෴������ֹͣ����
        break
    end
    
    if abs(last - x) < error   % ���������ﵽ�趨ֵ��ֹͣ����
        break
    end
    
end
    
        
