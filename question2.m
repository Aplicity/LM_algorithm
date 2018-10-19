clc;clear
fxy = @(x,y) ((1 - x)^2 + 100 * (y-x*x)^2)
Jf = @(x,y) ([2*(x-1)-400*x*(y-x*x),200*(y-x*x)])

x = 5; y = 5 % 初始化x，y
iter = 10000; % 迭代次数
error1 = 1e-6   % x的误差限
error2 = 1e-6   % y的误差限
u = 1   % 初始u

for i = 1:iter
    H = Jf(x,y)' .* Jf(x,y);  % Hessian 矩阵
    G = H + u*eye(size(H));  % LM算法修正Hessian矩阵 ,其中I为diag(1,1...)
    last_x = x; last_y = y;   % 记录上本次迭代的x,y值
    
    while det(G) == 0  % 当修正Hessian矩阵非正定
        u = u*4;      % 更新u值，直至修正Hessian矩阵正定
    end
    G = H + u*eye(size(H));   %修正Hessian矩阵
    r = (G)^(-1)* Jf(x,y)'* fxy(x,y);   % 计算迭代方向距离
    A = [x;y] - r  % 更新 x，y的值
    x = A(1)  % x的值
    y = A(2)  % y的值
    
    if 0 < r < 0.25   % 迭代方向距离太小的时候，在下一步迭代的时候适当加大迭代距离
        u = u * 4;
    end
    
    if r > 0.75 % 迭代方向距离太大的时候，在下一步迭代的时候适当减小迭代距离，避免越过最优点
        u = u / 2;
    end
    
    if 0.25 <= r <= 0.75 % 迭代方向距离合适，保持不变
        u = u;
    end
    
    if r < 0    %若rk≤0 ，说明函数值是向着上升而非下降的趋势变化了（与最优化的目标相反），则停止迭代
        break
    end
    
     if abs(last_x - x) < error1 && abs(last_y - y) < error2  % 当收敛误差达到设定值，停止迭代
         break
     end
 
end



