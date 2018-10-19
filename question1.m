clc;clear
fx = @(x) (x^2)
Jf = @(x) (2*x)
Hf = @(x) (2)

x = 10;  % 初值x0
iter = 1000; % 迭代次数
error = 1e-12   % 误差限
u = 1   % 初始u

for i = 1:iter
    H = Jf(x)' * Jf(x) ; % Hessian 矩阵
    G = H + u*eye(size(H));  % LM算法修正Hessian矩阵 ,其中I为diag(1,1...)
    last = x;   % 记录上本次迭代的x值
    while det(G) == 0  % 当修正Hessian矩阵非正定
        u = u*4        % 更新u值，直至修正Hessian矩阵正定
    end
    G = H + u*eye(size(H));   %修正Hessian矩阵
    r = (H)^(-1)* Jf(x)' * fx(x);   % 计算迭代方向距离
    x = x - r   % 更新x值
    
    
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
    
    if abs(last - x) < error   % 当收敛误差达到设定值，停止迭代
        break
    end
    
end


