function M_global = BeamMass(RHO, A, Iy, Iz, J, L, T)

% 计算3D欧拉-伯努利梁单元的一致质量矩阵 (12x12)
% 并将其转换到全局坐标系
%
% 输入:
%   RHO: 材料密度 (Mass density)
%   A:   截面积
%   Iy, Iz: 惯性矩 (用于考虑转动惯量，但在欧拉梁中通常主要用 A)
%   J:   抗扭常数 (用于计算极惯性矩 Ip)
%   L:   单元长度
%   T:   12x12 变换矩阵 (从 BeamStiff 或外部传入)
%
% 输出:
%   M_global: 全局坐标系下的质量矩阵

    % 初始化局部质量矩阵
    m = zeros(12, 12);
    
    % --- 系数定义 ---
    % 1. 平动质量系数 (Total Mass)
    m_const = RHO * A * L;
    
    % 2. 转动/扭转质量系数 (Polar Mass Moment of Inertia)
    % 对于非圆截面，Ip 近似取 J (或 Iy+Iz)，这里用 J 近似
    % 如果是圆截面 Ip = J; 
    rx2 = J / A; % 回转半径平方近似
    m_rot = RHO * A * L * rx2; 

    % ==========================================
    % 填充局部质量矩阵 (Euler-Bernoulli Consistent Mass)
    % ==========================================
    
    % --- 1. 轴向振动 (u) ---
    % 形函数为线性，系数为 1/6
    % 自由度: 1 (u1), 7 (u2)
    c1 = m_const / 3.0;
    c2 = m_const / 6.0;
    
    m(1,1) = c1; m(1,7) = c2;
    m(7,1) = c2; m(7,7) = c1;
    
    % --- 2. 扭转振动 (theta_x) ---
    % 形函数为线性，类似于轴向
    % 自由度: 4 (thx1), 10 (thx2)
    ct1 = m_rot / 3.0;
    ct2 = m_rot / 6.0;
    
    m(4,4) = ct1;   m(4,10) = ct2;
    m(10,4) = ct2;  m(10,10) = ct1;
    
    % --- 3. XY平面弯曲 (v, theta_z) ---
    % 对应位移: v (2, 8), 转角: theta_z (6, 12)
    % 标准一致质量矩阵系数 (420)
    
    val_156 = 156 * m_const / 420;
    val_54  = 54  * m_const / 420;
    val_22  = 22  * m_const * L / 420;
    val_13  = 13  * m_const * L / 420;
    val_4   = 4   * m_const * L^2 / 420;
    val_3   = 3   * m_const * L^2 / 420; % 注意这里是 -3 的绝对值部分
    
    % 对角块 1 (Node I)
    m(2,2) = val_156;
    m(6,6) = val_4;
    m(2,6) = val_22;  m(6,2) = val_22;
    
    % 对角块 2 (Node J)
    m(8,8) = val_156;
    m(12,12) = val_4;
    m(8,12) = -val_22; m(12,8) = -val_22; % 注意符号变化
    
    % 耦合块 (Node I - Node J)
    m(2,8) = val_54;  m(8,2) = val_54;
    m(2,12) = -val_13; m(12,2) = -val_13;
    
    m(6,8) = val_13;  m(8,6) = val_13;
    m(6,12) = -val_3; m(12,6) = -val_3;
    
    % --- 4. XZ平面弯曲 (w, theta_y) ---
    % 对应位移: w (3, 9), 转角: theta_y (5, 11)
    % 系数与 XY 平面相同，但注意转角符号约定
    % 此时 w 对应 v， -theta_y 对应 theta_z
    
    % 对角块 1
    m(3,3) = val_156;
    m(5,5) = val_4;
    m(3,5) = -val_22; m(5,3) = -val_22; % 符号与XY面相反
    
    % 对角块 2
    m(9,9) = val_156;
    m(11,11) = val_4;
    m(9,11) = val_22; m(11,9) = val_22;
    
    % 耦合块
    m(3,9) = val_54;  m(9,3) = val_54;
    m(3,11) = val_13; m(11,3) = val_13;
    
    m(5,9) = -val_13; m(9,5) = -val_13;
    m(5,11) = -val_3; m(11,5) = -val_3;
    
    % 坐标变换
    M_global = T' * m * T;

end