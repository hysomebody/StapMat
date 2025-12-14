function k = BeamLocalStiff(E, G, A, Iy, Iz, J, L)
% BEAMLOCALSTIFF 计算3D欧拉-伯努利梁单元的局部刚度矩阵 (12x12)
% 输入:
%   E, G: 弹性模量, 剪切模量
%   A: 截面积
%   Iy, Iz: 绕局部y轴和z轴的惯性矩
%   J: 抗扭常数
%   L: 单元长度
% 输出:
%   k: 12x12 局部刚度矩阵

    k = zeros(12, 12);
    
    % --- 刚度系数 ---
    A_L = A * E / L;
    GJ_L = G * J / L;
    
    % 绕 Z 轴弯曲 (在 xy 平面内) -> 相关参数 Iz
    EIz_L3 = 12 * E * Iz / L^3;
    EIz_L2 = 6 * E * Iz / L^2;
    EIz_L  = 4 * E * Iz / L;
    EIz_L_half = 2 * E * Iz / L;
    
    % 绕 Y 轴弯曲 (在 xz 平面内) -> 相关参数 Iy
    EIy_L3 = 12 * E * Iy / L^3;
    EIy_L2 = 6 * E * Iy / L^2;
    EIy_L  = 4 * E * Iy / L;
    EIy_L_half = 2 * E * Iy / L;
    
    % --- 填充矩阵 ---
    % 自由度顺序: u1, v1, w1, thx1, thy1, thz1, u2, v2, w2, thx2, thy2, thz2
    % 索引对应:    1   2   3    4     5     6    7   8   9    10    11    12
    
    % 1. 轴向 (u - x轴)
    k(1,1) = A_L;   k(1,7) = -A_L;
    k(7,1) = -A_L;  k(7,7) = A_L;
    
    % 2. 扭转 (theta_x - 绕x轴)
    k(4,4) = GJ_L;   k(4,10) = -GJ_L;
    k(10,4) = -GJ_L; k(10,10) = GJ_L;
    
    % 3. 弯曲 XY 平面 (v, theta_z) -> 使用 Iz
    % 剪力 v
    k(2,2) = EIz_L3;   k(2,8) = -EIz_L3;
    k(8,2) = -EIz_L3;  k(8,8) = EIz_L3;
    % 弯矩 Mz (theta_z)
    k(6,6) = EIz_L;    k(6,12) = EIz_L_half;
    k(12,6) = EIz_L_half; k(12,12) = EIz_L;
    % 耦合 v - theta_z
    k(2,6) = EIz_L2;   k(2,12) = EIz_L2;
    k(6,2) = EIz_L2;   k(12,2) = EIz_L2;
    k(8,6) = -EIz_L2;  k(8,12) = -EIz_L2;
    k(6,8) = -EIz_L2;  k(12,8) = -EIz_L2;
    
    % 4. 弯曲 XZ 平面 (w, theta_y) -> 使用 Iy
    % 注意：根据右手定则，Y轴弯曲的符号通常与Z轴弯曲有差异，具体取决于坐标系定义
    % 这里采用标准 FEM 教科书符号
    % 剪力 w
    k(3,3) = EIy_L3;   k(3,9) = -EIy_L3;
    k(9,3) = -EIy_L3;  k(9,9) = EIy_L3;
    % 弯矩 My (theta_y)
    k(5,5) = EIy_L;    k(5,11) = EIy_L_half;
    k(11,5) = EIy_L_half; k(11,11) = EIy_L;
    % 耦合 w - theta_y (注意符号变化)
    k(3,5) = -EIy_L2;  k(3,11) = -EIy_L2;
    k(5,3) = -EIy_L2;  k(11,3) = -EIy_L2;
    k(9,5) = EIy_L2;   k(9,11) = EIy_L2;
    k(5,9) = EIy_L2;   k(11,9) = EIy_L2;
    
end