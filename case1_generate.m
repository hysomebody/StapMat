% Generate_Verification_Case.m
% 修正版：修复了桁架节点转动自由度导致的矩阵奇异问题

clear; clc;

filename = 'Data/verification_100.in';
[pathstr, ~, ~] = fileparts(filename);
if ~exist(pathstr, 'dir'), mkdir(pathstr); end

fid = fopen(filename, 'w');

%% 1. 参数定义
L_beam = 10.0;          
n_beam_elems = 100;     
dx = L_beam / n_beam_elems;
L_truss = 2.0;          

% 材料与截面
E = 2.0e11;             
rho = 7850;             
G = 8.0e10;             
b = 0.1; h = 0.1;
A_beam = b*h;
Iz_beam = (b*h^3)/12;   
Iy_beam = (h*b^3)/12;
J_beam = 0.0001;        
A_truss = 0.001;        

% 载荷
P_load = -1000.0;       

%% 2. 理论解计算
disp_beam_tip = (P_load * L_beam^3) / (3 * E * Iz_beam);
disp_truss_ext = (P_load * L_truss) / (E * A_truss);
disp_total = disp_beam_tip + disp_truss_ext;

fprintf('======================================================\n');
fprintf('  STAPMAT 动力学验证算例 (修正版)\n');
fprintf('  理论静力学参考解:\n');
fprintf('  1. 梁末端(Node %d) Y位移: %12.6e m\n', n_beam_elems+1, disp_beam_tip);
fprintf('  2. 桁架底端(Node %d) Y位移: %12.6e m\n', n_beam_elems+2, disp_total);
fprintf('======================================================\n');

%% 3. 写输入文件

fprintf(fid, 'Verification Case: Fixed Truss Rotations\n');

% Control Line
% NUMNP, NUMEG, NLCASE, MODEX, ANALYSIS_TYPE
NUMNP = n_beam_elems + 2; 
NUMEG = 2;                
NLCASE = 1;
MODEX = 1;
ANALYSIS_TYPE = 1;        
fprintf(fid, '%5d %5d %5d %5d %5d\n', NUMNP, NUMEG, NLCASE, MODEX, ANALYSIS_TYPE);

% Dynamic Params
% DT, NSTEPS, RHO_INF, ALPHA, BETA
DT = 0.01;
NSTEPS = 500;             
RHO_INF = 0.8;
ALPHA = 5.0;              
BETA = 0.0;
fprintf(fid, '%10.4f %10d %10.2f %10.4f %10.4f\n', DT, NSTEPS, RHO_INF, ALPHA, BETA);

% Nodal Data
% Node 1 (Root): Fixed
fprintf(fid, '%5d %5d %5d %5d %5d %5d %5d %15.6f %15.6f %15.6f\n', ...
    1, 1, 1, 1, 1, 1, 1, 0.0, 0.0, 0.0);

% Beam Nodes (2 to N+1): All Free
for i = 1:n_beam_elems
    node_id = i + 1;
    x_coord = i * dx;
    fprintf(fid, '%5d %5d %5d %5d %5d %5d %5d %15.6f %15.6f %15.6f\n', ...
        node_id, 0, 0, 0, 0, 0, 0, x_coord, 0.0, 0.0);
end

% [修正点] Truss Bottom Node (N+2)
% 必须固定转动自由度 (BC 4,5,6 = 1)，只保留平动 (0 0 0 1 1 1)
node_truss_bottom = n_beam_elems + 2;
x_truss = L_beam;
y_truss = -L_truss; 
fprintf(fid, '%5d %5d %5d %5d %5d %5d %5d %15.6f %15.6f %15.6f\n', ...
    node_truss_bottom, 0, 0, 0, 1, 1, 1, x_truss, y_truss, 0.0);

% Load Data
fprintf(fid, '   1    1\n');
fprintf(fid, '%5d %5d %15.6f\n', node_truss_bottom, 2, P_load);

% Element Group 1: BEAM
fprintf(fid, '   2 %5d    1\n', n_beam_elems);
fprintf(fid, '   1 %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n', ...
    E, G, A_beam, Iy_beam, Iz_beam, J_beam, rho);
for i = 1:n_beam_elems
    n1 = i;
    n2 = i + 1;
    fprintf(fid, '%5d %5d %5d    1    0.0    1.0    0.0\n', i, n1, n2);
end

% Element Group 2: TRUSS
fprintf(fid, '   1    1    1\n');
fprintf(fid, '   1 %10.3e %10.3e %10.3e\n', E, A_truss, rho);
fprintf(fid, '   1 %5d %5d    1\n', n_beam_elems+1, n_beam_elems+2);

fclose(fid);
fprintf('修正后的输入文件生成成功: %s\n', filename);