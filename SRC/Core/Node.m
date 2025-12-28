% 
% 功能：储存节点信息，管理坐标和 6 个自由度的边界条件。
% 
% Call procedures:
%   None
%
% Called by:
%   ./Domain.m - Used to create node instances during input.
% Programmed by: 宋博,程云志

classdef Node < handle
    properties (Constant)
        NDF = 6; %  X, Y, Z, RX, RY, RZ
    end
    
    properties
        ID      
        XYZ     
        BCode         % 0 = Free, 1 = Fixed
        Displacement  % 存储节点的 6 自由度位移向量
        Temperature   % 节点温度，用于热应力
    end
    
    methods
        function obj = Node()
            obj.ID = 0;
            obj.XYZ = zeros(3, 1);
            obj.BCode = zeros(obj.NDF, 1);
            obj.Displacement = zeros(obj.NDF, 1);  % 初始化为 0
            obj.Temperature = 0.0; % 初始化温度
        end
        

        function Read(obj, fid, expectedID)
            % 循环以跳过空行和纯注释行
            data = [];
            while isempty(data)
                lineStr = fgetl(fid);
                if ~ischar(lineStr) % 检查是否读到了文件末尾
                    error('Unexpected End of File reading Node %d', expectedID);
                end
                
                data = str2num(lineStr); %#ok<ST2NM>
                % 如果 data 为空，说明是空行或纯注释行，循环继续读取下一行
                % 如果 data 不为空，说明读到了有效数据，跳出循环
            end
            
            if isempty(data)
                error('Error reading node data for expected ID: %d', expectedID);
            end
            
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Node ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            
            % 根据数据列数判断格式
            colCount = length(data);
            
            if colCount == 7 
                % 格式 A: 3 BCs (Truss单元)
                % ID, BC1, BC2, BC3, X, Y, Z
                
                obj.BCode(1:3) = round(data(2:4));
                
                % 只有 3 个 BC 时，默认固定转动自由度，防止 Truss 计算奇异
                obj.BCode(4:6) = 1; 
                
                obj.XYZ(1:3) = data(5:7);
                
            elseif colCount == 10 || colCount == 11
                % 格式 B: 6 BCs Beam
                % Base: ID, BC1..BC6, X, Y, Z (10 cols)
                % Ext : ID, BC1..BC6, X, Y, Z, Temp (11 cols)
                
                obj.BCode(1:6) = round(data(2:7));
                obj.XYZ(1:3) = data(8:10);
                
                if colCount == 11
                    obj.Temperature = data(11);
                end
                
            else
                error('Unknown Node Data Format. Expected 7 or 10 columns, got %d', colCount);
            end
        end
        
        % Debug output
        function PrintInfo(obj)
            fprintf('Node %d: [%.2f, %.2f, %.2f] T=%.2f BC:[%d %d %d %d %d %d]\n', ...
                obj.ID, obj.XYZ(1), obj.XYZ(2), obj.XYZ(3), obj.Temperature, obj.BCode');
        end
    end
end