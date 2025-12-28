% 
% 功能：增加实体单元需要的泊松比 Nu、密度 Rho 和热膨胀系数 Alpha
%
% Call procedures:
%   ./Material.m - Material() (Base Constructor)
%
% Called by:
%   ./Domain.m - ReadElements()
% 程云志

classdef TetraMaterial < Material
    properties
        Nu     
        Rho    
        Alpha  
    end
    
    methods

        function obj = TetraMaterial()
            obj@Material(); 
            obj.Nu = 0.3;  % 默认钢材泊松比
            obj.Rho = 0.0;
            obj.Alpha = 0.0;
        end
        
        function Read(obj, fid, expectedID)
            
            lineStr = fgetl(fid);
            data = str2num(lineStr); 
            
            if isempty(data)
                error('Error reading TetraMaterial data');
            end
            
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Material ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            
            % 顺序: ID, RHO, E, NU, ALPHA
            obj.Rho = data(2);
            obj.E   = data(3); 
            obj.Nu  = data(4);
            
            if length(data) >= 5
                obj.Alpha = data(5);
            else
                obj.Alpha = 0.0;
            end
        end
    end
end