%
% 功能：除Material类中共有的材料属性外，增加桁架单元需要的截面积 Area。
% Call procedures:
%   ./Material.m - Material() 
%
% Called by:
%   ./Domain.m 
    
classdef TrussMaterial < Material
    properties
        Area 
        Alpha % Thermal Expansion Coefficient
    end
    
    methods
        function obj = TrussMaterial()
            obj@Material(); 
            obj.Area = 0.0;
        end
        
        % Format: ID E Area Density
        function Read(obj, fid, expectedID)
            
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>
            
            if isempty(data)
                error('Error reading TrussMaterial data');
            end
            
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Material ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;          
            obj.E = data(2);    
            obj.Area = data(3); 

            if length(data) >= 4
                obj.Density = data(4);
            else
                obj.Density = 0.0;
                fprintf('Warning: No density found for TrussMaterial %d. Assuming 0.\n', obj.ID);
            end
            if length(data) >= 5
                obj.Alpha = data(5);
            else
                obj.Alpha = 0.0;
            end

        end
    end
end