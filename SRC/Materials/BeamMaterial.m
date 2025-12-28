% 
% 功能：增加梁单元需要的剪切模量G,截面积A,惯性矩Iy Iz, 扭转常数J
% Call procedures:
%   ./Material.m - Material() 
%
% Called by:
%   ./Domain.m - ReadElements()

classdef BeamMaterial < Material
    properties
        G       
        Area    
        Iy      
        Iz      
        J       
    end
    
    methods
        function obj = BeamMaterial()
            obj@Material();
            obj.G = 0.0;
            obj.Area = 0.0;
            obj.Iy = 0.0;
            obj.Iz = 0.0;
            obj.J = 0.0;
        end
        
        function Read(obj, fid, expectedID)
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>
            
            if isempty(data)
                error('Error reading BeamMaterial data');
            end
            
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Material ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            obj.E  = data(2);
            obj.G  = data(3);
            obj.Area = data(4);
            obj.Iy = data(5);
            obj.Iz = data(6);
            obj.J  = data(7);
            % 动力学部分需要密度，如果输入文件中没定义默认为0
            if length(data) >= 8
                obj.Density = data(8);
            else
                obj.Density = 0.0;
                fprintf('Warning: No density found for BeamMaterial %d. Assuming 0.\n', obj.ID);
            end

        end
    end
end