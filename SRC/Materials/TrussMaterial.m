% TRUSSMATERIAL Material definition for Truss elements
%
% Purpose:
%   Extends base Material to include Cross-Sectional Area.
% 功能：除Material类中共有的材料属性外，增加桁架单元需要的截面积 Area。
% Call procedures:
%   ./Material.m - Material() (Base Constructor)
%
% Called by:
%   ./Domain.m (when parsing Truss groups)
    
classdef TrussMaterial < Material
    properties
        Area % Cross-sectional Area
    end
    
    methods
        % Constructor
        function obj = TrussMaterial()
            obj@Material(); % Call base constructor
            obj.Area = 0.0;
        end
        
        % Read material data from file stream
        % Format: ID(int) E(double) Area(double)
        function Read(obj, fid, expectedID)
            % called by: Domain.ReadTrussElementData (to be implemented)
            
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>
            
            if isempty(data)
                error('Error reading TrussMaterial data');
            end
            
            % 1. Check ID
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Material ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            
            % 2. Read Properties
            obj.E = data(2);    % From Base Class
            obj.Area = data(3); % Specific to Truss
        end
    end
end