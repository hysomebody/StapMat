% TETRAMATERIAL Material definition for 3D Solid (Tetra) elements
%
% Purpose:
%   Extends base Material to include Poisson's ratio and Density.
%   Properties: E, Nu, Rho, Alpha
% 功能：增加实体单元需要的泊松比 Nu、密度 Rho 和热膨胀系数 Alpha
%
% Call procedures:
%   ./Material.m - Material() (Base Constructor)
%
% Called by:
%   ./Domain.m - ReadElements()

classdef TetraMaterial < Material
    properties
        Nu      % Poisson's Ratio (泊松比)
        Rho     % Density (密度)
        Alpha   % Thermal Expansion Coefficient (热膨胀系数)
    end
    
    methods
        % Constructor
        function obj = TetraMaterial()
            obj@Material(); % Call base constructor
            obj.Nu = 0.3;
            obj.Rho = 0.0;
            obj.Alpha = 0.0;
        end
        
        % Read material data from file stream
        % Format based on old ReadTe4.m: 
        % ID(int) Rho(double) E(double) Nu(double) Alpha(double)
        function Read(obj, fid, expectedID)
            
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>
            
            if isempty(data)
                error('Error reading TetraMaterial data');
            end
            
            % 1. Check ID
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Material ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            
            % 2. Read Properties
            % 顺序: ID, RHO, E, NU, ALPHA
            obj.Rho = data(2);
            obj.E   = data(3); % Inherited from Material
            obj.Nu  = data(4);
            
            if length(data) >= 5
                obj.Alpha = data(5);
            else
                obj.Alpha = 0.0;
            end
        end
    end
end