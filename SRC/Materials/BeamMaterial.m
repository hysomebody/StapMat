% BEAMMATERIAL Material definition for 3D Beam elements
%
% Purpose:
%   Extends base Material to include Sectional properties for Beam.
%   Properties: E, G, Area, Iy, Iz, J
% 功能：增加梁单元需要的剪切模量G,截面积A,惯性矩Iy Iz, 扭转常数J
% Call procedures:
%   ./Material.m - Material() (Base Constructor)
%
% Called by:
%   ./Domain.m - ReadElements()

classdef BeamMaterial < Material
    properties
        G       % Shear Modulus
        Area    % Cross-sectional Area
        Iy      % Moment of inertia about local y-axis
        Iz      % Moment of inertia about local z-axis
        J       % Torsional constant (Polar moment of inertia)
        Alpha
    end
    
    methods
        % Constructor
        function obj = BeamMaterial()
            obj@Material(); % Call base constructor
            obj.G = 0.0;
            obj.Area = 0.0;
            obj.Iy = 0.0;
            obj.Iz = 0.0;
            obj.J = 0.0;
        end
        
        % Read material data from file stream
        % Format: ID(int) E(double) G(double) A(double) Iy(double)
        % Iz(double) J(double) Density(double)
        function Read(obj, fid, expectedID)
            % Call procedures: None
            % Called by: ./Domain.m - ReadElements()
            
            lineStr = fgetl(fid);
            data = str2num(lineStr); %#ok<ST2NM>
            
            if isempty(data)
                error('Error reading BeamMaterial data');
            end
            
            % 1. Check ID
            inputID = round(data(1));
            if inputID ~= expectedID
                error('Material ID mismatch. Expected: %d, Read: %d', expectedID, inputID);
            end
            obj.ID = inputID;
            
            % 2. Read Properties (Reference C++ BeamMaterial::Read line 924)
            obj.E  = data(2);
            obj.G  = data(3);
            obj.Area = data(4);
            obj.Iy = data(5);
            obj.Iz = data(6);
            obj.J  = data(7);
            if length(data) >= 8
                obj.Density = data(8);
            else
                obj.Density = 0.0;
                fprintf('Warning: No density found for BeamMaterial %d. Assuming 0.\n', obj.ID);
            end
            if length(data) >= 9
                obj.Alpha = data(9);
            else
                obj.Alpha = 0.0; 
            end

        end
    end
end