% MATERIAL Base class for material properties
%
% Purpose:
%   Store basic material properties (like ID, Young's Modulus E).
%   Serve as a parent for TrussMaterial and BeamMaterial.
% 功能：材料基类。为了支持多种单元，将 Read 定义为抽象方法，或者作为通用接口。
%
% Call procedures:
%   None
%
% Called by:
%   ./TrussMaterial.m 
%   ./TrussMaterial.m 
%   ./Element.m 
    
classdef Material < handle
    properties
        ID  % Material Set Number (nset)
        E   % Young's Modulus
    end
    
    methods
        % Constructor
        function obj = Material()
            obj.ID = 0;
            obj.E = 0;
        end
        
        % Read material data (Virtual/Base implementation)
        % Subclasses should override this to read extra properties (Area, I, etc.)
        function Read(obj, fid, expectedID)
            % called by: Domain (or specific ElementGroup reading logic)
            
            % This is a placeholder base Read. 
            % In practice, subclasses will handle the specific format.
            % Example format: ID E [Other Props...]
            
            % Note: Actual reading logic usually happens in the subclass 
            % because the number of parameters varies.
        end
    end
end