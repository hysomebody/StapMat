% 材料基类
%
% 功能说明：
%   所有材料类的父类 (Base Class)。
%   1. 存储通用的材料属性（如 ID、杨氏模量 E、密度 Density）。
%   2. 定义标准接口（如 Read 方法），供子类继承和重写。
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
        ID  % 材料编号
        E   % Young's Modulus
        Density % Mass Density,
    end
    
    methods
        function obj = Material()
            obj.ID = 0;
            obj.E = 0;
            obj.Density = 0;
        end
        
        % 读取材料数据 (虚方法)
        % 实际的读取逻辑由子类实现，不同单元用到的材料需要不同的参数
        function Read(obj, fid, expectedID)
        end
    end
end