% 
% 功能：
%   定义所有单元类型（杆、梁、实体）通用的属性和接口。
%   1. 存储单元连接结构和材料属性
%   2. 实现刚度矩阵等的计算
% 
% Call procedures:
%   ./Node.m
%   ./Material.m 
%   
% Called by:
%   ./TrussElement.m (via inheritance)
%   ./Domain.m

classdef (Abstract) Element < handle
    properties
        ID          % 单元编号
        Nodes       % 节点对象数组
        Material    % 材料对象数组
    end
    
    methods
        function obj = Element()
            obj.ID = 0;
        end

        % 初始化节点连接关系
        function SetNodes(obj, globalNodeList, nodeIDs)
            numNodes = length(nodeIDs);
            obj.Nodes = Node.empty(0, numNodes); 
            
            for i = 1:numNodes
                obj.Nodes(i) = globalNodeList(nodeIDs(i));
            end
        end
    end
    
    methods (Abstract)
        k = CalcStiffness(obj)
        
        % 读取单元数据
        Read(obj, fid, expectedID, materialList, nodeList)

        % 计算单元应力
        results = CalcStress(obj, globalDisplacement)

        % 计算单元质量矩阵
        m = CalcMass(obj)
    end
end