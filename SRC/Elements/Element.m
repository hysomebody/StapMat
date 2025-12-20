% ELEMENT Abstract Base Class for all Finite Elements
%
% Purpose:
%   Define the common interface for all element types (Truss, Beam, etc.).
%   Store topology (Nodes) and property (Material) references.
% 功能：单元抽象基类。定义了所有单元（杆、梁）必须具备的通用属性（如节点连接、材料引用）
% 和抽象接口（如刚度矩阵计算）。
% 
% Call procedures:
%   ./Node.m
%   ./Material.m 
%   
%
% Called by:
%   ./TrussElement.m (via inheritance)
%   ./Domain.m

classdef (Abstract) Element < handle
    properties
        ID          % Element Number
        Nodes       % (Array of Node objects) Handle to connected nodes
        Material    % (Material object) Handle to assigned material
    end
    
    methods
        % Constructor
        function obj = Element()
            obj.ID = 0;
            % Nodes and Material will be assigned during Read
        end
        
        % Initialize nodes (Helper to link Node handles)
        function SetNodes(obj, globalNodeList, nodeIDs)
            % globalNodeList: Domain.NodeList (all nodes)
            % nodeIDs: Array of node IDs for this element
            
            numNodes = length(nodeIDs);
            obj.Nodes = Node.empty(0, numNodes); % Pre-allocate handle array
            
            for i = 1:numNodes
                % Note: Assuming nodeIDs match index for simplicity. 
                % In complex cases, we might search or use a map.
                obj.Nodes(i) = globalNodeList(nodeIDs(i));
            end
        end
    end
    
    methods (Abstract)
        % Abstract method: Calculate Element Stiffness Matrix
        % Must be implemented by subclasses (BarElement, BeamElement)
        % Returns: k (Matrix), usually 12x12 for 3D generic case
        k = CalcStiffness(obj)
        
        % Abstract method: Read element data from file
        Read(obj, fid, expectedID, materialList, nodeList)

        % Abstract method: Calculate Element Stress/Forces
        % Input: globalDisplacement (全场位移向量)
        % Output: stressResults (结构体或数组，视单元类型而定)
        results = CalcStress(obj, globalDisplacement)

        % Calculate Element Mass Matrix (Consistent or Lumped)
        m = CalcMass(obj)
    end
end