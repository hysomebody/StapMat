classdef (Abstract) Solver < handle
% 定义求解器的标准接口
    
% 抽象接口
    methods (Abstract)
        Solve(obj, domainObj)
    end
end