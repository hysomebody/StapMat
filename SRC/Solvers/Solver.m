classdef (Abstract) Solver < handle
    % SOLVER Abstract base class for all solvers
    %
    % Purpose:
    %   Define the standard interface that all solvers (Static, Dynamic) must follow.
    
    methods (Abstract)
        % Abstract method: Execute the solution process
        % Input: domainObj - The Domain instance containing K, M, F, etc.
        Solve(obj, domainObj)
    end
end