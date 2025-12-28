classdef StaticSolver < Solver
    % 求解静力学问题
    %
    % 仅求解 K * u = F
    
    methods
        function Solve(obj, domainObj)
            fprintf(' [StaticSolver] Start Solving...\n');
            
            % 先检查一些刚度矩阵
            if isempty(domainObj.GlobalK)
                error('Global Stiffness Matrix K is empty! Run Assemble() first.');
            end
            
            K = domainObj.GlobalK;
            nCases = domainObj.NLCASE;
            
            for lc = 1:nCases
                F = domainObj.AssembleForce(lc);
                fprintf('   Processing Load Case %d...\n', lc);
                U = K \ F;
                
                domainObj.UpdateNodalDisplacements(U);
            end
            
            fprintf(' [StaticSolver] Solution Finished.\n');
        end
    end
end