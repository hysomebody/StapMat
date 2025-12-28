classdef StaticSolver < Solver
    % STATICSOLVER Linear Static Solver
    %
    % Purpose:
    %   Solve Ku = F for multiple load cases.
    
    methods
        function Solve(obj, domainObj)
            fprintf(' [StaticSolver] Start Solving...\n');
            
            % 1. Check Stiffness Matrix
            if isempty(domainObj.GlobalK)
                error('Global Stiffness Matrix K is empty! Run Assemble() first.');
            end
            
            K = domainObj.GlobalK;
            nCases = domainObj.NLCASE;
            
            % 2. Loop over all Load Cases
            for lc = 1:nCases
                % Call Domain to assemble force vector F for this case
                F = domainObj.AssembleForce(lc);
                
                % Solve Linear System
                fprintf('   Processing Load Case %d...\n', lc);
                U = K \ F;
                
                % Update Domain with results (Optional: print or store)
                domainObj.UpdateNodalDisplacements(U);
            end
            
            fprintf(' [StaticSolver] Solution Finished.\n');
        end
    end
end