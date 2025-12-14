%* *****************************************************************
%* - Basic data class of STAPMAT                                   *
%*                                                                 *
%* - Purpose:                                                      *
%*     Storing variables which control the running of STAPMAT      *
%*                                                                 *
%* - Programmed by:                                                *
%*     LeiYang Zhao, Yan Liu,                                      *
%*     Computational Dynamics Group, School of Aerospace           *
%*     Engineering, Tsinghua University, 2019.02.20                *
%*                                                                 *
%* *****************************************************************
classdef ControlData
    properties
        NUMNP;         % Total number of nodal points
                       % = 0 : Program stop

        NPAR;          % Element group control data
                       %   NPAR(1) - Element type
                       %             1 : Truss element
                       %   NPAR(2) - Number of elements
                       %   NPAR(3) - Number of different sets of material
                       %             and cross-sectional constants
        NUMEG;         % Total number of element groups, > 0
        NLCASE;        % Number of load case (>0)
        LL;            % Load case number
        NLOAD;         % The number of concentrated loads applied in this load case

        MODEX;         % Solution mode: 0 - data check only; 1 - execution

        TIM;           % Timing information
        HED;           % Master heading information for usr in labeling the output

        IIN;           % file pointer used for input
        IOUT;          % file pointer used for output
    end
end
