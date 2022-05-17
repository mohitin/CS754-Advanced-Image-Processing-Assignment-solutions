classdef sensingMatrixTranspose
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        projections_size % size of the projections
        data_size % dimension of the 2D DCT
        projection_angles % angles used for projections
    end
    
    methods
        function obj = sensingMatrixTranspose(m,n,angles)
            %UNTITLED Construct an instance of this class
            %  Default constructor

            obj.projections_size = m;
            obj.data_size = n;
            obj.projection_angles = angles;
        end
        
        function output = mtimes(A, vect)
            %mtimes: Overriding the matrix multiplication in transpose
            %       computes iradon(dct2(vect)) ; backward model.

            projections_dim = double(A.projections_size/size(A.projection_angles,2));
            x = reshape(vect, projections_dim, []);
            projecs = iradon(x, A.projection_angles, 'linear', 'Ram-Lak', 1, sqrt(A.data_size));
            output = dct2(projecs);
            output = output(:);
        end
    end
end