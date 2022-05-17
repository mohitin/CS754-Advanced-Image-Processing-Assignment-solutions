classdef sensingMatrix
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        projections_size % size of the projections
        data_size % dimension of the 2D DCT
        projection_angles % angles used for projections
    end
    
    methods
        function obj = sensingMatrix(m,n,angles)
            %sensingMatrix: Construct an instance of this class
            %   Default constructor

            obj.projections_size = m;
            obj.data_size = n;
            obj.projection_angles = angles;
        end
        
        function output = mtimes(A, obj)
            %mtimes: Overriding the matrix multiplication
            %   Computed radon(idct2(obj)); forward model

            image_dimension = uint64(sqrt(A.data_size));
            x = reshape(obj, image_dimension, []);
            beta_val = idct2(x);
            output = radon(beta_val, A.projection_angles);
            output = output(:);
        end
    end
end
            
            