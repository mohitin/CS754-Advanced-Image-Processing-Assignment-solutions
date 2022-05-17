classdef sensingMatrixC
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        projection_size % size of each projection
        data_size % dimension of the 2D DCT
        projection_angles_1 % projected angles 1
        projection_angles_2 % projected angles 2
    end
    
    methods
        function obj = sensingMatrixC(projection_size,data_size, projection_angles_1, projection_angles_2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.data_size = data_size;
            obj.projection_size = projection_size;
            obj.projection_angles_1 = projection_angles_1;
            obj.projection_angles_2 = projection_angles_2;
        end
        
        function output = mtimes(A, obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            x1 = obj(1:length(obj)/2);
            delta_x1 = obj(0.5*length(obj)+1:end);
            x1 = reshape(x1, A.data_size, A.data_size);
            delta_x1 = reshape(delta_x1, A.data_size, A.data_size);
            beta1 = idct2(x1);
            delta_beta1 = idct2(delta_x1);
            radon1_beta1 = radon(beta1, A.projection_angles_1);
            radon2_beta1 = radon(beta1, A.projection_angles_2);
            radon2_delta_beta = radon(delta_beta1, A.projection_angles_2);
            output = [radon1_beta1(:); radon2_beta1(:) + radon2_delta_beta(:)];
        end
    end
end