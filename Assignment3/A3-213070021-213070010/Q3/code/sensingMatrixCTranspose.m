classdef sensingMatrixCTranspose
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties
        projection_size
        data_size
        projection_angles_1
        projection_angles_2
    end
    
    methods
        function obj = sensingMatrixCTranspose(projection_size, data_size, projection_angles_1, projection_angles_2)
             %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here

            obj.projection_size = projection_size;
            obj.data_size = data_size;
            obj.projection_angles_1 = projection_angles_1;
            obj.projection_angles_2 = projection_angles_2;
        end
        
        function output = mtimes(At, vect)
             %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here

            len_vec = length(vect);
            y1 = vect(1:0.5*len_vec);
            y2 = vect(0.5*len_vec+1:end);
            
            y1 = reshape(y1, At.projection_size, size(At.projection_angles_1,2));
            y2 = reshape(y2, At.projection_size, size(At.projection_angles_2,2));
            
            beta1 = iradon(y1, At.projection_angles_1, 'linear','Ram-Lak', 1, At.data_size);
            delta_beta1 = iradon(y2, At.projection_angles_2, 'linear','Ram-Lak', 1, At.data_size);
            
            x1 = dct2(beta1);
            delta_x1 = dct2(delta_beta1);
            
            output = [x1(:) + delta_x1(:); delta_x1(:)];
        end
    end
end
            
        