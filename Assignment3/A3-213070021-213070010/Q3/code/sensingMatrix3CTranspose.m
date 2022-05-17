classdef sensingMatrix3CTranspose
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        projected_angles_1
        projected_angles_2
        projected_angles_3
        projection_size
        data_size
    end
    
    methods
        function obj = sensingMatrix3CTranspose(projection_size, data_size, angles_1, angles_2, angles_3)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.projected_angles_1 = angles_1;
            obj.projected_angles_2 = angles_2;
            obj.projected_angles_3 = angles_3;
            obj.projection_size = projection_size;
            obj.data_size = data_size;
            
        end
        
        function output = mtimes(At,vect)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            y1 = vect(1:length(vect)/3);
            y2 = vect((length(vect)/3)+1:2*length(vect)/3);
            y3 = vect((2*length(vect)/3)+1:end);

            y1 = reshape(y1,At.projection_size, size(At.projected_angles_1,2));
            y2 = reshape(y2,At.projection_size, size(At.projected_angles_1,2));
            y3 = reshape(y3,At.projection_size, size(At.projected_angles_1,2));

            im_beta_1 = iradon(y1, At.projected_angles_1, 'linear', 'Ram-Lak', 1, At.data_size);
            im_delta_beta_1 = iradon(y2, At.projected_angles_2, 'linear', 'Ram-Lak', 1, At.data_size);
            im_delta_beta_2 = iradon(y3, At.projected_angles_3, 'linear', 'Ram-Lak', 1, At.data_size);

            beta_1 = dct2(im_beta_1);
            delta_beta_1 = dct2(im_delta_beta_1);
            delta_beta_2 = dct2(im_delta_beta_2);

            output = [beta_1(:) + delta_beta_1(:) + delta_beta_2(:); delta_beta_1(:); delta_beta_2(:)];
        end
    end
end

