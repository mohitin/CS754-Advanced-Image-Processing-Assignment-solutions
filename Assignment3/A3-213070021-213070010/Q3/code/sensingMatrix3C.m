classdef sensingMatrix3C
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
        function obj = sensingMatrix3C(projection_size, data_size, angles_1, angles_2, angles_3)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.projected_angles_1 = angles_1;
            obj.projected_angles_2 = angles_2;
            obj.projected_angles_3 = angles_3;
            obj.projection_size = projection_size;
            obj.data_size = data_size;
        end
        
        function output = mtimes(A,obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            beta_1 = obj(1:length(obj)/3);
            delta_beta_1 = obj((length(obj)/3)+1:2*length(obj)/3);
            delta_beta_2 = obj((2*length(obj)/3)+1:end);

            beta_1_dct2 = reshape(beta_1, A.data_size, A.data_size);
            delta_beta_1_dct2 = reshape(delta_beta_1, A.data_size, A.data_size);
            delta_beta_2_dct2 = reshape(delta_beta_2, A.data_size, A.data_size);
            
            im_beta_1 = idct2(beta_1_dct2);
            im_delta_beta_1 = idct2(delta_beta_1_dct2);
            im_delta_beta_2 = idct2(delta_beta_2_dct2);

            radon1_beta_1 = radon(im_beta_1, A.projected_angles_1);
            radon2_beta_1 = radon(im_beta_1, A.projected_angles_2);
            radon3_beta_1 = radon(im_beta_1, A.projected_angles_3);

            radon2_delta_beta_1 = radon(im_delta_beta_1, A.projected_angles_2);
            radon3_delta_beta_2 = radon(im_delta_beta_2, A.projected_angles_3);

            output = [radon1_beta_1(:); radon2_beta_1(:) + radon2_delta_beta_1(:); radon3_beta_1(:) + radon3_delta_beta_2(:)];
            
        end
    end
end

