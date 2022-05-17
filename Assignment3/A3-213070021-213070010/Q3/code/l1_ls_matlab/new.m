clear;

rand('state',0);
randn('state',0);

angles = 10:10:180;
% disp(angles)

[image_1,cmap] = imread('slice_50.png');
padded_image = zeros(220,220);
padded_image(110 - round(size(image_1,1)/2)+1: 110 + round(size(image_1,1)/2)-1, 110 - round(size(image_1,2)/2)+1: 110 + round(size(image_1,2)/2)-1) = image_1;
% disp(size(padded_image))

% imshow(image_1, cmap)

x = ISTA(padded_image, angles);

image_reconstructed = idct2(reshape(x,[220,220]));
disp(["Size of reconstructed image : ",size(image_reconstructed)])

imshow(image_reconstructed,[])

% A = @sensing_matrix;
% disp(size(x))

function x = ISTA(image, angles)

    projections = radon(image,angles);
    y = reshape(projections, [],1);

    m = 5670;
    n = 48400;

    J1 = randperm(48400);
    J2 = randperm(5670);

    A = sensingMatrix(J2,J1,angles);
    At = A';


%     disp(["Is At an object ? ",isobject(At)])
%     disp(["Is At a vector ?", ~isvector(At)])

    lambda = find_lambdamax_l1_ls(At,y)/9;
    disp(lambda)
    ret_tol = 0.001;
    
    disp([m,n])
    [x,status,history] = l1_ls(A,At,m,n,y,lambda,ret_tol);
%     disp(history)
%     [x,status] = l1_ls(A,y,lambda, ret_tol);
    disp(status)
end




