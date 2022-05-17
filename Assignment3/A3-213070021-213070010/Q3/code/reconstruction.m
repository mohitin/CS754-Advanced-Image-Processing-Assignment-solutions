clear;
close all;
clc;
%% 
% Creating projections at 18 uniformly spaced angles

% Equispaced Angles
angles = 10:10:180;

% Reading the image slices
image_slice_50 = im2double(imread('../slices/slice_50.png'));
image_slice_51 = im2double(imread('../slices/slice_51.png'));

% Padding the image to size 225x225
padded_slice_50 = padarray(image_slice_50, [22,4], 0,'both');
padded_slice_51 = padarray(image_slice_51, [22,4], 0,'both');

% Obtaining the projections using the Radon transform
projections_slice_50 = radon(padded_slice_50, angles);
projections_slice_51 = radon(padded_slice_51, angles);

%% a
% Tomographic reconstruction using filtered back projection
% FBP using the iradon matlab function with Ram-Lak filtering
[reconstructed_slice_50,H] = iradon(projections_slice_50, angles,'linear','Ram-Lak',1,225);
[reconstructed_slice_51,H] = iradon(projections_slice_51, angles,'linear','Ram-Lak',1,225);

% Displaying the reconstructed images
figure(1)
imshow(reconstructed_slice_50,[])
title('Reconstruction of Slice 50 using filtered back projection')
imwrite(reconstructed_slice_50,'../saved_outputs/a_recon_slice_50.jpg')


figure(2)
imshow(reconstructed_slice_51,[])
title('Reconstruction of Slice 51 using filtered back projection')
imwrite(reconstructed_slice_51,'../saved_outputs/a_recon_slice_51.jpg')


%% b 
% Independent CS based reconstruction of Slices 50 and 51

%Adding the l1_ls library to path
addpath('l1_ls_matlab')

% Vectorizing the projections of slice 50 and 51
y_slice_50 = projections_slice_50(:);
y_slice_51 = projections_slice_51(:);

% Sizes of input DCT (n), and output projection (m)
n = size(padded_slice_50(:),1);
m = size(y_slice_50,1);

% Initiating the sensing matrix radon*idct2 and its transpose iradon*ct2
A = sensingMatrix(m,n,angles);
At = sensingMatrixTranspose(m,n,angles);

% Parameters for the L1 solver
lambda = 0.01;
rel_tol = 1e-7;
quiet = true;

% Obtaining the vectorized DCT from the L1 solver for slices 50 and 51
[csbased_dct_50, status] = l1_ls(A, At, m, n, y_slice_50, lambda, rel_tol, quiet);
disp(["Part b) Status of l1 solver for Slice 50 : ", status])
[csbased_dct_51, status] = l1_ls(A, At, m, n, y_slice_51, lambda, rel_tol, quiet);
disp(["Part b) Status of l1 solver for Slice 51 : ", status])

% Obtaining the reconstructed image from the vectorized DCT
reconstructed_slice_50_idcs = idct2(reshape(csbased_dct_50, 225, 225));
reconstructed_slice_51_idcs = idct2(reshape(csbased_dct_51, 225, 225));

% Displaying the images
figure(3)
imshow(reconstructed_slice_50_idcs, [])
title('CS based reconstruction of Slice 50')
imwrite(reconstructed_slice_50_idcs,'../saved_outputs/b_recon_idcs_slice_50.jpg')

figure(4)
imshow(reconstructed_slice_51_idcs, [])
title('CS based reconstruction of Slice 51')
imwrite(reconstructed_slice_51_idcs,'../saved_outputs/b_recon_idcs_slice_51.jpg')

%% c 
% Coupled CS Reconstruction using two slices : Slices 50 and 51

%Adding the l1_ls library to path
addpath('l1_ls_matlab')

% Different sets of equispaced angles
angles_1 = 10:10:180;
angles_2 = 7:10:177;

% Obtaining projections for slice 50 using set1 and slice 52 using set 2
projections_slice_50 = radon(padded_slice_50, angles_1);
projections_slice_51 = radon(padded_slice_51, angles_2);

% Vectorizing the projections of slice 50 and 51
y_slice_50 = projections_slice_50(:);
y_slice_51 = projections_slice_51(:);

% Concatenating (row-wise) the vectorized projections
y_50_51 = [y_slice_50; y_slice_51];

% Size of the output concatenated projection (m), and the concatenated
% vectorized DCT (n)
m = size(y_50_51, 1);
n = size(padded_slice_50(:), 1) + size(padded_slice_51(:), 1);

% Size of a single radon transform
projection_size = size(projections_slice_50, 1);

% Initiating the sensing matrix radon*idct2 and its transpose iradon*ct2
A = sensingMatrixC(projection_size,225,angles_1,angles_2);
At = sensingMatrixCTranspose(projection_size,225,angles_1,angles_2);

% Parameters for the L1 solver
lambda = 0.001;
rel_tol = 1e-6;
quiet = true;

% Obtaining the vectorized coupled DCT from the L1 solver for slices 50 and
% 51
[coupled_dct, status] = l1_ls(A, At, m, n, y_50_51, lambda, rel_tol, quiet);
disp(["Part c) Status of l1 solver for Coupled Reconstruction of Slice 50 and 51: ", status])

% Obtaining individual vectorized DCT of slices 50 and 51
beta1 = coupled_dct(1:(n/2));
delta_beta1 = coupled_dct(1+(n/2):end);
beta2 = beta1 + delta_beta1;

% Obtaining the reconstructed image from the vectorized DCT
reconstructed_slice_50_coupled = idct2(reshape(beta1,225,225));
reconstructed_slice_51_coupled = idct2(reshape(beta2,225,225));

% Displaying the images
figure(5)
imshow(reconstructed_slice_50_coupled)
title('Coupled Reconstruction of Slice 50')
imwrite(reconstructed_slice_50_coupled,'../saved_outputs/c_recon_slice_coupled_50.jpg')

figure(6)
imshow(reconstructed_slice_51_coupled)
title('Coupled Reconstruction of Slice 51')
imwrite(reconstructed_slice_51_coupled,'../saved_outputs/c_recon_slice_coupled_51.jpg')

%% d 
% Coupled CS reconstruction using three slices : Slices 50,51 and 52

% Different sets of equispaced angles
angles_1 = 10:10:180;
angles_2 = 7:10:177;
angles_3 = 4:10:174;

% Reading Slice 52
image_slice_52 = im2double(imread('../slices/slice_52.png'));

% Padding Slice 52 to 225x225
padded_slice_52 = padarray(image_slice_52,[22,4],0,'both');

% Obtaining projections for slice 50 using set1, slice 52 using set 2 and
% slice 53 using set 3
projections_slice_50 = radon(padded_slice_50, angles_1);
projections_slice_51 = radon(padded_slice_51, angles_2);
projections_slice_52 = radon(padded_slice_52, angles_3);

% Vectorizing the projections of slice 50, 51 and 52
y_50 = projections_slice_50(:);
y_51 = projections_slice_51(:);
y_52 = projections_slice_52(:);

% Concatenating (row-wise) the vectorized projections
y_50_51_52 = [y_50; y_51; y_52];

% Size of the output concatenated projection (m), and the concatenated
% vectorized DCT (n)
m = size(y_50_51_52,1);
n = size(padded_slice_50(:),1) + size(padded_slice_51(:),1) + size(padded_slice_52(:),1);

% Size of a single radon transform
projection_size = size(projections_slice_50,1);

% Initiating the sensing matrix radon*idct2 and its transpose iradon*ct2
A = sensingMatrix3C(projection_size,225,angles_1,angles_2,angles_3);
At = sensingMatrix3CTranspose(projection_size,225,angles_1,angles_2,angles_3);

% Parameters for the L1 solver
lambda = 0.001;
rel_tol = 1e-6;
quiet = true;

% Obtaining the vectorized coupled DCT from the L1 solver for slices 50, 51
% and 52
[coupled_dct_3, status] = l1_ls(A, At, m, n, y_50_51_52, lambda, rel_tol, quiet);
disp(["Status of l1 solver for Coupled Reconstruction of Slices 50, 51 and 52 : ", status])

% Obtaining individual vectorized DCT of slices 50, 51 and 52
beta1 = coupled_dct_3(1:n/3);
delta_beta1 = coupled_dct_3(1+n/3:2*n/3);
delta_beta2 = coupled_dct_3(1+2*n/3:end);
beta2 = beta1 + delta_beta1;
beta3 = beta1 + delta_beta1 + delta_beta2;

% Obtaining the reconstructed image from the vectorized DCT
reconstructed_slice_50_coupled3 = idct2(reshape(beta1,225,225));
reconstructed_slice_51_coupled3 = idct2(reshape(beta2,225,225));
reconstructed_slice_52_coupled3 = idct2(reshape(beta3,225,225));

% Displaying the images
figure(7)
imshow(reconstructed_slice_50_coupled3)
title('Coupled Reconstruction of Slice 50 using threes slices')
imwrite(reconstructed_slice_50_coupled3,'../saved_outputs/d_recon_slice_coup3_50.jpg')

figure(8)
imshow(reconstructed_slice_51_coupled3)
title('Coupled Reconstruction of Slice 51 using three slices')
imwrite(reconstructed_slice_51_coupled3,'../saved_outputs/d_recon_slice_coup3_51.jpg')

figure(9)
imshow(reconstructed_slice_52_coupled3)
title('Coupled Reconstruction of Slice 52 using three slices')
imwrite(reconstructed_slice_52_coupled3,'../saved_outputs/d_recon_slice_coup3_52.jpg')









