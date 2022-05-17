%%
addpath('l1_ls_matlab');

%%
clear;
clc;
%% a)

%%%%%% NOTE : Keep Running this Block until we get solved status for all lambdas 
running_status = true;
r_iter = 0;
while running_status
    % Generating X and Y from a Uniform random variable
    n = 500;
    m = 200;
    x_0 = 18;
    lambda_set = [0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1, 2, 5, 10, 15, 20, 30, 50, 100];
    
    % Generating the sensing matrix
    p=0.5;
    A=(rand(200,500)<p);
    A = 2*A/sqrt(m) - 1/sqrt(m);
    A_t = A';
    
    % Selecting the reconstruction and valdation set sizes
    R = 0.9*m;
    V = 0.1*m;
    
    %  Generating the values of x from a Uniform Distritbution of
    %  range(0,1000), and computing Y correspondingly
    x = zeros(n,1);
    indices = randperm(n, x_0);
    x(indices) = randi([0,1000],x_0,1);
    y = A*x;
    sigma = 0.05*mean(abs(y));
    y = y + (sigma.^(2)).*randn(m,1);
    
    % Selecting R random indices from m to create the reconstruction set
    arr_indices = 1:m;
    R_indices = randperm(m,0.9*m);
    R_y = y(R_indices);
    
    % Considering the remaining as validation set
    V_indices = setdiff(arr_indices, R_indices);
    V_y = y(V_indices);
    
    % Considering corresponding rows of A for the reconstruction and
    % validations set
    A_R = A(R_indices,:);
    A_V = A(V_indices,:);
    
    
    % Using the l1_ls solver, computing the estimated x from the
    % reconstruction set for different values of lambda
    disp(["Running iteration : ", r_iter])
    validation_array = zeros(17,1);
    rmse_array = zeros(17,1);
    failed_set = 0;
    for i=1:length(lambda_set)
        disp(["Estimating X for lambda : ", lambda_set(:,i)])
        quiet = true;
        [x_estimated, status] = l1_ls(A_R, R_y, lambda_set(:,i), 0.01, quiet);
        disp(["Execution Status : ", status])
        if status == "Failed"
            failed_set = failed_set + 1;
        end
        validation_error = (V_y - A_V*x_estimated)'*(V_y - A_V*x_estimated)/length(V_y);
        disp(["Validation Error : ",validation_error])
        validation_array(i) = validation_error;
        difference = x_estimated - x;
        rmse = sqrt(difference'*difference)/sqrt(x'*x);
        disp(["RMSE Error : ", rmse])
        rmse_array(i) = rmse;
        disp("---------------------------------------------------------------")
    end
    if failed_set == 0
        running_status = false;
    end

   r_iter = r_iter + 1;
end

%%
disp(["Number of running iterations required to solve l1-ls for all lambdas : ", r_iter])

% Plotting the figure for lambda vs validation scores
figure(1)
plot(log(lambda_set), validation_array)
xticks(log(lambda_set))
xlabel('log(\lambda)')
ylabel('Val error')
title('Validation Error vs log(\lambda)')


% Plotting the figure for lambda vs RMSE scores
figure(2)
plot(log(lambda_set), rmse_array)
xticks(log(lambda_set))
xlabel('log(\lambda)')
ylabel('RMSE')
title('RMSE vs log(\lambda)')

% Optimal value of lambda obtained from both the methods
[~,I1] = min(validation_array);
[~,I2] = min(rmse_array);
disp(["Best value of log lambda for minimum cross validation is", log(lambda_set(:,I1)), "and the corres. lambda value is", lambda_set(:,I1)])
disp(["Best value of log lambda for minimum rmse is", log(lambda_set(:,I2)), "and the corres. lambda value is", lambda_set(:,I2)])




%% b) When V and R are not disjoint but coincident sets

%%%%%% NOTE : Keep Running this Block until we get solved status for all lambdas 
running_status = true;
r_iter = 0;
while running_status
    arr_indices = 1:m;
    R_indices = randperm(m,0.9*m);
    R_y = y(R_indices);
    
    V_indices = randperm(m,0.1*m);
    V_y = y(V_indices);
    
    % Considering corresponding rows of A for the reconstruction and
    % validations set
    A_R = A(R_indices,:);
    A_V = A(V_indices,:);

    % Using the l1_ls solver, computing the estimated x from the
    % reconstruction set for different values of lambda
    disp(["Running iteration : ", r_iter])
    validation_array = zeros(17,1);
    rmse_array = zeros(17,1);
    failed_set = 0;
    for i=1:length(lambda_set)
        disp(["Estimating X for lambda : ", lambda_set(:,i)])
        quiet = true;
        [x_estimated, status] = l1_ls(A_R, R_y, lambda_set(:,i), 0.01, quiet);
        disp(["Execution Status : ", status])
        if status == "Failed"
            failed_set = failed_set + 1;
        end
        validation_error = (V_y - A_V*x_estimated)'*(V_y - A_V*x_estimated)/length(V_y);
        disp(["Validation Error : ",validation_error])
        validation_array(i) = validation_error;
        difference = x_estimated - x;
        rmse = sqrt(difference'*difference)/sqrt(x'*x);
        disp(["RMSE Error : ", rmse])
        rmse_array(i) = rmse;
        disp("---------------------------------------------------------------")
    end
    if failed_set == 0
        running_status = false;
    end

   r_iter = r_iter + 1;
end

disp(["Number of running iterations required to solve l1-ls for all lambdas : ", r_iter])
% Plotting the figure for lambda vs validation scores
figure(3)
plot(log(lambda_set), validation_array)
xticks(log(lambda_set))
xlabel('log(\lambda)')
ylabel('Val error')
title('Validation Error vs log(\lambda)')

% Plotting the figure for lambda vs RMSE scores
figure(4)
plot(log(lambda_set), rmse_array)
xticks(log(lambda_set))
xlabel('log(\lambda)')
ylabel('RMSE')
title('RMSE vs log(\lambda)')

% Optimal value of lambda obtained from both the methods
% Optimal value of lambda obtained from both the methods
[~,I1] = min(validation_array);
[~,I2] = min(rmse_array);
disp(["Optimal log lambda with coincident sets for minimum cross validation is", log(lambda_set(:,I1)), "and the corres. lambda value is", lambda_set(:,I1)])
disp(["Optimal log lambda with coincident sets for minimum rmse is", log(lambda_set(:,I2)), "and the corres. lambda value is", lambda_set(:,I2)])

%% 
m =200;
p=0.5;
A=(rand(200,500)<p);
A = 2*A/sqrt(m) - 1/sqrt(m);
A_t = A';