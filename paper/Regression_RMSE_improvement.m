%% this is the code for experiments in Tree Representations of Brain Structural Connectivity via Persistent Homology by D. Li, P. Nguyen, Z. Zhang and D. Dunson
% Creater: Didong Li, Jan. 2021
% We focus on 45 cognitive traits in HCP data
% Apply linear regression and GP regression to adjacency matrix and tree
% representation, and calculate the root of mean squared error (RMSE) on
% test data
% baseline method: use the mean of training samples for prediction
% Comare the improvement of two methods over basedline method on two data representation
% Note: we also checked other regression methods including decision tree,
% boosting, SVM, etc. GP regression and Linear regression turn out to be the
% most representative
% All repression methods are based on MATLAB regression toolbox.


close all
clear

load('data/HCP_cortical_TensorData_desikan.mat')
load('data/HCP_Covariates.mat')
load('data/Trees_vector.mat')
load('data/Adjacency_vector.mat')
load('data/ind_ctns')

[A,idx1] = sort(allsubject_id);
[~,idx2] = histc(all_id,A);

[~,Adjacency_all] = pca(Adjacency_vector); 
Adjacency_vector = Adjacency_all(:,1:23);

N = size(ind_ctns,2);
REP = 10;

RMSE_Baseline = zeros(REP,N);


RMSE_Adj = zeros(1,N);
RMSE_Tree = zeros(1,N);


RMSE_GPR_Adj = zeros(REP,N);
RMSE_GPR_Tree = zeros(REP,N);


RMSE_LR_Adj = zeros(REP,N);
RMSE_LR_Tree = zeros(REP,N);


for i = 1: N
    Y = cog_measure(ind_ctns(i),idx2).';
    ind_null = find(isnan(Y));
    Y(ind_null)=[];
    
    
    [~,D] = size(Adjacency_vector);
    [~,d] = size(Trees_vector);
    n = size(Y,1);
    n_train = floor(0.8*n);
    n_test = n-n_train;
    
    X_Adj = Adjacency_vector;
    X_Adj(ind_null,:) = [];
    X_Tree = Trees_vector;
    X_Tree(ind_null,:) = [];
    
    
    display(['Checking the ', num2str(i),'-th trait......'])
    
    
    % Regression for randomely permuted data, repeat REP times
    
    for rep = 1:REP
        
        ind = randperm(n);
        X_test_Tree = X_Tree(ind(n_train+1:n),:);
        X_test_Adj = X_Adj(ind(n_train+1:n),:);
        Y_test = Y(ind(n_train+1:n));
        
        
        %% randomely permute the data, use 80% as training and the rest for
        % testing
        
        X_train_Adj = X_Adj(ind(1:n_train),:);
        
        X_train_Tree = X_Tree(ind(1:n_train),:);
        
        Y_train = Y(ind(1:n_train));
        
        
        
        train_Adj = [X_train_Adj, Y_train];
        train_Tree = [X_train_Tree, Y_train];
        test_Adj = [X_test_Adj, Y_test];
        test_Tree = [X_test_Tree, Y_test];
        
        
    
        %% Baseline: use mean of training response as prediction
        y_pred = mean(Y_train)*ones(n_test,1);
        RMSE_Baseline(rep,i) = (Y_test-y_pred).'*(Y_test-y_pred)/n_test;
        display(['Baseline: RMSE = ',num2str(RMSE_Baseline(rep,i)),', # trait = ',num2str(ind_ctns(i)),', rep = ', num2str(rep)])
        
        
        %% Linear Regression
        %Adjacency matrix representation
        mdl_Adj = fitlm(X_train_Adj,Y_train);
        y_pred = predict(mdl_Adj,X_test_Adj);
        RMSE_LR_Adj(rep,i) = (Y_test-y_pred).'*(Y_test-y_pred)/n_test;
        display(['Linear Regression_Adjacency: RMSE = ',num2str(RMSE_LR_Adj(rep,i)),', # trait = ',num2str(ind_ctns(i)),', rep = ', num2str(rep)])
        
        
        % Tree representation
        mdl_Tree = fitlm(X_train_Tree,Y_train);
        y_pred = predict(mdl_Tree,X_test_Tree);
        RMSE_LR_Tree(rep,i) = (Y_test-y_pred).'*(Y_test-y_pred)/n_test;
        display(['Linear Regression_Tree: RMSE = ',num2str(RMSE_LR_Tree(rep,i)),', # trait ',num2str(ind_ctns(i)),', rep = ', num2str(rep)])
        
        %% Exponential Gaussian Process Regression
        % Adjacency matrix representation
        [trained_GPR_Adj,validationRMSE_GPR_Adj] = trainRegressionModel_GPR(train_Adj);
        y_pred = trained_GPR_Adj.predictFcn(test_Adj(:,1:D));
        RMSE_GPR_Adj(rep,i) = (test_Adj(:,D+1)-y_pred).'*(test_Adj(:,D+1)-y_pred)/n_test;
        display(['Exponential Gaussian Process_Adjacency: RMSE = ',num2str(RMSE_GPR_Adj(rep,i)),', # trait = ',num2str(ind_ctns(i)),', rep = ', num2str(rep)])

        
        % Tree representation
        [trained_GPR_Tree,validationRMSE_GPR_Tree] = trainRegressionModel_GPR(train_Tree);
        y_pred = trained_GPR_Tree.predictFcn(test_Tree(:,1:d));
        RMSE_GPR_Tree(rep,i) = (test_Tree(:,d+1)-y_pred).'*(test_Tree(:,d+1)-y_pred)/n_test;
        display(['Exponential Gaussian Process_Tree: RMSE = ',num2str(RMSE_GPR_Tree(rep,i)),', # trait = ',num2str(ind_ctns(i)),', rep = ', num2str(rep)])

       
    end
end

RMSE_Baseline(:,138) = [];
RMSE_LR_Adj(:,138) = [];
RMSE_LR_Tree(:,138) = [];
RMSE_GPR_Adj(:,138) = [];
RMSE_GPR_Tree(:,138) = [];


Improve_LR_Adj = -(mean(RMSE_LR_Adj,1) - mean(RMSE_Baseline,1))./mean(RMSE_Baseline,1);
Improve_LR_Tree= -(mean(RMSE_LR_Tree,1) - mean(RMSE_Baseline,1))./mean(RMSE_Baseline,1);

Improve_GPR_Adj = -(mean(RMSE_GPR_Adj,1) - mean(RMSE_Baseline,1))./mean(RMSE_Baseline,1);
Improve_GPR_Tree= -(mean(RMSE_GPR_Tree,1) - mean(RMSE_Baseline,1))./mean(RMSE_Baseline,1);



figure
hold on
plot(Improve_LR_Adj,Improve_LR_Tree,'bo')
line([min(min(Improve_LR_Adj),min(Improve_LR_Tree)),max(max(Improve_LR_Adj),max(Improve_LR_Tree))],[min(min(Improve_LR_Adj),min(Improve_LR_Tree)),max(max(Improve_LR_Adj),max(Improve_LR_Tree))])
hold off
xlabel('Adjacency')
ylabel('Tree')
title('Linear Regression Improvement of MSE over Baseline')


figure
hold on
plot(Improve_GPR_Adj,Improve_GPR_Tree,'bo')
line([min(min(Improve_GPR_Adj),min(Improve_GPR_Tree)),max(max(Improve_GPR_Adj),max(Improve_GPR_Tree))],[min(min(Improve_GPR_Adj),min(Improve_GPR_Tree)),max(max(Improve_GPR_Adj),max(Improve_GPR_Tree))])
hold off
xlabel('Adjacency')
ylabel('Tree')
title('GP Regression Improvement of MSE over Baseline')


save('RMSE_Baseline.mat','RMSE_Baseline')
save('RMSE_LR_Adj.mat','RMSE_LR_Adj')
save('RMSE_LR_Tree.mat','RMSE_LR_Tree')
save('RMSE_GPR_Adj.mat','RMSE_GPR_Adj')
save('RMSE_GPR_Tree.mat','RMSE_GPR_Tree')

save('Improve_LR_Adj.mat','Improve_LR_Adj')
save('Improve_LR_Tree.mat','Improve_LR_Tree')
save('Improve_GPR_Adj.mat','Improve_GPR_Adj')
save('Improve_GPR_Tree.mat','Improve_GPR_Tree')
