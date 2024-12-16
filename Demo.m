clc;clear;
warning('off');
tic

%% 1. import data

Wdd = load('./Dataset1/Disease_similarity.csv');            %218*218 
WDD = load('./Dataset1/LncRNA_similarity.csv');             %447*447 
WdD = load('./Dataset1/Disease_LncRNA_association.csv');    %218*447 


%{
load('./Dataset2/wholeset_dissim37.mat'); 
Wdd = wholeset_dissim37;
load('./Dataset2/wholeset_lncdismatrix.mat');  
WdD = wholeset_lncdismatrix';
load('./Dataset2/wholeset_lncsim28.mat');  
WDD = wholeset_lncsim28;
%}

%{
WDD = xlsread('./Dataset3/04-lncRNA-lncRNA.xlsx'); 
WdD = xlsread('./Dataset3/05-lncRNA-disease.xlsx'); 
WdD = WdD';
Wdd = xlsread('./Dataset3/07-disease-disease.xlsx'); 
%}

[dn,dr] = size(WdD);
%% 2. structure KCKA
WdD_temp1 = WdD;
gamma_net = 2^0;
ks = 0.6;
[KD_COM,Kd_COM] = CT(WdD_temp1', WDD, Wdd, gamma_net,ks);

%% 3. structure heterogeneous matrix
WdD_ori = WdD;
maxiter = 300;
alpha = 0.5;   
beta = 5; 
tol1 = 2*1e-3;
tol2 = 1*1e-5;
C = 4;
T = [KD_COM, WdD_temp1'; WdD_temp1, Kd_COM];
[t1, t2] = size(T);
trIndex = double(T ~= 0);
[WW,~] = WNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1, C);
Score = WW((t1-dn+1) : t1, 1 : dr);  

%% 
Score_ori = Score;
index=find(WdD_ori==1);
auc = zeros(1,10);

%%  10-fold CV

for i = 1:5
    i
    indices = crossvalind('Kfold', length(index), 10);
    WdD_temp = WdD_ori;
    Score = Score_ori;
    for j = 1:10
        j
        index_2 = find(j == indices);
        WdD_temp(index(index_2)) = 0; 
        gamma_net = 2^0;
        %[KD_COM_temp, Kd_COM_temp] = CT(WdD_temp',WDD,Wdd,gamma_net,ks);
        %T_1 = [KD_COM_temp, WdD_temp'; WdD_temp, Kd_COM_temp];
        T_1 = [KD_COM, WdD_temp'; WdD_temp, Kd_COM];
        [t1_1, t2_1] = size(T_1);
        trIndex_1 = double(T_1 ~= 0);
        [WW_1,iter_1] = WNNR(alpha, beta, T_1, trIndex_1, tol1, tol2, maxiter, 0, 1, C);
        M_recovery_1 = WW_1((t1_1-dn+1) : t1_1, 1 : dr);  
        Score(index(index_2)) = M_recovery_1(index(index_2));
        WdD_temp = WdD;
    end
        pre_label_score = Score(:);
        label_y = WdD_ori(:);
        auc(i) = roc_1(pre_label_score,label_y,'red');
end

auc_ave = mean(auc);
auc_std = std(auc);

%% global loocv
%{
step = sum(WdD(:));
WdD_temp = WdD_ori;
for i = 1 : step
   i
   WdD_temp(index(i))=0;
   [KD_COM,Kd_COM] = CT(WdD_temp', WDD, Wdd, gamma_net,ks);
   T_1 = [KD_COM, WdD_temp'; WdD_temp, Kd_COM]; 
   [t1_1, t2_1] = size(T_1);
   trIndex_1 = double(T_1 ~= 0);
   [WW_1,iter_1] = WNNR(alpha, beta, T_1, trIndex_1, tol1, tol2, maxiter, 0, 1, C);
   M_recovery_1 = WW_1((t1_1-dn+1) : t1_1, 1 : dr);
   Score(index(i)) = M_recovery_1(index(i));
   WdD_temp = WdD;
end
pre_label_score = Score(:);
label_y = WdD_ori(:);
auc=roc_1(pre_label_score,label_y,'red');

auc_ave = mean(auc);
auc_std = std(auc);
%}

%%
toc
