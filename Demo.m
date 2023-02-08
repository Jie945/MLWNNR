clc;clear;
warning('off');
tic

%% 1. 导入数据
Wdd = load('./dataset/Disease_similarity.csv');   %218*218 疾病
WDD = load('./dataset/LncRNA_similarity.csv');  %447*447  RNA
WdD = load('./dataset/Disease_LncRNA_association.csv');     %218*447 疾病*RNA
[dn,dr] = size(WdD);

%% 2. 计算相似性和核融合
WdD_temp1 = WdD;
gamma_net = 2^0;
ks = 0.6;
[KD_COM,Kd_COM] = CT(WdD_temp1', WDD, Wdd, gamma_net,ks);

%% 3. 计算一次预测
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

%% 测试前准备
Score_ori = Score;
index=find(WdD_ori==1);
auc = zeros(1,10);

%%  5-fold CV
%{
for i = 1:1
    i
    indices = crossvalind('Kfold', length(index), 5);
    WdD_temp = WdD_ori;
    Score = Score_ori;
    for j = 1:5
        j
        index_2 = find(j == indices);
        WdD_temp(index(index_2)) = 0; 
        gamma_net = 2^0;
        [KD_COM_temp, Kd_COM_temp] = CT(WdD_temp',WDD,Wdd,gamma_net,ks);
        T_1 = [KD_COM_temp, WdD_temp'; WdD_temp, Kd_COM_temp];
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
%}
%% global loocv

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


%%
toc
