function [K_COM1,K_COM2] = CT(WDd,WDD,Wdd,gamma_net,ks)

K1 = [];
K1(:,:,1) = (WDD);                         % RNA表达相似性
K1(:,:,2) = kernel_gip(WDd,1, gamma_net);  % Gaussion
K1(:,:,3) = LNS(WDd,0,90,'regulation2');   % Linear neighborhood similarity
K1(:,:,4) = lncRNAfunsim(Wdd,WDd);         % function similarity


K2 = [];
K2(:,:,1) = (Wdd);                          % 疾病功能相似性
K2(:,:,2) = kernel_gip(WDd,2, gamma_net);   % Gaussion
K2(:,:,3) = LNS(WDd',0,170,'regulation2');  % Linear neighborhood similarity


%%
[weight_v1] = cka_kernels_weights(K1,WDd,1,ks);
K_COM1 = combine_kernels(weight_v1, K1);		

[weight_v2] = cka_kernels_weights(K2,WDd,2,ks);
K_COM2 = combine_kernels(weight_v2, K2);

end
