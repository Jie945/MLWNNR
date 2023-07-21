 function [knn_network] = knn(P,k)
    [~,n] = size(P);
    knn_network = zeros(n, n);
    [sort_network,idx]=sort(P,2,'descend');
    for i = 1 : n
       knn_network(i,idx(i,1:k))=sort_network(i,1:k);
    end
 
end