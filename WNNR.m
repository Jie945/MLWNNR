function [T_recovery, iter] = WNNR(alpha, beta, T, trIndex, tol1, tol2, maxiter, a, b, c)

X = T;
W = X;
Y = X;
i = 1;
stop1 = 1;
stop2 = 1;

while(stop1 > tol1 || stop2 > tol2)
    
    %the process of computing W
    tran = (1/beta) * (Y + alpha * (T.* trIndex)) + X;
    W = tran - (alpha/ (alpha + beta)) * (tran.* trIndex);
    W(W < a) = a;
    W(W > b) = b;
%%
    
    [U, D, V, flag] = svds(W- 1/beta* Y);
    if flag == 1
        flag
    end
    m = 0.00000001;
    c1 = diag(D) - m;
    c2 = (diag(D) + m).^2 - (4 * c);
    E = ((c1 + sqrt(abs(c2))) / 2) .* (c2 > 0); 
    E = max(E - 1/beta, 0);
    X_1 = U * diag(E) * V';
    
 %%  
    %the process of computing Y
    Y = Y + beta * (X_1 - W);
    
    
    stop1_0 = stop1;
    stop1 = norm(X_1 - X, 'fro') / norm(X, 'fro');
    stop2 = abs(stop1 - stop1_0)/ max(1, abs(stop1_0));
    
    X = X_1;
    i = i+1;
    
    if i < maxiter
        iter = i - 1;
    else
        iter = maxiter;
        warning('reach maximum iteration~~do not converge!!!');
        break
    end
    
end

T_recovery = W;

end
