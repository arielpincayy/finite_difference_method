function [x, iter] = gaussSeidel(A, b, tol, maxIter)
    
    n = length(b);           
    x = zeros(n,1);                 
    iter = 0;                
    error = tol + 1;         

    while error > tol && iter < maxIter
        x_old = x;           
        for i = 1:n
            sum1 = A(i, 1:i-1) * x(1:i-1);  
            sum2 = A(i, i+1:n) * x_old(i+1:n); 
            x(i) = (b(i) - sum1 - sum2) / A(i, i);
        end
        iter = iter + 1;        
        error = norm(x - x_old); 
    end
end
