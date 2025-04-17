function [u_exact] = exact(c,l_i,l_j,delta_t,n,m,u_i,g)

    u_exact = zeros(n,m);

    X = linspace(0,l_j,m);
    Y = linspace(0,l_i,n);

    h = X(2) - X(1);
    k = Y(2) - Y(1);
    
    for i=1:n
        for j=1:m
            u_exact(i,j) = u_i(X(j),Y(i),delta_t);
        end
    end

    if(g)
        figure;
        surf(X,Y,u_exact);
        set(gca, 'FontSize', 16);
        title('Exact solution');
        xlabel('X axis');
        ylabel('Y axis');
        zlabel('u(x, y, t)');
        view(45, 30);  
        colorbar;
    end
end
