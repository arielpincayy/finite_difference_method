function [error_mesh] = error_graph(u_exact, u_aprox, n, m, l_i, l_j)
    
    error_mesh = abs(u_exact - u_aprox);

    X = linspace(0,l_j,m);
    Y = linspace(0,l_i,n);


    figure;
    surf(X,Y,error_mesh);
    set(gca, 'FontSize', 16);
    title('Error');
    xlabel('Eje X');
    ylabel('Eje Y');
    zlabel('u(x, y, t)');
    view(45,30)
    colorbar;

end