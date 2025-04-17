function [mesh] = fdm_explicit_animation(top, bottom, left, right, c, u_i, f, l_i, l_j, delta_t, t_end, z_range)

    n = size(left, 1) + 2;  
    m = size(top, 2);   

    num_steps = floor(t_end / delta_t);

    X = linspace(0, l_j, m);
    Y = linspace(0, l_i, n); 
    T = linspace(delta_t, t_end, num_steps);

    z_min = z_range(1);
    z_max = z_range(2);
    
    k = X(2) - X(1);
    h = Y(2) - Y(1);

    errors = zeros(1,num_steps);

    rx = (c * delta_t / h)^2;
    ry = (c * delta_t / k)^2;
        
    mesh = zeros(n, m);  % u^0
    df_mesh = zeros(n, m); % u^1 (derivada inicial)

    % Condiciones iniciales y derivadas en t = 0
    if(isa(u_i, 'function_handle'))
        for i = 2:n-1
            for j = 2:m-1
                mesh(i,j) = u_i(X(j),Y(i),0);  % u_0
                df_mesh(i,j) = (u_i(X(j),Y(i),0.00000001) - u_i(X(j),Y(i),0))/0.00000001;  
            end
        end
    else
        mesh = u_i;
    end

    % Condiciones de frontera
    mesh(1, :) = top;            
    mesh(end, :) = bottom;       
    mesh(2:end-1, 1) = left;     
    mesh(2:end-1, end) = right;  

    % Crear el primer paso (t = delta_t)
    new_mesh = zeros(n, m);  % u^1 (primer paso temporal)

    % Fórmula especial para el primer paso usando la derivada inicial
    for i = 2:n-1
        for j = 2:m-1
            new_mesh(i,j) = mesh(i,j) + delta_t * df_mesh(i,j) + 0.5 * (rx * (mesh(i+1,j) + mesh(i-1,j)) + ry * (mesh(i,j+1) + mesh(i,j-1)) - 2 * (rx + ry) * mesh(i,j)) + 0.5 * delta_t^2 * f(X(j),Y(i),T(1));
        end
    end

    [u_exact] = exact(c,l_i,l_j,delta_t,n,m,u_i,0);
    errors(1,1) = max(max((abs(u_exact - new_mesh))));

    % Dibujar la solución para el primer paso
    figure;
    surf(X,Y,new_mesh);
    set(gca, 'FontSize', 16);
    title(['Wave equation on t = ', num2str(delta_t)]);
    xlabel('Eje X');
    ylabel('Eje Y');
    zlabel('u(x, y, t)');
    view(45,30)
    colorbar;
    zlim([z_min, z_max]);
    pause(0.01);

    % Actualizar para los siguientes pasos
    old_mesh = mesh;  % u^0 -> old_mesh
    mesh = new_mesh;  % u^1 -> mesh (será el actual en el siguiente paso)

    % Iterar sobre los siguientes pasos temporales
    for t = 2:num_steps
        new_mesh = zeros(n, m);

        % Fórmula explícita para los pasos siguientes (t > delta_t)
        for i = 2:n-1
            for j = 2:m-1
                new_mesh(i,j) = 2*mesh(i,j) - old_mesh(i,j) + rx*(mesh(i+1,j) + mesh(i-1,j)) + ry*(mesh(i,j+1) + mesh(i,j-1)) - 2*(rx + ry)*mesh(i,j) + delta_t^2 * f(X(j),Y(i),T(t));
            end
        end

        [u_exact] = exact(c,l_i,l_j,T(t),n,m,u_i,0);
        errors(1,t) = max(max((abs(u_exact - new_mesh))));

        % Dibujar la solución
        surf(X,Y,new_mesh);
        set(gca, 'FontSize', 16);
        title(['Wave equation on t = ', num2str(t * delta_t)]);
        xlabel('Eje X');
        ylabel('Eje Y');
        zlabel('u(x, y, t)');
        view(45,30)
        colorbar;
        zlim([z_min, z_max]);
        pause(0.01);

        % Actualizar para el siguiente paso
        old_mesh = mesh;  % u^k -> old_mesh
        mesh = new_mesh;  % u^{k+1} -> mesh
    end


    figure;
    plot(T,errors);
    set(gca, 'FontSize', 16);
    title('Error curve. Explicit method');
    xlabel('Time (seconds)');
    ylabel('Error');
    colorbar;
    
    %Grafica la solucion exacta, y el error.
    [u_exact] = exact(c,l_i,l_j,t_end,n,m,u_i,1);
    [error_mesh] = error_graph(u_exact, mesh, n, m, l_i, l_j);
end
