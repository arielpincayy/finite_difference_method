function [mesh] = fdm_implicit_animation(top, bottom, left, right, c_w, u_i, f, l_i, l_j, delta_t, t_end, z_range)
    
    n = size(left, 1) + 2;  
    m = size(top, 2);  

    num_steps = t_end / delta_t;

    z_min = z_range(1);
    z_max = z_range(2);

    X = linspace(0,l_j,m);
    Y = linspace(0,l_i,n);
    T = linspace(delta_t,t_end, num_steps);
    its = zeros(1,num_steps-1);

    h = X(2) - X(1);
    k = Y(2) - Y(1);


    errors = zeros(1,num_steps);

    rx = (c_w * delta_t / h)^2;
    ry = (c_w * delta_t / k)^2;
    
    a = 1/(c_w*delta_t)^2;
    b = 1/(4*h^2);
    c = 1/(4*k^2);

    s = a+(2*b)+(2*c);
    
    mov = [0,0; 1,0; 0,1; -1,0; 0,-1];
    terms = [s,-b,-c,-b,-c];

    stencil = zeros(n,m);

    pos = 1;
    for i=2:n-1
        for j=2:m-1
            stencil(i,j) = pos;
            pos = pos + 1;
        end
    end

    for i=1:m
        stencil(1,i) = pos;
        stencil(end,i) = (m-1) + pos;
        pos = pos + 1; 
    end

    for i=2:n-1
        stencil(i,1) = pos;
        stencil(i,end) = (n-2) + pos;
        pos = pos + 1;
    end

    mesh = zeros(n, m);
    df_mesh = zeros(n,m);
    if(isa(u_i, 'function_handle'))
        for i=2:n-1
            for j=2:m-1
                mesh(i,j) = u_i(X(j),Y(i),0);
                df_mesh(i,j) = df_mesh(i,j) = (u_i(X(j),Y(i),0.00000001) - u_i(X(j),Y(i),0))/0.00000001;
            end
        end
    else
        mesh = u_i;
    end

    mesh(1, :) = top;            
    mesh(end, :) = bottom;       
    mesh(2:end-1, 1) = left;     
    mesh(2:end-1, end) = right;  

    new_mesh = zeros(n,m);
 
    % Fórmula especial para el primer paso usando la derivada inicial
    for i = 2:n-1
        for j = 2:m-1
            new_mesh(i,j) = mesh(i,j) + delta_t * df_mesh(i,j) + 0.5 * (rx * (mesh(i+1,j) + mesh(i-1,j)) + ry * (mesh(i,j+1) + mesh(i,j-1)) - 2 * (rx + ry) * mesh(i,j)) + 0.5 * delta_t^2 * f(X(j),Y(i),T(1));
        end
    end


    [u_exact] = exact(c,l_i,l_j,delta_t,n,m,u_i,0);
    errors(1,1) = norm(abs(u_exact - new_mesh));

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


    for iter=2:num_steps

        matrix = zeros((n-2)*(m-2),(n-2)*(m-2));
        alpha = zeros((n-2)*(m-2),1);

        counter = 1;
        for i=2:n-1
            for j=2:m-1
                alpha_n = (2*a - 4*b - 4*c)*mesh(i,j) - s*(old_mesh(i,j)) + 2*b*mesh(i-1,j) + 2*b*mesh(i+1,j) + 2*c*mesh(i,j-1) + 2*c*mesh(i,j+1) + c*(old_mesh(i,j-1)) + c*(old_mesh(i,j+1)) + b*(old_mesh(i-1,j)) + b*(old_mesh(i+1,j)) + f(X(j),Y(i),T(iter));
                for k_i=1:5
                    pos = stencil(i + mov(k_i,1),j + mov(k_i,2));
                    num = terms(k_i);

                    if(pos > (n-2)*(m-2))
                        alpha_n = alpha_n - num*mesh(i + mov(k_i,1),j + mov(k_i,2));
                    else
                        matrix(counter,pos) = num;
                    end
                end
                alpha(counter,1) = alpha_n;
                counter = counter + 1;
            end
        end

        [sol, it] = gaussSeidel(matrix, alpha, 1e-6, 200);
        its(1,iter-1) = it;
    

        p = 1;
        for i=2:n-1
            for j=2:m-1
                new_mesh(i,j) = sol(p,1);
                p = p + 1;
            end
        end
        
        old_mesh = mesh;  % u^k -> old_mesh
        mesh = new_mesh;  % u^{k+1} -> mesh

        [u_exact] = exact(c,l_i,l_j,T(iter),n,m,u_i,0);
        errors(1,iter) = norm(abs(u_exact - new_mesh));

        if(mod(iter,5)==0 || iter==num_steps || iter==1)
            surf(X,Y,new_mesh);
            set(gca, 'FontSize', 16);
            title(['Wave equation on t = ', num2str(iter * delta_t)]);
            xlabel('Eje X');
            ylabel('Eje Y');
            zlabel('u(x, y, t)');
            view(45,30);
            colorbar;
    
            zlim([z_min, z_max]);  
    
            pause(0.01);
        end

    end

    figure;
    plot(T,errors);
    set(gca, 'FontSize', 16);
    title('Error curve. Implicit method');
    xlabel('Time (seconds)');
    ylabel('Error');
    colorbar;

    figure'
    plot(its);
    set(gca, 'FontSize', 16);
    title('Number of iterations');
    xlabel('Time (Seconds)');
    ylabel('Iterations');
    colorbar;

    %Grafica la solucion exacta, y el error.
    [u_exact] = exact(c,l_i,l_j,t_end,n,m,u_i,1);
    [error_mesh] = error_graph(u_exact, mesh, n, m, l_i, l_j);
end