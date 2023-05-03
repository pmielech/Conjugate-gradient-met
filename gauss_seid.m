function [x, iter] = gauss_seid(A, b, x0, max_iter,f , tol, to_plot)

n = length(b);
x = x0;
err = zeros(max_iter, 2);
iter = 0;
converged = false;
plot_p = {};
surf_plot = false;

while ~converged && iter < max_iter
    
    plot_p = [plot_p, [x(1), x(2)]];
    iter = iter + 1;
    for i = 1:n
        sigma = A(i,1:i-1) * x(1:i-1) + A(i,i+1:n) * x0(i+1:n);
        x(i) = (b(i) - sigma) / A(i,i);
    end
    
    if surf_plot
        axis square
        hold on
        f(x(1), x(2))
        scatter3(x(1),x(2),f(x(1),x(2)),'fill', "DisplayName", "Gauss-Seidle, iter" + num2str(iter))         
    end

    err(iter) = norm(A*x - b, 2);   
    if err(iter) < tol              % the convergence condition
        converged = true;
    end
    x0 = x; 
end

if to_plot
    X = cell2mat(cellfun(@(p) p(1), plot_p, 'UniformOutput', false));
    Y = cell2mat(cellfun(@(p) p(2), plot_p, 'UniformOutput', false));
    plot(X, Y, '-o');
    xlabel("x1")
    ylabel("x2")
    legend("Conjugate gradient", "Gauss-Seidle")
end
end

