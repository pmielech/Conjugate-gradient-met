function [x, k]= conjugate_grad(A, b, x_0, iter, to_plot, f, tol)

col = ["black", "red", "yellow"];
r_0 = b - A * x_0;
p_0 = r_0;
p1 = zeros(1, 3);
p2 = 0;
p3 = 0;

p = p_0;
x = x_0;
r = r_0;
for k = 0:iter

    if(to_plot)
        if(k == 0)
            p1 = [x(1),x(2),f(x(1),x(2))];
            text(p1(1),p1(2), f(x(1),x(2)) + 2, num2str(k), "Color", 'r', "HorizontalAlignment","left", ...
                'BackgroundColor',"w")
        end
        
        if(k == 1)
            p2 = [x(1),x(2),f(x(1),x(2))];
            text(p2(1) + 0.5,p2(2), f(x(1),x(2)) + 2, num2str(k), "Color", 'r', "HorizontalAlignment","left", ...
                'BackgroundColor',"w")
        end

        axis square
        hold on
        z = f(x(1), x(2))
        scatter3(x(1),x(2),f(x(1),x(2)),col(k+1),'fill', "DisplayName", "Conjugate Gradient, iter" + num2str(k)) 

        if k == 2
        quiver3(p1(1),p1(2),p1(3),p2(1)-p1(1),p2(2),p2(3)-p1(3),'r','linewidth',2, "DisplayName","[iter0, iter1]");
        quiver3(p2(1),p2(2),p2(3),p3(1)-p2(1),p3(2)-p2(2),p3(3)-p2(3),'r','linewidth',2, "DisplayName","[iter1, iter2]");
        end 
    end
    a_k = (r.' * r)/(p.' * A * p);
    x_next = x + a_k * p;
    r_next = r - a_k * A * p;
    
    if norm(r_next) <= tol   % the convergence condition
        k = k +1;
        if to_plot
        p3 = [x_next(1),x_next(2),f(x_next(1),x_next(2))];
        text(p3(1) - 0.5,p3(2), f(x_next(1),x_next(2)) + 2, num2str(k), "Color", 'r', "HorizontalAlignment","right", ...
            'BackgroundColor',"w")
        z = f(x_next(1), x_next(2))
        scatter3(x_next(1),x_next(2),z,col(k+1),'fill', "DisplayName", "Conjugate Gradient, iter" + num2str(k)) 
        quiver3(p1(1),p1(2),p1(3),p2(1)-p1(1),p2(2),p2(3)-p1(3),'r','linewidth',2, "DisplayName","[iter0, iter1]");
        quiver3(p2(1),p2(2),p2(3),p3(1)-p2(1),p3(2)-p2(2),p3(3)-p2(3),'r','linewidth',2, "DisplayName","[iter1, iter2]");
        
        figure
        plot([p1(1), p2(1), p3(1)],[p1(2), p2(2), p3(2)], '-o')
        hold on 
        end
        break;
    end

    B_k = (r_next.' * r_next)/(r.' * r);
    p_next = r_next + B_k * p;

    p = p_next
    r = r_next
    x = x_next
    
end
    
end
