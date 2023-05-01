%% Init 

A = [4, 2, -1;
     2, 5, 2;
     -1, 2 10];
b = [12; -8; 4];
x_0 = zeros(3, 1)
iter = 10;
tol = 1e-6;
%% test 

r_0 = b - A * x_0;
p_0 = r_0;

p_k = p_0;
x_k = x_0;
r_k = r_0
for k = 0:iter
    a_k = (r_k.' * r_k)/(p_k.' * A * p_k)
    x_k_next = x_k + a_k * p_k
    r_k_next = r_k - a_k * A * p_k
    B_k = (r_k_next.' * r_k_next)/(r_k.' * r_k)
    p_k_next = r_k_next + B_k * p_k

    p_k = p_k_next
    r_k = r_k_next
    x_k = x_k_next
    if norm(p_k) <= tol
        disp(k)
        break
    end
end




