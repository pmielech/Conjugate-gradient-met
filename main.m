%% 3x3 matrix 
close all;
f = @(x1,x2) 0.5 .* [x1;x2].' * A * [x1;x2] - b.' * [x1;x2] + 0.5 .* b.' * A^(-1) * b;

A = [4, 2, -1;
    2, 5, 2;
    -1, 2, 10];
b = [12; -8; 4];
c = [5.8321; -4.7023; 1.7982];

tol = 1e-6;
x_0 = zeros(size(b));
iter = 1000;
to_plot = false;

disp("Matrix 3x3 calculation - Conjugate gradient method:")
[x_conj, iter_conj] = conjugate_grad(A, b, x_0, iter, to_plot, f, tol);
disp("Number of iterations:")
disp(iter_conj)

%% 5x5 matrix 

A = [5 -1 0 0 0; 
     -1 4 -1 0 0; 
     0 -1 3 -1 0; 
     0 0 -1 4 -1; 
     0 0 0 -1 2];

b = [1; 2; 3; 4; 5];

x_0 = zeros(size(b));

disp("Matrix 5x5 calculation - Conjugate gradient method:")
[x_conj, iter_conj] = conjugate_grad(A, b, x_0, iter, to_plot, f, tol);
disp("Number of iterations:")
disp(iter_conj)
%% 16x16 matrix

A = gallery('poisson', 4) * 100;
b = ones(16, 1);

x_0 = zeros(size(b));

disp("Matrix 16x16 calculation - Conjugate gradient method:")
[x_conj, iter_conj] = conjugate_grad(A, b, x_0, iter, to_plot, f, tol);
disp("Number of iterations:")
disp(iter_conj)

%% 2x2 matrix with plot 

A = [5, 2;
     2, 3];
b = [-2; 4];

x_0 = zeros(size(b));
to_plot = true;

fsurf(f);
legend("f(x)")
xlabel("x1");
ylabel("x2");
zlabel("z");

% conjugate_grad
disp("Matrix 2x2 calculation - Conjugate gradient method:")
[x_conj, iter_conj] = conjugate_grad(A, b, x_0, iter, to_plot, f, tol);
disp("Number of iterations:")
disp(iter_conj)

% gauss
disp("Matrix 2x2 calculation - Gauss-Seidl method:")
[x_gauss, iter_gauss] = gauss_seid(A, b, x_0, iter, f, tol,to_plot)
disp("Number of iterations:")
disp(iter_gauss)
