%% % Plot u(x) n=1000

data_1000 = readmatrix('exact_solutions_1000.csv');
x_1000 = data_1000(:,1);
u_1000 = data_1000(:,2);

%continuos function
figure();
plot(x_1000, u_1000)
grid on
xlabel('x');
ylabel('u(x)');
title('Plot of u(x)');
