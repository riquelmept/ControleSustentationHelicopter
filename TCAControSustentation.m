A = [0 1; -81 -54]
B = [0; 1]
C = [54 9/10]
polos = eig(A);
polos_desejados = 5 * polos;
L = place(A', C', polos_desejados)'
D = 0
sys = ss(A, B, C, D);
t = 0:0.01:10 % Definindo o vetor de tempo
u = ones(size(t)) % Entrada de degrau unitário

% Simulando a resposta do sistema à entrada de degrau
[y, t, x] = lsim(sys, u, t)
% Inicialização do estado estimado
x_hat = zeros(2, length(t))
y_hat = zeros(1, length(t))

for i = 1:length(t) - 1
    % Atualizando o estado estimado com base na entrada e na saída medida
    dx_hat = (A - L*C)*x_hat(:, i) + B*u(i) + L*y(i);
    x_hat(:, i+1) = x_hat(:, i) + dx_hat * (t(i+1) - t(i));
    % Calculando a saída estimada
    y_hat(i) = C*x_hat(:, i);
end
y_hat(end) = C*x_hat(:, end) % Calculando a última saída estimada
figure;
subplot(2,1,1);
plot(t, x);
title('Estados Reais');
xlabel('Tempo (s)');
ylabel('Estados');

subplot(2,1,2);
plot(t, x_hat');
title('Estados Estimados pelo Observador');
xlabel('Tempo (s)');
ylabel('Estimados');

legend({'Estado 1', 'Estado 2'}, 'Location', 'best');