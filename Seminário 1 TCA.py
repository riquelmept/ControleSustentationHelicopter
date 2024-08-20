import numpy as np
from scipy.signal import place_poles, StateSpace, step, lti
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

"""
Sistemas de controle de sustentação em helicópteros tem por objetivo manter o direcionamento da aeronave frente a presença de distúrbios,
tendo como entrada o nível de acionamento do pedal (NP) e a saída a velocidade angular no rotor (Wr).
Considerando os estados velocidade do vento (Vv) e frequência da rajada (fv), 𝑥=[𝑉𝑣; 𝑓𝑣], uma parte do movimento específico desse controle possui a dinâmica dada por,
𝐴=[0 1; −81 −54], 𝐵=[0; 1], 𝐶=[54 9/10]
"""

"""
Projete um observador de estados capaz de emular o sistema 5x mais rápido. 
Demonstre, de forma comparativa entre o sistema original e o observador, o comportamento de cada estado;
"""
# Definição das Matrizes do Sistema
A = np.array([[0, 1], [-81, -54]])
B = np.array([[0], [1]])
C = np.array([[54, 9/10]])
D = np.array([[0]])

# Calcular os autovalores do sistema
polos = np.linalg.eigvals(A)
print('Polos do Sistema')
print(polos)

# Projeto do Observador (polos desejados 5x mais rápidos)
polos_desejados = 5 * polos
print('Polos Desejados')
print(polos_desejados)

# Calcular a matriz de ganhos do observador L usando place_poles
place_obj = place_poles(A.T, C.T, polos_desejados)
L_1 = place_obj.gain_matrix.T
print(L_1)

# Nova matriz A do observador
AOB_1 = A - L_1 @ C
print(AOB_1)

# Simulação da Resposta ao Degrau
t = np.linspace(0, 10, 1000)

# Sistema original
sys_GMA = StateSpace(A, B, C, D)
t_GMA, y_GMA = step(sys_GMA, T=t)

# Sistema observador
sys_GOB_1 = StateSpace(AOB_1, B, C, D)
t_GOB_1, y_GOB_1 = step(sys_GOB_1, T=t)

# Sistema observador compensado
C_compensado = 5 * C
sys_GOB = StateSpace(AOB_1, B, C_compensado, D)
t_GOB, y_GOB = step(sys_GOB, T=t)

# Plotando as respostas ao degrau em um único gráfico para comparar diretamente

plt.figure(figsize=(10, 6))

# Plot da resposta do sistema Malha Aberta
plt.plot(t_GMA, y_GMA, label="Sistema Malha Aberta (GMA)", color="blue")

# Plot da resposta do sistema Observador
plt.plot(t_GOB_1, y_GOB_1, label="Sistema Observador (GOB_1)", color="red")

# Plot da resposta do sistema Observador Compensado
plt.plot(t_GOB, y_GOB, label="Sistema Observador Compensado (GOB)", color="green")

plt.xlabel('Tempo (s)')
plt.ylabel('Resposta ao Degrau')
plt.title('Comparação entre o Sistema Original, Observador, e Observador Compensado')
plt.legend()
plt.grid(True)
plt.show()

"""
Através da fórmula de Ackermann, obtenha o vetor de ganhos K capaz de fazer o tempo de acomodação (critério de 2%) ser o mais rápido possível,
com sobressinal limitado a 10%, tanto para o sistema original quando para sistema do observador;
"""
#Letra B

# 1. Ganhos Iniciais (MP = 0.1, Ts = 1s)

# Definindo os requisitos de desempenho iniciais
Mp_initial = 0.10  # Sobressinal máximo de 10%
zeta_initial = -np.log(Mp_initial) / np.sqrt(np.pi**2 + np.log(Mp_initial)**2)  # Fator de amortecimento associado a Mp

Ts_initial = 1.0  # Tempo de acomodação desejado em segundos
omega_n_initial = 4 / (zeta_initial * Ts_initial)  # Frequência natural associada

# Calculando os polos desejados para essa condição inicial
p1_initial = -zeta_initial * omega_n_initial + 1j * omega_n_initial * np.sqrt(1 - zeta_initial**2)
p2_initial = -zeta_initial * omega_n_initial - 1j * omega_n_initial * np.sqrt(1 - zeta_initial**2)
polos_desejados_initial = [p1_initial, p2_initial]

# Aplicando Ackermann para obter o K correspondente
K_original_initial = place_poles(A, B, polos_desejados_initial).gain_matrix
print("K Original Inicial:")
print(K_original_initial)

# Para o sistema com observador
L_1_initial = place_poles(A.T, C.T, 5 * np.array(polos_desejados_initial)).gain_matrix.T
A_observer_initial = A - L_1_initial @ C
K_observer_initial = place_poles(A_observer_initial, B, polos_desejados_initial).gain_matrix
print("K Observador Inicial:")
print(K_observer_initial)

# Simulação da Resposta ao Degrau com os ganhos iniciais

# Sistema original com ganho inicial
A_K_initial = A - B @ K_original_initial
sys_GMA_initial = StateSpace(A_K_initial, B, C, 0)
t_GMA_initial, y_GMA_initial = step(sys_GMA_initial, T=np.linspace(0, 10, 1000))

print("AK Original Inicial:")
print(A_K_initial)

# Sistema observador com ganho inicial
AOB_K_initial = A_observer_initial - B @ K_observer_initial
sys_GOB_initial = StateSpace(AOB_K_initial, B, C, 0)
t_GOB_initial, y_GOB_initial = step(sys_GOB_initial, T=np.linspace(0, 10, 1000))
print("A K Observador Inicial:")
print(AOB_K_initial)

# 2. Ganhos Novos (MP = 0.05, Ts = 0.5s)

# Ajustando o fator de amortecimento para reduzir o Mp
Mp_new = 0.05  # Reduzir o Mp para 5%
zeta_new = -np.log(Mp_new) / np.sqrt(np.pi**2 + np.log(Mp_new)**2)

# Ajustando a frequência natural para reduzir o Ts
Ts_new = 0.5  # Tempo de acomodação menor
omega_n_new = 4 / (zeta_new * Ts_new)

# Recalculando os polos desejados
p1_new = -zeta_new * omega_n_new + 1j * omega_n_new * np.sqrt(1 - zeta_new**2)
p2_new = -zeta_new * omega_n_new - 1j * omega_n_new * np.sqrt(1 - zeta_new**2)
polos_desejados_new = [p1_new, p2_new]

# Aplicação da Fórmula de Ackermann para o sistema original com novos parâmetros
K_original_new = place_poles(A, B, polos_desejados_new).gain_matrix
print('K Novo Original')
print(K_original_new)

# Aplicação da Fórmula de Ackermann para o sistema com observador
L_1_new = place_poles(A.T, C.T, 5 * np.array(polos_desejados_new)).gain_matrix.T  # Observador 5x mais rápido
A_observer_new = A - L_1_new @ C
K_observer_new = place_poles(A_observer_new, B, polos_desejados_new).gain_matrix
print('K Novo Observador')
print(K_observer_new)

# Simulação da Resposta ao Degrau com os ganhos novos

# Sistema original com ganho novo
A_K_new = A - B @ K_original_new
sys_GMA_new = StateSpace(A_K_new, B, C, 0)
t_GMA_new, y_GMA_new = step(sys_GMA_new, T=np.linspace(0, 10, 1000))
print('A_K_new')
print(A_K_new)

# Sistema observador com ganho novo
AOB_K_new = A_observer_new - B @ K_observer_new
sys_GOB_new = StateSpace(AOB_K_new, B, C, 0)
t_GOB_new, y_GOB_new = step(sys_GOB_new, T=np.linspace(0, 10, 1000))
print('AOB_K_new')
print(AOB_K_new)

# 3. Plotagem das respostas

plt.figure(figsize=(10, 6))

# Plot da resposta do sistema original com ganhos iniciais e novos
plt.plot(t_GMA_initial, y_GMA_initial, label="Original com Ganho Inicial", color="blue")
plt.plot(t_GMA_new, y_GMA_new, label="Original com Ganho Novo", color="cyan")

# Plot da resposta do sistema observador com ganhos iniciais e novos
plt.plot(t_GOB_initial, y_GOB_initial, label="Observador com Ganho Inicial", color="red")
plt.plot(t_GOB_new, y_GOB_new, label="Observador com Ganho Novo", color="orange")

plt.xlabel('Tempo (s)')
plt.ylabel('Resposta ao Degrau')
plt.title('Comparação entre Ganhos Iniciais e Novos: Sistema Original vs Observador')
plt.legend()
plt.grid(True)
plt.show()


"""
Simule o sistema (original e do observador) antes e depois da aplicação do controle, 
excitando-o com um degrau unitário. Obs: 𝐺𝑀𝐴=𝐶(𝑠𝐼−𝐴)−1𝐵 , 𝐺𝑀𝐹=𝐶(𝑠𝐼−𝐻)−1𝐵; onde 𝐻=𝐴−𝐵𝐾;
"""
#Letra C

H_initial = A - B @ K_original_initial
H_new = A - B @ K_original_new

# Funções de Transferência
# Sistema Original
sys_GMA_initial = lti(A, B, C, 0)
sys_GMF_initial = lti(H_initial, B, C, 0)

sys_GMA_new = lti(A, B, C, 0)
sys_GMF_new = lti(H_new, B, C, 0)

# Simulação da Resposta ao Degrau
t = np.linspace(0, 10, 1000)

# Simulação para os ganhos iniciais
t_GMA_initial, y_GMA_initial = step(sys_GMA_initial, T=t)
t_GMF_initial, y_GMF_initial = step(sys_GMF_initial, T=t)

# Simulação para os ganhos novos
t_GMA_new, y_GMA_new = step(sys_GMA_new, T=t)
t_GMF_new, y_GMF_new = step(sys_GMF_new, T=t)

# Plotagem das Respostas
plt.figure(figsize=(12, 8))

# Ganhos Iniciais
plt.subplot(2, 1, 1)
plt.plot(t_GMA_initial, y_GMA_initial, label="Sistema Original (GMA) - Ganho Inicial", color="blue")
plt.plot(t_GMF_initial, y_GMF_initial, label="Sistema com Controle (GMF) - Ganho Inicial", color="cyan")
plt.title('Resposta ao Degrau - Ganhos Iniciais (MP = 0.1, Ts = 1s)')
plt.xlabel('Tempo (s)')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

# Ganhos Novos
plt.subplot(2, 1, 2)
plt.plot(t_GMA_new, y_GMA_new, label="Sistema Original (GMA) - Ganho Novo", color="red")
plt.plot(t_GMF_new, y_GMF_new, label="Sistema com Controle (GMF) - Ganho Novo", color="orange")
plt.title('Resposta ao Degrau - Ganhos Novos (MP = 0.05, Ts = 0.5s)')
plt.xlabel('Tempo (s)')
plt.ylabel('Amplitude')
plt.legend()
plt.grid(True)

plt.tight_layout()
plt.show()


"""
Obtenha os índices de desempenho no transiente e o erro de regime para o sistema antes e depois da aplicação do controle.
"""
#Letra D

# Função para calcular os índices de desempenho
def performance_indices(t, y):
    y_final = y[-1]
    
    # Sobressinal máximo (Mp)
    Mp = (np.max(y) - y_final) / y_final * 100  # em porcentagem
    
    # Tempo de pico (Tp)
    Tp_index = np.argmax(y)
    Tp = t[Tp_index]
    
    # Tempo de acomodação (Ts) - Considerando a faixa de 2%
    Ts_indices = np.where(np.abs(y - y_final) <= 0.02 * y_final)[0]
    Ts = t[Ts_indices[-1]] if len(Ts_indices) > 0 else np.nan
    
    # Erro em regime permanente
    erro_regime = 1 - y_final  # Degrau unitário, então o erro é 1 - valor final
    
    return Mp, Tp, Ts, erro_regime

# Função para simular e calcular os índices de desempenho
def simulate_and_evaluate(sys):
    t = np.linspace(0, 10, 1000)
    t_out, y_out = step(sys, T=t)
    Mp, Tp, Ts, erro_regime = performance_indices(t_out, y_out)
    return Mp, Tp, Ts, erro_regime

# Simulação para Ganhos Iniciais

# Sistema original com ganho inicial
sys_GMA_initial = lti(A, B, C, 0)
H_initial = A - B @ K_original_initial
sys_GMF_initial = lti(H_initial, B, C, 0)

# Sistema observador com ganho inicial
sys_GMA_initial_perf = simulate_and_evaluate(sys_GMA_initial)
sys_GMF_initial_perf = simulate_and_evaluate(sys_GMF_initial)

# Simulação para Ganhos Novos

# Sistema original com ganho novo
sys_GMA_new = lti(A, B, C, 0)
H_new = A - B @ K_original_new
sys_GMF_new = lti(H_new, B, C, 0)

# Sistema observador com ganho novo
sys_GMA_new_perf = simulate_and_evaluate(sys_GMA_new)
sys_GMF_new_perf = simulate_and_evaluate(sys_GMF_new)

# Exibindo os resultados
print("Índices de Desempenho para Ganhos Iniciais (MP = 0.1, Ts = 1s):")
print(f"Sistema Original: Mp = {sys_GMA_initial_perf[0]:.2f}%, Tp = {sys_GMA_initial_perf[1]:.2f}s, Ts = {sys_GMA_initial_perf[2]:.2f}s, Erro de Regime = {sys_GMA_initial_perf[3]:.2f}")
print(f"Sistema com Controle: Mp = {sys_GMF_initial_perf[0]:.2f}%, Tp = {sys_GMF_initial_perf[1]:.2f}s, Ts = {sys_GMF_initial_perf[2]:.2f}s, Erro de Regime = {sys_GMF_initial_perf[3]:.2f}")

print("\nÍndices de Desempenho para Ganhos Novos (MP = 0.05, Ts = 0.5s):")
print(f"Sistema Original: Mp = {sys_GMA_new_perf[0]:.2f}%, Tp = {sys_GMA_new_perf[1]:.2f}s, Ts = {sys_GMA_new_perf[2]:.2f}s, Erro de Regime = {sys_GMA_new_perf[3]:.2f}")
print(f"Sistema com Controle: Mp = {sys_GMF_new_perf[0]:.2f}%, Tp = {sys_GMF_new_perf[1]:.2f}s, Ts = {sys_GMF_new_perf[2]:.2f}s, Erro de Regime = {sys_GMF_new_perf[3]:.2f}")


"""
Para o sistema de controle original, adote a transformação P diagonal, considerando o estado inicial 𝑥0𝑇=[0 1],
obtenha 𝐽𝑚𝑖𝑛 e os novos ganhos 𝐾=[𝑘 0,5𝑘], para os casos: 𝜆=0.2 e 𝜆=1.4; Observe que, com 𝐻=𝐴−𝐵𝐾, apenas soluções onde k>0 serão estáveis.
"""
#Letra E

# 1. Transformação Diagonal P
eigvals, eigvecs = np.linalg.eig(A)
P = eigvecs  # Matriz P de autovetores
P_inv = np.linalg.inv(P)
Lambda = np.diag(eigvals)  # Matriz diagonal dos autovalores
x0 = np.array([[0], [1]])

# 2. Cálculo de J_min
J_min = x0.T @ P @ x0

# 3. Cálculo dos Novos Ganhos K para λ = 0.2 e λ = 1.4
lambdas = [0.2, 1.4]
K_values = []

for lamb in lambdas:
    k = lamb
    K = np.array([[k, 0.5 * k]])
    H = A - B @ K
    # Garantindo que a solução é estável (k > 0 e autovalores negativos)
    eig_H = np.linalg.eigvals(H)
    if np.all(np.real(eig_H) < 0):
        K_values.append((lamb, K))

# Exibindo os resultados
print(f"J_min = {J_min[0,0]:.4f}")
for lamb, K in K_values:
    print(f"\nλ = {lamb}: K = {K[0]}")


"""
Considerando a entrada 𝑢=−𝐾𝑥, determine seu valor para o início da operação em cada caso estável de e);
"""
#Letra F
    
# Definição dos ganhos K obtidos anteriormente
# Exemplo: Para λ = 0.2 e λ = 1.4, usamos os ganhos obtidos no código anterior
K_02 = np.array([[0.2, 0.1]])  # Ganho para λ = 0.2
K_14 = np.array([[1.4, 0.7]])  # Ganho para λ = 1.4

# Cálculo da entrada u no início da operação para λ = 0.2
u_02 = -K_02 @ x0

# Cálculo da entrada u no início da operação para λ = 1.4
u_14 = -K_14 @ x0

# Exibindo os resultados
print(f"Valor de u para λ = 0.2: u = {u_02[0, 0]:.4f}")
print(f"Valor de u para λ = 1.4: u = {u_14[0, 0]:.4f}")


"""
Qual o fator de amortecimento e frequência natural em cada caso estável de e)?
Aplique excitação do tipo degrau nos dois casos e discuta o comportamento dinâmico do sistema em cada caso
"""
#Letra G

# Matriz H para cada caso
H_02 = A - B @ K_02
H_14 = A - B @ K_14

# Função para calcular o fator de amortecimento (zeta) e a frequência natural (omega_n)
def calculate_damping_and_natural_frequency(H):
    eigenvalues = np.linalg.eigvals(H)
    real_part = np.real(eigenvalues[0])
    imag_part = np.imag(eigenvalues[0])
    
    omega_n = np.sqrt(real_part**2 + imag_part**2)
    zeta = -real_part / omega_n
    
    return zeta, omega_n

# Cálculo de zeta e omega_n para cada caso
zeta_02, omega_n_02 = calculate_damping_and_natural_frequency(H_02)
zeta_14, omega_n_14 = calculate_damping_and_natural_frequency(H_14)

# Simulação da resposta ao degrau
t = np.linspace(0, 10, 1000)

# Sistema com λ = 0.2
sys_H_02 = lti(H_02, B, C, 0)
t_H_02, y_H_02 = step(sys_H_02, T=t)

# Sistema com λ = 1.4
sys_H_14 = lti(H_14, B, C, 0)
t_H_14, y_H_14 = step(sys_H_14, T=t)

# Plotando as respostas ao degrau
plt.figure(figsize=(10, 6))

plt.plot(t_H_02, y_H_02, label="λ = 0.2", color="blue")
plt.plot(t_H_14, y_H_14, label="λ = 1.4", color="red")

plt.xlabel('Tempo (s)')
plt.ylabel('Resposta ao Degrau')
plt.title('Resposta ao Degrau para os Casos λ = 0.2 e λ = 1.4')
plt.legend()
plt.grid(True)
plt.show()

# Exibindo os resultados
print(f"Para λ = 0.2: ζ = {zeta_02:.4f}, ω_n = {omega_n_02:.4f} rad/s")
print(f"Para λ = 1.4: ζ = {zeta_14:.4f}, ω_n = {omega_n_14:.4f} rad/s")