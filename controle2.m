clear;
clc;

%%
syms s ka1 ka2 ka3 ka4

Ka = [ka1 ka2 ka3 ka4];

k1 = 1.0; k2 = 1.0;
b1 = 1.0; b2 = 1.0; 
m1 = 1.0; m2 = 1.0;

A = [
        0.0           1.0        0.0        0.0; 
    -(k1+k2)/m1     -b1/m1      k2/m1       0.0;
        0.0           0.0        0.0        1.0; 
      k2/m2         0.0       -k2/m2      -b2/m2
    ];

B = [
            0.0; 
           k1/m1; 
            0.0; 
            0.0
    ];

C = [1.0     0.0     0.0     0.0]

D = 0;

A
B
C
D

% Matriz de controlabilidade
Ctrb = ctrb(A, B);
disp('Matriz de controlabilidade:');
disp(Ctrb);

% Verifica se o sistema é controlável
disp('Posto: ');
disp(rank(Ctrb));

% Equacao caracteristica
equacao_caracteristica = det(s*eye(size(A)) - A);
raizes_ec = solve(equacao_caracteristica, s);

% Polos do sistema
disp('E.C: ');
disp(equacao_caracteristica);
disp('Raizes da E.C: ');
disp(double(raizes_ec));

% Matriz da Malha Fechada
Amf = A - B*Ka;

% E.C da Malha Fechada
equacao_caracteristica_mf = det(s*eye(size(Amf)) - Amf);
raizes_ec_mf = solve(equacao_caracteristica_mf, s);
disp('det(A - BK):');
disp(equacao_caracteristica_mf);
%%
%ESCOLHENDO OS POLOS COMO -3, -4, -20 E -100
At = [A zeros(4,1); -C 0];
Bt = [B;0];
Et = [zeros(4,1);1];
Ct = [C 0];

P = [-3 -4 -20 -50 -100];
Kt = acker(At, Bt, P);
Ka = Kt(1:4);
Ki = -Kt(5);

Ka
Ki
Amf = At - Bt*Kt;
Amf
%%
t = 0:0.1:20; % Tempo de 0 a 10 unidades de tempo, com intervalos de 1 unidade

% Defina a entrada do sistema (por exemplo, uma entrada de degrau)
u = ones(size(t)); % Degrau unitário para cada ponto de tempo

yMA = lsim(A, B, C, D, u, t);

% Simule a resposta do sistema usando lsim
yMF = lsim(Amf, Et, Ct, 0, u, t);

% Plote a resposta do sistema
figure(1);
plot(t, yMA);
xlabel('Tempo (unidades de tempo)');
ylabel('Saída do Sistema');
title('Resposta do Sistema em Malha Aberta');
grid on;

% Plote a resposta do sistema
figure(2);
plot(t, yMF);
xlabel('Tempo (unidades de tempo)');
ylabel('Saída do Sistema');
title('Resposta do Sistema em Malha Fechada');
grid on;


%%
% Escolhendo os polos desejados para o observador
Po = [-1, -2, -40, -50]; 

% Calculando matriz L de observabilidade usando o método de Ackermann
L = acker(A', C', Po);
L = L';

% A matriz L agora contém os ganhos do observador
disp('Matriz L de Observabilidade:');
disp(L);