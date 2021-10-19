clear all
clc

% Variables
MATERIAL_RHO = 2000;
AIR_RHO = 1.225;
c = 0.5;
b = c/2;
t = 0.02;
Xf = 0.45 * c;
LENGTH = 2.0;
MASS = MATERIAL_RHO * t * c * LENGTH;
Wh = 2 * 2 * pi;
Wa = 6 * 2 * pi;

a = Xf - b;
e = (Xf/c) - 0.25;

S = -MASS * a;
I_a = ((1/12)*MASS*c^2)+(MASS*a^2);

Kh = MASS*(Wh^2);
Ka = I_a*(Wa^2);

% Forming Matrix A, B, C, D , E, F
A = [ MASS S;
      S I_a ];

B = pi * (b^2) * [
      1 (b - Xf);
      (b - Xf) (((b-Xf)^2)+((b^2)/8)) ];

C = zeros(2);

D = pi * c * [ 
      1 (((3*c/4)-Xf)+c/4);
      (-e*c) ((b-Xf)^2+((3*c/4)-Xf)*(c/4)) ];

E = [ Kh 0;
      0 Ka ];

F = pi* c * [ 
      0 1;
      0 (-e*c) ];

U = linspace(1, 70);
QS_freq = zeros(70, 2);
QS_damp = zeros(70, 2);
I = eye(2);
ZEROS = zeros(2);

for i=1:length(U)
  
    M = A + AIR_RHO * B;
    H = C + AIR_RHO * U(i) * D;
    K = E + AIR_RHO * (U(i)^2) * F;

    
    Q = [ I ZEROS; ZEROS M ] \ [ ZEROS I; -K -H ];
    [vector, value] = eig(Q);
    sorted_value = sort(diag(value));
    
    QS_freq(i, 1) = abs(sorted_value(1)) / (2*pi); 
    QS_freq(i, 2) = abs(sorted_value(3)) / (2*pi);
    
    QS_damp(i, 1) = -real(sorted_value(1)) / abs(sorted_value(1));
    QS_damp(i, 2) = -real(sorted_value(3)) / abs(sorted_value(3));
    
    
end


figure
title('Frequency vs Speed')
hold on 
grid on
plot(U, QS_freq(:, 1), '--')
plot(U, QS_freq(:, 2), '--')


figure
hold on
grid on
title('Damping vs Speed')
plot(U, QS_damp(:, 1), '--')
plot(U, QS_damp(:, 2), '--')






