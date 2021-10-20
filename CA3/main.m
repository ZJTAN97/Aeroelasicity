clc
clear all

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

U = linspace(1, 70, 7000);
QS_freq = zeros(70, 2);
QS_damp = zeros(70, 2);
I = eye(2);
ZEROS = zeros(2);

QS_flutter = 0;
stopper = 0;
idx = 0;

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
    
    if(QS_damp(i, 2) < 0 && stopper == 0)
        QS_flutter = U(i);
        idx = i;
        stopper = stopper + 1;
    end

end


% Unsteady Aerodynamics
sai_1 = 0.165;
sai_2 = 0.335;
epsilon_1 = 0.0455;
epsilon_2 = 0.30;

US_freq = zeros(70, 2);
US_damp = zeros(70, 2);

US_flutter = 0;
US_stopper = 0;
US_idx = 0;

for j = 1:length(U)
    
    phi_0 = 0.5; % slide 29
    phi_d_0 = (epsilon_1 * U(j) * sai_1)/b + (epsilon_2 * U(j) * sai_2)/b;
    
    M = [
        MASS + AIR_RHO*pi*(b^2), S - AIR_RHO*pi*(b^2)*(Xf-b);
        S - AIR_RHO*pi*(b^2)*(Xf-b),  I_a + AIR_RHO*pi*(b^2)*((Xf-b)^2 + b^2/8)
    ];

    C = AIR_RHO*pi*U(j)*c*[
        phi_0,  c/4 + phi_0 * (3*c/4 - Xf);
        -e*c*phi_0,  (3*c/4 - Xf)*(c/4 - e*c*phi_0);
    ];

    K = [
        Kh + AIR_RHO*pi*U(j)*c*phi_d_0,  AIR_RHO*pi*U(j)*c*(U(j)*phi_0 + (3*c/4 - Xf)*phi_d_0);
        -AIR_RHO*pi*U(j)*e*(c^2)*phi_d_0,  Ka - AIR_RHO*pi*U(j)*e*(c^2)*(U(j)*phi_0 + (3*c/4 - Xf) * phi_d_0)
    ];
    
    W = 2*AIR_RHO*pi*(U(j)^3)*[
        -sai_1*(epsilon_1^2)/b, -sai_2*(epsilon_2^2)/b,  sai_1*epsilon_1*(1-epsilon_1*(1-2*e)),  sai_2*epsilon_2*(1-epsilon_2*(1-2*e));
        e*c*sai_1*(epsilon_1^2)/b, e*c*sai_2*(epsilon_2^2)/b,  -e*c*sai_1*epsilon_1*(1-epsilon_1*(1-2*e)), -e*c*sai_2*epsilon_2*(1-epsilon_2*(1-2*e))
    ];
        
    B = [1 0; 1 0; 0 1; 0 1];
    
    G = [
        -epsilon_1*U(j)/b 0 0 0;
        0 -epsilon_2*U(j)/b 0 0;
        0 0 -epsilon_1*U(j)/b 0;
        0 0 0 -epsilon_2*U(j)/b;
    ];
    
    Q = [
        (-inv(M)*C), (-inv(M)*K), (-inv(M)*W);
        eye(2), zeros(2), zeros(2,4);
        zeros(4,2), B, G  
    ];

    [vector, value] = eig(Q);
    sorted_value = value(imag(value)~=0);
        
    US_freq(j, 1) = abs(sorted_value(1)) / (2*pi); 
    US_freq(j, 2) = abs(sorted_value(3)) / (2*pi);
    
    US_damp(j, 1) = -real(sorted_value(1)) / abs(sorted_value(1));
    US_damp(j, 2) = -real(sorted_value(3)) / abs(sorted_value(3));
    
    if(US_damp(j, 1) < 0 && US_stopper == 0)
        US_flutter = U(j);
        US_idx = j;
        US_stopper = US_stopper + 1;
    end
    
end


figure
title('Natural Frequency vs Velocity')
ylabel('Natural Frequency \omega (Hz)')
xlabel('Velocity, U(m/s)')
grid on
hold on 
plot(U, QS_freq(:, 1), '--')
plot(U, QS_freq(:, 2), '--')
plot(QS_flutter, QS_freq(idx, 2), 'r.', 'MarkerSize', 18)
plot(U, US_freq(:, 1))
plot(U, US_freq(:, 2))
plot(US_flutter, US_freq(US_idx, 1), 'b.', 'MarkerSize', 18)
legend('QS', 'QS', 'QS Flutter Point', 'US', 'US', 'US Flutter Point')
text(QS_flutter-20, QS_freq(idx,2)-0.15, sprintf('Velocity = %.2f m/s', QS_flutter))
text(US_flutter-20, US_freq(US_idx,1)-0.15, sprintf('Velocity = %.2f m/s', US_flutter))


figure
title('Damping Ratio vs Speed')
ylabel('Damping Ratio \xi (Hz)')
xlabel('Velocity, U(m/s)')
grid on
hold on
title('Damping vs Speed')
plot(U, QS_damp(:, 1), '--')
plot(U, QS_damp(:, 2), '--')
plot(QS_flutter, QS_damp(idx, 2), 'r.', 'MarkerSize', 18)
plot(U, US_damp(:, 1))
plot(U, US_damp(:, 2))
plot(US_flutter, US_damp(US_idx, 1), 'b.', 'MarkerSize', 18)
legend('QS', 'QS', 'QS Flutter Point', 'US', 'US', 'US Flutter Point')
text(QS_flutter, QS_damp(idx,2)-0.10, sprintf('Velocity = %.2f m/s', QS_flutter))
text(US_flutter, US_damp(US_idx,1)+0.05, sprintf('Velocity = %.2f m/s', US_flutter))


