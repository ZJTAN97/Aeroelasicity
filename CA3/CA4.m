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

% Unsteady Aerodynamics
sai_1 = 0.165;
sai_2 = 0.335;
epsilon_1 = 0.0455;
epsilon_2 = 0.30;

VELOCITY_CAP = 70;
US_freq = zeros(VELOCITY_CAP, 2);
US_damp = zeros(VELOCITY_CAP, 2);

US_flutter = 0;
US_stopper = 0;
US_idx = 0;


U = linspace(1, VELOCITY_CAP, VELOCITY_CAP);

tic
% Unsteady Aero using Wagners Function
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
disp("--Wagner--");
toc

% P-k Method
pk_omegas = zeros(VELOCITY_CAP, 2);
pk_damp = zeros(VELOCITY_CAP, 2);

pk_flutter = 0;
pk_stopper = 0;
pk_idx = 0;

tic
for iter = 1:length(U)
    
    error = 1;
    DOF = [1, 2]; % 1 for pitch; 2 for plunge
    omega_initials = [(Kh/MASS)^0.5, (Ka/I_a)^0.5];

    for item = 1:length(DOF)
        omega = omega_initials(item);
        while error > 0.01
            k = omega*b/U(iter);
            C_k = 1 - (0.165/(1-(0.0455i/k))) - (0.335/(1-(0.3i/k)));
            
            M = [ MASS, S;
                  S, I_a ];
              
            K = [ Kh, 0;
                  0, Ka ];
              
            D_k = [ -AIR_RHO*pi*(b^2),  AIR_RHO*pi*(b^2)*a;
                    AIR_RHO*pi*(b^2)*a,  -AIR_RHO*pi*(b^2)*(a^2 + b^2/8) ];
                
            E_k = [ 
                -2*pi*AIR_RHO*b*U(iter)*C_k,  -AIR_RHO*pi*(b^2)*U(iter) + (2*pi*AIR_RHO*b*U(iter)*C_k*(a-b/2));
                AIR_RHO*pi*(b^2)*U(iter) - (AIR_RHO*pi*b*U(iter)*(b-(2*a + b)*C_k)),  -AIR_RHO*pi*(b^2)*U(iter)*(c/4) + AIR_RHO*pi*b*U(iter)*(b-(2*a + b)*C_k)*(a - b/2) 
            ];
               
            F_k = [
                0,  -2*pi*AIR_RHO*b*(U(iter)^2)*C_k;
                0,  AIR_RHO*pi*(b^2)*(U(iter)^2) - AIR_RHO*pi*b*(U(iter)^2)*(b-(2*a+b)*C_k)
            ];

            A_k = [
                zeros(2),  eye(2);
                (inv(M - D_k))*(F_k - K),  (inv(M - D_k))*E_k
            ];

            [vectors, values] = eig(A_k);
            sorted_values = values(imag(values)~=0);
            
%             disp(sorted_values);
            
            if item == 1
                omega_new = abs(imag(sorted_values(2)));
            end
            if item == 2
               omega_new = abs(imag(sorted_values(2)));
            end         
            
            error = calculateError(omega, omega_new);
            omega = omega_new;
        end
        
        if item == 1
            pk_omegas(iter, 1) = abs(imag(sorted_values(1))) / (2*pi);
            pk_damp(iter, 1) = -real(sorted_values(1)) / abs(sorted_values(1));
        end
        if item == 2
           pk_omegas(iter, 2) = abs(imag(sorted_values(2))) / (2*pi);
           pk_damp(iter, 2) = -real(sorted_values(2)) / abs(sorted_values(2));
        end      
    end
    
    if(pk_damp(iter, 1) < 0 && pk_stopper == 0)
        pk_flutter = U(iter);
        pk_idx = iter;
        pk_stopper = pk_stopper + 1;
    end
    
end
disp("--Pk Method--");
toc


figure
title('Natural Frequency vs Velocity')
ylabel('Natural Frequency \omega (Hz)')
xlabel('Velocity, U(m/s)')
grid on
hold on 
plot(U, US_freq(:, 1), '--')
plot(U, US_freq(:, 2), '--')
plot(US_flutter, US_freq(US_idx, 1), 'b.', 'MarkerSize', 20)
text(US_flutter, US_freq(US_idx,1)+0.20, sprintf('Flutter Velocity (Wagner) = %.1f m/s', US_flutter))

plot(U, pk_omegas(:, 1))
plot(U, pk_omegas(:, 2))
plot(pk_flutter, pk_omegas(pk_idx, 1), 'r.', 'MarkerSize', 20)
text(pk_flutter-10, pk_omegas(US_idx,1)-0.2, sprintf('Flutter Velocity (P-k) = %.1f m/s', pk_flutter))

legend('Wagner Function', 'Wagner Function', 'Wagner Flutter Point', 'P-k Method', 'P-k Method', 'P-k Flutter Point')


figure
title('Damping Ratio vs Velocity')
ylabel('Damping Ratio \xi (Hz)')
xlabel('Velocity, U(m/s)')
grid on
hold on
plot(U, US_damp(:, 1), '--')
plot(U, US_damp(:, 2), '--')
plot(US_flutter, US_damp(US_idx, 1), 'b.', 'MarkerSize', 18)
text(US_flutter, US_damp(US_idx,1)+0.05, sprintf('Flutter Velocity (Wagner) = %.1f m/s', US_flutter))

plot(U, pk_damp(:, 1))
plot(U, pk_damp(:, 2))
plot(pk_flutter, pk_damp(pk_idx, 1), 'r.', 'MarkerSize', 18)
text(pk_flutter-10, pk_damp(pk_idx,1), sprintf('Flutter Velocity (P-k) = %.1f m/s', pk_flutter))

legend('Wagner Function', 'Wagner Function', 'Wagner Flutter Point', 'P-k Method', 'P-k Method', 'P-k Flutter Point')



%--Helper Functions--%

% Error Function
function [error] = calculateError(prev, next)
    error = (abs((next - prev))/prev)*100;
end
