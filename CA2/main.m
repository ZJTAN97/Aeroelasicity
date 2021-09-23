clear all 
clc

% Equations for phi (Power Series) %
LENGTH = 44;
phi = @(y, i)(((y./LENGTH).^(i+1)).*(2+i-i.*(y./LENGTH)))/(i.*(i+1).*(i+2));
phi_dd = @(y, i) (1./(LENGTH.^3)).*(LENGTH-y).*((y./LENGTH).^(i-1));

% Equations for phi (Sine Series) %
phi_sine = @(y, i) ((-pi.*LENGTH.*i.*sin((i.*pi.*y)./LENGTH) + pi.*y.*i.*sin((i.*pi.*y)./LENGTH) + 2.*LENGTH.*i.*cos((i.*pi.*y)./LENGTH)) ./ ((pi.^3).*(i.^3).*LENGTH));
phi_dd_sine = @(y, i) ( ((LENGTH-y)./(LENGTH.^3)) .* (sin((i.*pi.*y)./LENGTH)));

error = 1;
n_conv = 4;
while error > 0.01
    n_conv = n_conv + 1;
    % to get first frequency so that can compare error.
    if error == 1
        [~, ~, ~, firstFreq, ~] = FormRitz(phi, phi_dd, 4, LENGTH);
        prevFreq = firstFreq(4);
    end
    [M, K, V, Freq, idx] = FormRitz(phi, phi_dd, n_conv, LENGTH);
    error = calculateError(prevFreq, Freq(4));
    prevFreq = Freq(4);
    %fprintf("Previous Frequency: %.5f\n", firstFreq(4));
    % fprintf("Current Frequency: %.5f\n", Freq(4));
    %fprintf("Current mode shape: %d\n", n_conv); 
    fprintf("Current error: %.3f\n" , error);
end


y_vals = linspace(0, LENGTH, 100);
eigen_vector = V(:, idx);

figure(1)
title('First 4 Associated Mode Shapes')
hold on
for i=1:4
    mode_shape = 0;
    for j=1:n_conv
        mode_shape = mode_shape + phi(y_vals, j).*eigen_vector(j,i);
    end
    [~, max_i] = max(abs(mode_shape));
    plot(y_vals/LENGTH, mode_shape/mode_shape(max_i));
end
yline(0, 'k--');
legend('Mode Shape 1', 'Mode Shape 2', 'Mode Shape 3', 'Mode Shape 4');
hold off

function [M, K, V, Freq, idx] = FormRitz(phi, phi_dd, n_terms, LENGTH)

MASS_ENGINE = 4100;
E = 2e12;
RHO = 300;

% Empty Matrices for Mass and Stifness
M = zeros(length(n_terms), length(n_terms));
K = zeros(length(n_terms), length(n_terms));

% Equation for second moment of Area %
I = @(y) 0.001.*(1 - (1./3).*(y./LENGTH) + (1./2).*((y./LENGTH).^2) - ((y./LENGTH).^3));
% Equation for Area %
A = @(y) 0.1.*(3 - (y./LENGTH) + 2.*((y./LENGTH).^2));

% Equations for Mass Matrix %
wing_contrib = @(y, i, j) RHO.*A(y).*phi(y, i).*phi(y, j);
engine1_contrib = @(i, j) 0.5.* MASS_ENGINE.*phi(LENGTH/5, i).*phi(LENGTH/5, j);
engine2_contrib = @(i, j) 0.5.* MASS_ENGINE.*phi((2*LENGTH/5), i).*phi((2*LENGTH)/5, j);
engine3_contrib = @(i, j) 0.5.* MASS_ENGINE.*phi((3*LENGTH/5), i).*phi((3*LENGTH)/5, j);

% Equation for Stiffness Matrix %
stiffness_contrib = @(y, i, j) E.*I(y).*phi_dd(y, i).*phi_dd(y, j);

for i = 1:n_terms
    for j = 1:n_terms
        wing = integral(@(y) wing_contrib(y, i, j), 0, LENGTH);
        engine1 = engine1_contrib(i, j);
        engine2 = engine2_contrib(i, j);
        engine3 = engine3_contrib(i, j);
        stiffness = integral(@(y) stiffness_contrib(y, i, j), 0, LENGTH);
        M(i, j) = wing + engine1 + engine2 + engine3;
        K(i, j) = stiffness;
    end
end

% D --> Eigenv Value
% V --> Eigen Vector: 
[V, D] = eig(K, M);
[Freq, idx] = sort(diag(sqrt(D)/(2*pi)));
end


function [error] = calculateError(prev, next)
% fprintf("Original: %.4f ", prev);
% fprintf("Next: %.4f\n", next);
error = (abs((next - prev))/prev)*100;

end
