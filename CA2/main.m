MASS_ENGINE = 4100;
E = 2e12;
RHO = 300;
LENGTH = 44;

% Number of Modes to Keep 
n = 4;

% Empty Matrices for Mass and Stifness
M = zeros(length(n), length(n));
K = zeros(length(n), length(n));


% Equation for second moment of Area %
I = @(y) 0.001.*(1 - (1./3).*(y./LENGTH) + (1./2).*((y./LENGTH).^2) - ((y./LENGTH).^3));

% Equation for Area %
A = @(y) 0.1.*(3 - (y./LENGTH) + 2.*((y./LENGTH).^2));


% Equations for phi %
phi = @(y, i)((y./LENGTH).^(i+1).*(2+i-i.*(y./LENGTH)))/(i.*(i+1).*(i+2));
phi_dd = @(y, i) (1./(LENGTH.^3)).*(LENGTH-y).*((y./LENGTH).^(i-1));


% Equations for Mass Matrix %
wing_contrib = @(y, i, j) RHO.*A(y).*phi(y, i).*phi(y, j);
engine1_contrib = @(i, j) 0.5.* MASS_ENGINE.*phi(LENGTH/5, i).*phi(LENGTH/5, j);
engine2_contrib = @(i, j) 0.5.* MASS_ENGINE.*phi((2*LENGTH/5), i).*phi((2*LENGTH)/5, j);
engine3_contrib = @(i, j) 0.5.* MASS_ENGINE.*phi((3*LENGTH/5), i).*phi((3*LENGTH)/5, j);

% Equation for Stiffness Matrix %
stiffness_contrib = @(y, i, j) E.*I(y).*phi_dd(y, i).*phi_dd(y, j);

for i = 1:n
    for j = 1:n
        wing = integral(@(y) wing_contrib(y, i, j), 0, 44);
        engine1 = engine1_contrib(i, j);
        engine2 = engine2_contrib(i, j);
        engine3 = engine3_contrib(i, j);
        stiffness = integral(@(y) stiffness_contrib(y, i, j), 0, 44);
        M(i, j) = wing + engine1 + engine2 + engine3;
        K(i, j) = stiffness;
    end
end


disp(M);
disp(K);





