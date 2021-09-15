MASS_ENGINE = 4100;
E = 2e12;
RHO = 300;
LENGTH = 44;

% Number of Modes to Keep 
n = 4;

% Empty Matrices for Mass and Stifness
M = zeros(length(n), length(n));
K = zeros(length(n), length(n));

syms y;
% Equation for second moment of Area %
I(y) = 0.001 * (1 - (1/3)*(y/l) + (1/2)*((y/l)^2) - ((y/l)^3));


% Equation for Area %
A(y) = 0.1 * (3 - (y/LENGTH) + 2*((y/LENGTH)^2));


% Equation for phi %
syms y;
syms l;
syms i;
phi(y,l,i) = (((y/l)^(i+1)) * (2+i-i*(y/l))) / (i*(i+1)(i+2));
phi_dd = ();

% Equation for Mass Matrix %
wing_portion = rho * A(y) * phi(y,LENGTH,i) * phi(y,LENGTH,j);
engine1_portion = 0.5 * MASS_ENGINE * phi(y, LENGTH/5, i) * phi(y, LENGTh/5, j);
engine2_portion = 0.5 * MASS_ENGINE * phi(y, (2*LENGTH)/5, i) * phi(y, (2*LENGTH)/5, j);
engine3_portion = 0.5 * MASS_ENGINE * phi(y, (3*LENGTH)/5, i) * phi(y, (3*LENGTH)/5, j);


% Equation for Stiffness Matrix %
K = 0.5 * E * I(y) *


for i = 1:n
    for j = 1:n
         disp(phi(10, 5, j));
    end
end

