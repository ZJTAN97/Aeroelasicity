c = 1;
a = c/4;
k1 = 500;
k2 = 2000;
m = 200;

K = [
  k1 0;
  0 k2
]

% mass of p range 1kg to 100kg
m_p = reshape(linspace(1,100,100), [100,1]);

% p distance from 0.75 to -0.25
p = reshape(linspace(0.75, -0.25, 100), [100,1]);

% Empty matrices of heave and pitch
heave = zeros(length(m_p), length(p));
pitch = zeros(length(m_p), length(p));

disp(m_p(2));

for i=1:length(m_p)
    for j=1:length(p)
        mass_term = m + m_p(i);
        S_term = (-1)*(m*a + m_p(i)*p(j));
        Ia_term = m*(((c^2)/12) + (a^2)) + (m_p(i) * (p(j)^2));
        
        M = [
           mass_term, S_term;
           S_term, Ia_term
        ];
    
        eigen_values = sqrt(eig(K, M));
        heave(i, j) = eigen_values(1);
        pitch(i, j) = eigen_values(2);
    end
end

disp(heave);
disp(pitch);