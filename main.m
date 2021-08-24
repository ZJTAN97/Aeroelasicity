c = 1;
a = c/4;
k1 = 500;
k2 = 2000;
m = 200;

K = [
  k1 0;
  0 k2
];

% mass of p range 1kg to 100kg
m_p = reshape(linspace(1,100,100), [100,1]);

% p distance from 0.75 to -0.25
p = reshape(linspace(0.75, -0.25, 100), [100,1]);

% Empty matrices of heave and pitch
heave = zeros(length(m_p), length(p));
pitch = zeros(length(m_p), length(p));

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
        % divide by 2pi to convert to get Hertz (7:45 of week 2 revision lecture)
        heave(i, j) = eigen_values(1) / (2*pi);
        pitch(i, j) = eigen_values(2) / (2*pi);
    end
end

plotter(p, m_p, heave, 'Heave')
plotter(p, m_p, pitch, 'Pitch')


function plotter(p, m_p, freq, name)
    figure
    mesh(p, m_p, freq)
    title(sprintf('Natural Frequences for %s', name))
    xlabel('p (m)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b')
    ylabel('Mass of p (kg)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b')
    zlabel('Freq (Hz)', 'FontSize', 12, 'FontWeight', 'bold', 'Color', 'b')
end





