% SNE - Numerical Simulation in Engineering (MCA - University of Minho (Braga))
% Authors-> Diogo Silva and Tomás Pereira 

% Solving Boundary Value Problem (BVP)
% The equation to solve is defined here

% This problem is to find the temperature distribution along a bare wire
% INFO: RANGE -> 0 m to Length/2 and temperature distribution is symmetrical at L/2
% This problem has two boundaries -> One of the boundaries is a Neumman
% Boundary and the other one is a Drichlet boundary -> This makes the
% problem mixed boundary

% We are solving the equation 11.28 in book Numerical Methods for Engineers
% and Scientists

function solve_heat()

    % Define constants
    k = 72.0; % Thermal conductivity
    h = 2000.0; % Convective heat
    e = 0.1; % Radiative emissivity
    sigma = 5.67 * 10^(-8); % Stefan-Boltzmann constant
    I = 2; % Current on the wire
    p = 32 * 10^(-8); % Electrical resistivity
    Tamb = 300; % Ambient temperature
    D = 7.62 * 10^(-5); % Wire diameter
    L = 4*10^(-3); % Length of copper wire
    
    % Constants calculated from parameters (you can update as needed)
    C1 = 1.45 * 10^6; % m⁻²K⁻¹
    C2 = 4.13 * 10^(-6); % m⁻²K⁻³
    C3 = 8.54 * 10^8; % m⁻²K
    
    % Define the domain and initial guess
    T0 = Tamb; % Initial temperature equal to ambient temperature
    domain = bvpinit(linspace(0, L/2, 50), [T0, 0]);
    
    % Solve the BVP
    sol = bvp4c(@bvsolve, @bvbound, domain);
    
    % Extract the solution
    x = sol.x;
    y = sol.y;
    
    % Plot results
    plot(x, y(1,:), '-o')
    xlabel('x (m)')
    ylabel('Temperature (K)')
    title('Temperature Distribution Along the Wire')

    maxTemp = max(y(1,:));
    fprintf('Maximum temperature along the wire: %.2f K\n', maxTemp);
    
    function dxdy = bvsolve(x, y)
        dxdy = [y(2); C1*(y(1) - Tamb) + C2*(y(1)^4 - Tamb^4) - C3];
    end

    function res = bvbound(ya, yb)
        res = [ya(1) - Tamb; yb(2)];  % Dirichlet at left, Neumann at right
    end

end

solve_heat();

