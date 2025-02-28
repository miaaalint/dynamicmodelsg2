% Here we simulate a three-species predator-prey system consisting of plants  hares, and lynx. 
% This model follows a Lotka-Volterra formulation with functional responses that regulate predation rates.
% 
% **Model Behavior:**
% - Plants grow logistically and are consumed by hares.
% - Lynx preys on hares.
% - Lynx rely on hares as their primary food source and decline in population when hares are scarce.
% - The system changes depending on parameter choices: stable oscillations, extinction scenarios, and chaotic dynamics.
%
% **Model Parameters:**
% - a1: (Plant consumption rate by hares): How efficiently hares feed on plants.
% - a2: (Hare consumption rate by lynx): Hpow effectively lynx preys on hares.
% - b1: (Plant growth limitation): How available plants are.
% - b2: (Hare predation saturation): How hares avoid over-predation.
% - d1: (Hare mortality rate): Hare lifespan in the absence of predation.
% - d2: (Lynx mortality rate): Rate of lynx decline without hares.
%
% Scenario 1: The system exhibits stable oscillations with a periodicity of approximately 70 months
% Scenario 2: Setup represents a scenario where lynx go extinct after about one year, while plants and hares stabilize within 120 months.
% Scenario 3: The system shows chaotic behavior and the lynx population peaks twice within the 200-month period. What makes this behavior chaotic and not oscillating or random? Hint: Chaotic behavior is deterministic, sensitive to initial conditions, bounded and irregular.
%Within this assignment, MATLAB version R2023B, if laptop being used to run program does not have this update, we advise. The code provided does not require any additional packages, however throughout our code we have used aspects from the MATLAB toolboxes, such as ode23 which is a MATLAB function used to solve a system of differential equations. 
% Simulation time parameters
clear; clc;
t0 = 0;      %time initial              
tfinal = 200;  % Run long enough to observe stable oscillations
tspan = [t0 tfinal];     

% Initial populations: [Plants, Hares, Lynx] (Away from equilibrium)
y0 = [0.9, 0.1, 8];  

% Corrected Parameters for Stable 70-Month Oscillations
a1 = 5;  % Increased plant growth rate
b1 = 3;  % Slightly reduced plant-hare interaction strength
a2 = 0.1;  % Increased hare reproduction rate
b2 = 1;  % Slightly higher hare-lynx interaction strength
d1 = 0.4;  % Reduced natural death rate of hares
d2 = 0.01; % Reduced natural death rate of lynx

% Solve the ODE system
[t, y] = ode23(@(t,y) plants_hare_lynx(t, y, a1, a2, b1, b2, d1, d2), tspan, y0);

% Plot population dynamics
figure;
plot(t, y(:,1), 'g', 'LineWidth', 2); hold on;
plot(t, y(:,2), 'b', 'LineWidth', 2);
plot(t, y(:,3), 'r', 'LineWidth', 2);
xlabel('Time [months]'); ylabel('Population');
legend('Plants', 'Hares', 'Lynx');
title('Stable 70-Month Oscillations in Plants-Hares-Lynx System');
grid on;

% Compute oscillation period using Fourier Transform
Y = fft(y(:,2)); % FFT of hare population
freq = (0:length(Y)-1) / (t(end) - t(1)); % Frequency in cycles per month
[~, idx] = max(abs(Y(2:end))); % Find dominant frequency
dominant_period = 1 / freq(idx+1); % Convert frequency to period

fprintf('Estimated Oscillation Period: %.2f months\n', dominant_period);

% 3D State-Space Plot (Plants vs. Hares vs. Lynx)
figure;
plot3(y(:,1), y(:,2), y(:,3), 'b-', 'LineWidth', 2);
xlabel('Plant Population');
ylabel('Hare Population');
zlabel('Lynx Population');
title('3D State-Space Plot (Plants, Hares, Lynx)');
grid on;



%b2 represents the lynx's predation rate on the hares. A low b2 basically
%represents the simulatanous relationship between the increase of the hare
%population and lynxs themselves becoming better predators. B2 in this
%system ultimately moderates the predator, prey interaction and controls
%the stability of it.

%% Simulation settings

t0 = 0; % start time of simulation (months)
tfinal = 200; % total simulation time
tspan = [t0 tfinal];
options = odeset('RelTol',1e-5,'AbsTol',1e-10, 'InitialStep', 1e-3, 'MaxStep', 2, 'NonNegative', [1 2 3]);

%% Scenario 2: Lynx Extinction (~1 year)
y0 = [0.9, 0.1, 8];
a1 = 3.7; a2 = 0.01;
b1 = 4; b2 = 0.1;
d1 = 0.6; d2 = 0.07;

[t,y] = ode23t(@(t,y) plants_hare_lynx(t,y,a1,a2,b1,b2,d1,d2), tspan, y0, options);
figure;
plot(t, y);
title('Lynx Extinction (~1 year)');
xlabel('Time (months)');
ylabel('Population Size');
legend('Plants','Hares','Lynx');
%% State Space Plot for Lynx Extinction Scenario
figure; plot(y(:,1), y(:,2));
title('State Space: Plants vs Hares after Lynx Extinction');
xlabel('Plant Population'); ylabel('Hare Population');

%scenario 3


% Simulation time parameters
t0 = 0;                   
tfinal = 1000;  % Run long enough to observe stable oscillations
tspan = [t0 tfinal];     

% Initial populations: [Plants, Hares, Lynx] (Away from equilibrium)
y0 = [0.9, 0.1, 8];  

% Corrected Parameters for Stable 70-Month Oscillations
a1 = 5;  % Increased plant growth rate
b1 = 3;  % Slightly reduced plant-hare interaction strength
a2 = 0.1;  % Increased hare reproduction rate
b2 = 2;  % Slightly higher hare-lynx interaction strength
d1 = 0.4;  % Reduced natural death rate of hares
d2 = 0.01; % Reduced natural death rate of lynx

% Solve the ODE system
[t, y] = ode23(@(t,y) plants_hare_lynx(t, y, a1, a2, b1, b2, d1, d2), tspan, y0);

% Plot population dynamics
figure;
plot(t, y(:,1), 'g', 'LineWidth', 2); hold on;
plot(t, y(:,2), 'b', 'LineWidth', 2);
plot(t, y(:,3), 'r', 'LineWidth', 2);
xlabel('Time [months]'); ylabel('Population');
legend('Plants', 'Hares', 'Lynx');
title('Stable 70-Month Oscillations in Plants-Hares-Lynx System');
grid on;

%Prints the oscillation period
fprintf('Estimated Oscillation Period: %.2f months\n', dominant_period);

% 3D State-Space Plot (Plants vs. Hares vs. Lynx)
figure;
plot3(y(:,1), y(:,2), y(:,3), 'b-', 'LineWidth', 2);
xlabel('Plant Population');
ylabel('Hare Population');
zlabel('Lynx Population');
title('3D State-Space Plot (Plants, Hares, Lynx)');
grid on;

% Deterministic Trajectory Plot (Hares vs. Lynx)
figure;
plot(y(:,2), y(:,3), 'r-', 'LineWidth', 2);
xlabel('Hare Population');
ylabel('Lynx Population');
title('Deterministic Trajectory: Hares vs. Lynx');
grid on;

%Graphing Poncare plot for hare population
figure;
scatter(y(1:end-1, 2), y(2:end, 2), 'go');  % 'go' for green circles
xlabel('Hare Population at Time t(i)');
ylabel('Hare Population at Time t(i+1)');
title('Scatter Plot: Hare Population at t(i) vs. t(i+1)');
grid on;

% Differential Equations Model
function dydt = plants_hare_lynx(~, y, a1, a2, b1, b2, d1, d2)
    dydt = zeros(3,1);
    dydt(1) = y(1)*(1 - y(1)) - ((a1*y(1)) / (1 + b1*y(1))) * y(2); % Plants
    dydt(2) = ((a1*y(1)) / (1 + b1*y(1))) * y(2) - d1*y(2) - ((a2*y(2)) / (1 + b2*y(2))) * y(3); % Hares
    dydt(3) = ((a2*y(2)) / (1 + b2*y(2))) * y(3) - d2*y(3); % Lynx
end

%deterministic, sensitive to initial conditions, clear shape 
%when looking at the Poncare plot, there is a clear shape, thus revealing
%the chaotic nature of the system and classifying it as determinsitic. 
%Sensitive to initial conditions - when changing the initial population of the lynx from 8 to 7.8, there was a noticeable leftward shift. Such high sensitivity to slight changes in the initial conditions is characteristic for chaotic systems and is a defining behaviour 
