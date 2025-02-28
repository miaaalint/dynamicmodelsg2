# dynamicmodelsg2
Here we simulate a three-species predator-prey system consisting of plants  hares, and lynx. This model follows a Lotka-Volterra formulation with functional responses that regulate predation rates.

**Model Behavior:**
- Plants grow logistically and are consumed by hares.
- Lynx preys on hares.
- Lynx rely on hares as their primary food source and decline in population when hares are scarce.
- The system changes depending on parameter choices: stable oscillations, extinction scenarios, and chaotic dynamics.
 **Model Parameters:**
 - a1: (Plant consumption rate by hares): How efficiently hares feed on plants.
 - a2: (Hare consumption rate by lynx): Hpow effectively lynx preys on hares.
 - b1: (Plant growth limitation): How available plants are.
- b2: (Hare predation saturation): How hares avoid over-predation.
 - d1: (Hare mortality rate): Hare lifespan in the absence of predation.
 - d2: (Lynx mortality rate): Rate of lynx decline without hares.

 Scenario 1: The system exhibits stable oscillations with a periodicity of approximately 70 months
 Scenario 2: Setup represents a scenario where lynx go extinct after about one year, while plants and hares stabilize within 120 months.
 Scenario 3: The system shows chaotic behavior and the lynx population peaks twice within the 200-month period. What makes this behavior chaotic and not oscillating or random? Hint: Chaotic behavior is deterministic, sensitive to initial conditions, bounded and irregular.
 Within this assignment, MATLAB version R2023B, if laptop being used to run program does not have this update, we advise. The code provided does not require any additional packages, however throughout our code we have used aspects from the MATLAB toolboxes, such as ode23 which is a MATLAB function used to solve a system of differential equations. 
