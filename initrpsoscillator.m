%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Shalin Shah
% Email: shalin.shah@duke.edu
% Last edited: 03/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ideal CRN network
crnidealrpsoscillator(0.5, 0.4, 1.0, 50000, 1);
% DNA polymerase-based CRN network
crnrpsoscillator(0.5, 0.4, 1.0, 50000, 1);

function crnidealrpsoscillator(concA, concB, concC, stopTime, figNo)
    % enter A, B, C 
    %
    % stopTime is the total simulation time
    %
    % figNo is the display figure number
    
    % rate of ideal reactions 
    rate = 1e-3;
    
    % Create the amplification model
    model = sbiomodel('Ideal rps CRN network');

    r = cell(1, 3);
    k = cell(1, length(r));
    p = cell(1, length(r));
         
    % Add reaction set to the solver
    r{1} = addreaction(model,'A + B -> B + B');
    r{2} = addreaction(model,'B + C -> C + C');
    r{3} = addreaction(model,'C + A -> A + A');
        
    % Set the Kinetic Law for Reactions.
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
    end
    
    p{1} = addparameter(k{1}, 'c1', 'Value', rate);
    p{2} = addparameter(k{2}, 'c2', 'Value', rate);
    p{3} = addparameter(k{3}, 'c3', 'Value', rate);
    
    % Set the Kinetic Law for Reactions.
    for i = 1:length(r)
        k{i}.ParameterVariableNames = {strcat('c',num2str(i))};
    end
    
    % Set initial amounts for species 
    r{1}.Reactants(1).InitialAmount = concA;        % A
    r{2}.Reactants(1).InitialAmount = concB;        % B
    r{3}.Reactants(1).InitialAmount = concC;        % C

    % Display the Completed Model Objects
    model

    % Display the Reaction Objects
    model.Reactions

    % Display the Species Objects
    model.Species

    % Simulate ODE and plot for totalTime
    cs = getconfigset(model,'active');
    % allowed error rate tolerance value
    absTol = 1e-5;
    relTol = 1e-5;
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);

    % plot the phase of specie 
    figure(figNo);
    hold on; box on;
    plot3(X(:,1)./1e-1, X(:,2)./1e-1, X(:, 3)./1e-1, ':', 'LineWidth', 2.0);
    zlabel('Z'); ylabel('Y'); xlabel('X');
    set(gca, 'LineWidth', 2.0); 
    
    % plot the specie concentration over time
    figure(figNo+1);
    box on; hold on;
    plot(t./3600, X(:,1)./1e-1, 'b:', 'LineWidth', 2.0);
    plot(t./3600, X(:,2)./1e-1, 'r:', 'LineWidth', 2.0);
    plot(t./3600, X(:,3)./1e-1, 'y:', 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
end

function crnrpsoscillator(concA, concB, concC, stopTime, figNo)
    % enter A, B and C concentration (in 10 nM) assuming all the gates are
    % in excess concentration
    %
    % stopTime is the total simulation time (s)
    %
    % figNo is the display figure number
    
    k = 1.0;
    alpha = 1e-3; beta = 1e-8; gamma = 2;
    % slow down all reactions by multiplying with a scaling factor 1e-3
    rateConst = k * (alpha);
    % fastest reaction are assumed to occur at this rate
    infRate = 1e3 * (k);
    % inf concentration of gates will be 1 mM
    infConc = 1e5 * (beta);
        
    % Create the amplification model
    model = sbiomodel('DNA polymerase-based CRN network');

    r = cell(1, 12);
    k = cell(1, length(r));
    p = cell(1, length(r));
    
    % Add reaction set A+B -> 2B to the solver
    r{1} = addreaction(model,'A + Gb -> I1');
    r{2} = addreaction(model,'I1 + Gbo -> B + B');
    r{3} = addreaction(model,'B + BfB1 -> Gb');
    r{4} = addreaction(model,'Gb + BfB2 -> B');

    % Add reaction set C+A -> 2A to the solver
    r{5} = addreaction(model,'C + Ga -> I2');
    r{6} = addreaction(model,'I2 + Gao -> A + A');
    r{7} = addreaction(model,'A + BfA1 -> Ga');
    r{8} = addreaction(model,'Ga + BfA2 -> A');

    % Add reaction set B+C -> 2C to the solver
    r{9} = addreaction(model,'B + Gc -> I3');
    r{10} = addreaction(model,'I3 + Gco -> C + C');
    r{11} = addreaction(model,'C + BfC1 -> Gc');
    r{12} = addreaction(model,'Gc + BfC2 -> C');
    
    % Set the Kinetic Law for Reactions
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction');
    end
    
    % set the rate constants for each reaction
    q1 = rateConst*(gamma)*(1/beta); 
    q2 = infRate; q3 = infRate; q4 = infRate;
    p{1} = addparameter(k{1}, 'c1', 'Value', q1);
    p{2} = addparameter(k{2}, 'c2', 'Value', q2);
    p{3} = addparameter(k{3}, 'c3', 'Value', q3);
    p{4} = addparameter(k{4}, 'c4', 'Value', q4);
    
    q5 = rateConst*(gamma)*(1/beta); 
    q6 = infRate; q7 = infRate; q8 = infRate;
    p{5} = addparameter(k{5}, 'c5', 'Value', q5);
    p{6} = addparameter(k{6}, 'c6', 'Value', q6);
    p{7} = addparameter(k{7}, 'c7', 'Value', q7);
    p{8} = addparameter(k{8}, 'c8', 'Value', q8);
    
    q9 = rateConst*(gamma)*(1/beta); 
    q10 = infRate; q11 = infRate; q12 = infRate;
    p{9} = addparameter(k{9}, 'c9', 'Value', q9);
    p{10} = addparameter(k{10}, 'c10', 'Value', q10);
    p{11} = addparameter(k{11}, 'c11', 'Value', q11);
    p{12} = addparameter(k{12}, 'c12', 'Value', q12);

    % Set the Kinetic Law for Reactions.
    for i = 1:length(r)
        k{i}.ParameterVariableNames = {strcat('c',num2str(i))};
    end

    % Set initial amounts for species 
    r{1}.Reactants(1).InitialAmount = concA*beta;      % A
    r{1}.Reactants(2).InitialAmount = concB*beta;      % Gb
    r{5}.Reactants(1).InitialAmount = concC*beta;      % C
    r{5}.Reactants(2).InitialAmount = concA*beta;      % Ga
    r{9}.Reactants(1).InitialAmount = concB*beta;      % B
    r{9}.Reactants(2).InitialAmount = concC*beta;      % Gc
    

    % set gate concentrations in excess
    r{2}.Reactants(2).InitialAmount = infConc;              % Gbo
    r{3}.Reactants(2).InitialAmount = infConc;              % BfB1
    r{4}.Reactants(2).InitialAmount = infConc;              % BfB2
    
    r{6}.Reactants(2).InitialAmount = infConc;              % Gao
    r{7}.Reactants(2).InitialAmount = infConc;              % BfA1
    r{8}.Reactants(2).InitialAmount = infConc;              % BfA2
    
    r{10}.Reactants(2).InitialAmount = infConc;              % Gco
    r{11}.Reactants(2).InitialAmount = infConc;              % BfC1
    r{12}.Reactants(2).InitialAmount = infConc;              % BfC1

    % Display the Completed Model Objects
    model

    % Display the Reaction Objects
    model.Reactions

    % Display the Species Objects
    model.Species

    % Simulate ODE and plot for 20 minutes
    cs = getconfigset(model,'active');
    % allowed error rate tolerance value in pico molars
    absTol = 1e-6;
    relTol = 1e-6;
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    % ode simulation type and stop time
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);

    % plot the phase of specie 
    figure(figNo); 
    hold on; box on;
    plot3(X(:,1)./1e-9, X(:,5)./1e-9,  X(:,8)./1e-9, 'LineWidth', 2.0);
    zlabel('C (nM)'); ylabel('A (nM)'); xlabel('B (nM)');
    legend('ideal CRN', 'DNA CRN');
    set(gca, 'LineWidth', 2.0); 
    
    % plot the specie concentration over time
    figure(figNo+1); hold on; box on;
    plot(t./3600, X(:,1)./1e-9, 'LineWidth', 2.0);
    plot(t./3600, X(:,5)./1e-9, 'LineWidth', 2.0);
    plot(t./3600, X(:,8)./1e-9, 'LineWidth', 2.0);
    legend('ideal A', 'ideal B', 'ideal C', 'A', 'B', 'C');
    ylabel('Concentration (nM)'); xlabel('Time (hours)');
    set(gca, 'LineWidth', 2.0); 
end