% crnidealamp(10, 16000, 1);
crnautocatalyticamp(10, 7200, 2);

function crnidealamp(concA, stopTime, figNo)
    % enter A and B concentration (in nM) assuming all gates are several
    % micro molars
    %
    % stopTime is the total simulation time
    %
    % figNo is the display figure number
    
    % rate of reactions 
    rate = 1e-3;
    % allowed error rate tolerance value in pico molars
    absTol = 1e-8;
    relTol = 1e-8;
    
    % Create the amplification model
    model = sbiomodel('Polymerase-based strand displacement CRN');

    r = cell(1, 1);
    k = cell(1, length(r));
    p = cell(1, length(r));
         
    % Add reaction set X -> 2X to the solver
    r{1} = addreaction(model,'A -> A + A');
        
    % Set the Kinetic Law for Reactions.
    for i = 1:length(r)
        k{i} = addkineticlaw(r{i}, 'MassAction'); 
    end
    
    p{1} = addparameter(k{1}, 'c1', 'Value', rate);
    
    % Set the Kinetic Law for Reactions.
    for i = 1:length(r)
        k{i}.ParameterVariableNames = {strcat('c',num2str(i))};
    end
    
    % Set initial amounts for species 
    r{1}.Reactants(1).InitialAmount = concA*1e-8;        % A

    % Display the Completed Model Objects
    model

    % Display the Reaction Objects
    model.Reactions

    % Display the Species Objects
    model.Species

    % Simulate ODE and plot for totalTime
    cs = getconfigset(model,'active');
    set(cs.SolverOptions, 'AbsoluteTolerance', absTol);
    set(cs.SolverOptions, 'RelativeTolerance', relTol);
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);

    
    % plot the specie concentration over time
    figure(figNo); 
    box on; hold on;
    plot(t./3600, X(:,1), ':r', 'LineWidth', 2.0);
    ylabel('Concentration'); xlabel('Time');
    set(gca, 'LineWidth', 2.0); 
end

function crnautocatalyticamp(concA, stopTime, figNo)
    
    k = 1.0;
    alpha = 1e-3; beta = 1e-8;
    % slow down all reactions by multiplying with a scaling factor 1e-3
    rateConst = k * (alpha);
    % fastest reaction are assumed to occur 
    infRate = 1e5 * (k);
    % inf concentration of gates will be 100 uM
    infConc = 1e3 * (beta);
    % allowed error rate tolerance value in pico molars
    absTol = 1e-12;
    relTol = 1e-12;
    
    % Create the amplification model
    model = sbiomodel('Polymerase-based strand displacement amplification CRN');

    % Add reaction set X -> kX to the solver
    % CHANGE k here to make it different
    r1 = addreaction(model,'A + Ga0 -> Ia');
    r2 = addreaction(model,'Ia + Ga1 -> A + A + A + A + A + A ');
    
    % Set the Kinetic Law for Reaction 1.
    kl1 = addkineticlaw(r1, 'MassAction');
    kl2 = addkineticlaw(r2, 'MassAction');
    
    % Add rate constant parameters
    % PMSD assume atleast 10x faster than TMSD
    p1 = addparameter(kl1, 'c1', 'Value', rateConst/infConc);
    p2 = addparameter(kl2, 'c2', 'Value', infRate);
    
    kl1.ParameterVariableNames = {'c1'};
    kl2.ParameterVariableNames = {'c2'};
    
    % Set initial amounts for species in Reaction 1
    r1.Reactants(1).InitialAmount = concA*beta;               % A
    r1.Reactants(2).InitialAmount = infConc;             % Gate ab
    r2.Reactants(2).InitialAmount = infConc;             % Gate ba

    % Display the Completed Model Objects
    model

    % Display the Reaction Objects
    model.Reactions

    % Display the Species Objects
    model.Species

    % Simulate ODE and plot for 20 minutes
    cs = getconfigset(model,'active');
    cs.SolverType = 'ode15s';
    cs.StopTime = stopTime;
    [t,X] = sbiosimulate(model);

    % plot X over the time
    figure(figNo); hold on; box on;
    plot(t./3600, X(:, 1), 'LineWidth', 2.0);
    title('Polymerase-based autocatalytic CRN');
    ylabel('Concentration (nM)'); xlabel('Time (s)')
    set(gca, 'LineWidth', 2.0)
end