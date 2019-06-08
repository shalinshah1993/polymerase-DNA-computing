function crnnoncatalytic(concA, concB, stopTime, figNo)
    % Adopted from Kevin et al. (2011) Molecular System Biology and Teruo et
    % al. (2013) ACS Nano
    %
    % It assumes excess nicking enzyme and uses Bst I polymerase at 37 deg C 
    % This will be faster in practice with Bst II at 55 deg C. Rate is
    % 0.0075 sec^-1
    rateK = 7.5 * 1e-4;
    % 1x = 50 nM used for the gate concentration
    baseConc = 50;

    % Create the amplification model
    model = sbiomodel('Polymerase-based strand displacement amplification CRN');

    % Add reaction set X + Y -> Z to the solver
    r1 = addreaction(model,'A + Gab -> Ia');
    r2 = addreaction(model,'B + Gba -> Ib');
    r3 = addreaction(model,'A + Ib -> O + W');
    r4 = addreaction(model,'B + Ia -> O + W');

    % Add waste reactions to the solver
    r5 = addreaction(model,'A + Gba -> W');
    r6 = addreaction(model,'B + Gab -> W');

    % Set the Kinetic Law for Reaction 1.
    kl1 = addkineticlaw(r1, 'MassAction');
    kl2 = addkineticlaw(r2, 'MassAction');
    kl3 = addkineticlaw(r3, 'MassAction');
    kl4 = addkineticlaw(r4, 'MassAction');
    kl5 = addkineticlaw(r5, 'MassAction');
    kl6 = addkineticlaw(r6, 'MassAction');

    % Add rate constant parameters
    % PMSD assume atleast 10x faster than TMSD
    p1 = addparameter(kl1, 'c1', 'Value', rateK);
    p2 = addparameter(kl2, 'c2', 'Value', rateK);
    p3 = addparameter(kl3, 'c3', 'Value', rateK);
    p4 = addparameter(kl4, 'c4', 'Value', rateK);
    p5 = addparameter(kl5, 'c5', 'Value', rateK);
    p6 = addparameter(kl6, 'c6', 'Value', rateK);

    kl1.ParameterVariableNames = {'c1'};
    kl2.ParameterVariableNames = {'c2'};
    kl3.ParameterVariableNames = {'c3'};
    kl4.ParameterVariableNames = {'c4'};
    kl5.ParameterVariableNames = {'c5'};
    kl6.ParameterVariableNames = {'c6'};

    % Set initial amounts for species in Reaction 1
    r1.Reactants(1).InitialAmount = concA*baseConc;             % A
    r2.Reactants(1).InitialAmount = concB*baseConc;            % B
    r1.Reactants(2).InitialAmount = baseConc;             % Gate ab
    r2.Reactants(2).InitialAmount = baseConc;             % Gate ba


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
    figure(figNo); hold on; box on
    plot(t, X(:,7), 'LineWidth', 2.0);
    title('Polymerase-based non-catalytic CRN');
    ylabel('Concentration (nM)'); xlabel('Time (s)')
    set(gca, 'LineWidth', 2.0);
end