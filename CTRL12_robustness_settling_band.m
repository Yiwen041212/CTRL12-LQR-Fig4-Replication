% CTRL12: Robustness test for the LQR tuning result
% Sensitivity of the selected Q and K to the position settling-time band.
% This version does NOT require Control System Toolbox.
% The LQR gain is obtained from the continuous-time algebraic Riccati
% equation using the Hamiltonian matrix, and the closed-loop response is
% simulated with the base-MATLAB solver ode45.
%
% Baseline reproduction settings:
%   Q0 = diag([800, 200, 10, 100]), R = 0.1
%   Five cumulative Q increments: [0, 1000, 3000, 6000, 10000]
%   xref = [0.05, 0, 0, 0]'
%
% Robustness settings:
%   Position settling bands = 1%, 2%, and 5% of the 0.05 m reference.
%   The angle reference is zero, so a percentage band is undefined.
%   Therefore, the angle band is fixed at 5e-4 rad for all tests.

clear; clc; close all;

%% 1. Linearised inverted-pendulum model
M = 2.4; % Cart mass [kg]
m = 0.23;   % Pendulum mass [kg]
l = 0.3; % Pendulum half-length [m]
g = 9.81; % Gravitational acceleration [m/s^2]
I = 0.0017;  % Moment of inertia [kg.m^2]
b = 0.001;  % Friction coefficient [N.s/m]

q = (I + m*l^2)*(M + m) - (m*l)^2;

A = [0, 1, 0, 0;
     0, -(I + m*l^2)*b/q, -m^2*l^2*g/q, 0;
     0, 0, 0, 1;
     0, m*l*b/q, (M + m)*m*l*g/q, 0];

B = [0;
     (I + m*l^2)/q;
     0;
     -m*l/q];

C = [1, 0, 0, 0;
     0, 0, 1, 0];
D = zeros(2, 1);

openLoopPoles = eig(A);
assert(any(real(openLoopPoles) > 0), ...
    'The open-loop model should contain an unstable pole.');

fprintf('Open-loop poles:\n');
disp(openLoopPoles);

%% 2. Reproduce the five LQR candidates once
Q0 = diag([800, 200, 10, 100]);
R = 0.1;
xref = [0.05; 0; 0; 0];
x0 = zeros(4, 1);
t = (0:0.005:5)';

qIncrements = [0, 1000, 3000, 6000, 10000];
nTrials = numel(qIncrements);

Qall = zeros(4, 4, nTrials);
Kall = zeros(nTrials, 4);
xall = zeros(numel(t), nTrials);
thetaall = zeros(numel(t), nTrials);
closedLoopMaxRealPole = zeros(nTrials, 1);

for trial = 1:nTrials
    Q = Q0;
    Q(1,1) = Q0(1,1) + qIncrements(trial);
    Q(3,3) = Q0(3,3) + qIncrements(trial);

    K = continuous_lqr_no_toolbox(A, B, Q, R);
    Acl = A - B*K;
    referenceInput = B*(K*xref);

    % Simulate xdot = (A-BK)x + BK*xref. This is equivalent to the
    % tracking law u = -K(x-xref). ode45 is part of base MATLAB.
    stateEquation = @(~, state) Acl*state + referenceInput;
    [~, stateResponse] = ode45(stateEquation, t, x0);
    y = stateResponse*C';

    Qall(:,:,trial) = Q;
    Kall(trial,:) = K;
    xall(:,trial) = y(:,1);
    thetaall(:,trial) = y(:,2);
    closedLoopMaxRealPole(trial) = max(real(eig(Acl)));
end

assert(all(closedLoopMaxRealPole < 0), ...
    'At least one LQR candidate is not closed-loop stable.');

%% 3. Robustness test: vary the position settling band
bandPercent = [1, 2, 5];
xTolerance = (bandPercent/100)*abs(xref(1));
thetaTolerance = 5e-4;  % Absolute band because theta_ref = 0 rad
nBands = numel(bandPercent);

% Columns:
% BandPercent, xTolerance, thetaTolerance, Trial, Q11, Q33,
% xTs, thetaTs, xMax, xMin, thetaMax, thetaMin, S, IsSelected
robustnessData = zeros(nBands*nTrials, 14);
SbyBand = zeros(nBands, nTrials);
selectedTrial = zeros(nBands, 1);

row = 0;
for bandIdx = 1:nBands
    for trial = 1:nTrials
        row = row + 1;

        x = xall(:,trial);
        theta = thetaall(:,trial);

        xTs = settling_time_from_band(t, x, xref(1), ...
            xTolerance(bandIdx));
        thetaTs = settling_time_from_band(t, theta, 0, ...
            thetaTolerance);

        if isnan(xTs) || isnan(thetaTs)
            S = Inf;
        else
            S = xTs^2 + max(x)^2 + min(x)^2 + ...
                thetaTs^2 + max(theta)^2 + min(theta)^2;
        end

        SbyBand(bandIdx,trial) = S;
        robustnessData(row,:) = [bandPercent(bandIdx), ...
            xTolerance(bandIdx), thetaTolerance, trial, ...
            Qall(1,1,trial), Qall(3,3,trial), ...
            xTs, thetaTs, max(x), min(x), ...
            max(theta), min(theta), S, 0];
    end

    [~, selectedTrial(bandIdx)] = min(SbyBand(bandIdx,:));
    selectedRow = (bandIdx - 1)*nTrials + selectedTrial(bandIdx);
    robustnessData(selectedRow,14) = 1;
end

variableNames = {'BandPercent','xTolerance','thetaTolerance', ...
    'Trial','Q11','Q33','xTs','thetaTs','xMax','xMin', ...
    'thetaMax','thetaMin','S','IsSelected'};
robustnessTable = array2table(robustnessData, ...
    'VariableNames', variableNames);

disp('Robustness results for all tolerance settings:');
disp(robustnessTable);

%% 4. Compact summary of the selected controller under each band
selectedQ11 = zeros(nBands,1);
selectedQ33 = zeros(nBands,1);
selectedS = zeros(nBands,1);

for bandIdx = 1:nBands
    trial = selectedTrial(bandIdx);
    selectedQ11(bandIdx) = Qall(1,1,trial);
    selectedQ33(bandIdx) = Qall(3,3,trial);
    selectedS(bandIdx) = SbyBand(bandIdx,trial);
end

summaryTable = table(bandPercent', xTolerance', ...
    thetaTolerance*ones(nBands,1), selectedTrial, ...
    selectedQ11, selectedQ33, selectedS, ...
    'VariableNames', {'BandPercent','xTolerance','thetaTolerance', ...
    'SelectedTrial','SelectedQ11','SelectedQ33','MinimumS'});

disp('Selected controller under each position settling band:');
disp(summaryTable);

fprintf('\nInterpretation:\n');
for bandIdx = 1:nBands
    fprintf(['%g%% position band: Trial %d selected, ', ...
        'Q11 = %.0f, Q33 = %.0f, minimum S = %.6f\n'], ...
        bandPercent(bandIdx), selectedTrial(bandIdx), ...
        selectedQ11(bandIdx), selectedQ33(bandIdx), ...
        selectedS(bandIdx));
end

%% 5. Export tables for the Week 8 report
writetable(robustnessTable, 'CTRL12_robustness_full_results.csv');
writetable(summaryTable, 'CTRL12_robustness_summary.csv');

%% 6. Plot penalty S against trial for all settling bands
figure('Color','w','Name','Settling-band robustness test');
hold on;
colors = lines(nBands);
markers = {'o','s','^'};

for bandIdx = 1:nBands
    plot(1:nTrials, SbyBand(bandIdx,:), ...
        'Color', colors(bandIdx,:), ...
        'Marker', markers{bandIdx}, ...
        'LineWidth', 1.6, ...
        'MarkerSize', 7, ...
        'DisplayName', sprintf('%g%% position band', ...
        bandPercent(bandIdx)));

    trial = selectedTrial(bandIdx);
    plot(trial, SbyBand(bandIdx,trial), 'p', ...
        'Color', colors(bandIdx,:), ...
        'MarkerFaceColor', colors(bandIdx,:), ...
        'MarkerSize', 12, ...
        'HandleVisibility', 'off');
end

hold off;
grid on;
xticks(1:nTrials);
xlim([0.75, nTrials + 0.25]);
xlabel('LQR candidate trial');
ylabel('Penalty S');
title('Sensitivity of LQR Selection to the Position Settling Band');
legend('Location','best');

print(gcf, 'CTRL12_robustness_S_vs_trial.png', '-dpng', '-r300');

%% 7. Save closed-loop stability evidence
stabilityTable = table((1:nTrials)', closedLoopMaxRealPole, ...
    'VariableNames', {'Trial','MaximumRealPartOfClosedLoopPoles'});
disp('Closed-loop stability check:');
disp(stabilityTable);
writetable(stabilityTable, 'CTRL12_closed_loop_stability.csv');

%% Local function
function settlingTime = settling_time_from_band(t, response, target, tolerance)
% Return the first time after which the response remains inside the band.
% NaN means that the response has not settled within the simulation window.

    lastOutside = find(abs(response - target) > tolerance, 1, 'last');

    if isempty(lastOutside)
        settlingTime = 0;
    elseif lastOutside == numel(t)
        settlingTime = NaN;
    else
        settlingTime = t(lastOutside + 1);
    end
end

function K = continuous_lqr_no_toolbox(A, B, Q, R)
%CONTINUOUS_LQR_NO_TOOLBOX Compute continuous-time LQR gain without lqr().
%
% The stabilising solution P of
%   A'P + PA - PBR^(-1)B'P + Q = 0
% is recovered from the stable invariant subspace of the Hamiltonian
% matrix. The feedback law is u = -Kx, where K = R^(-1)B'P.

    n = size(A,1);
    hamiltonian = [A, -B*(R\B');
                   -Q, -A'];

    [eigenvectors, eigenvalueMatrix] = eig(hamiltonian);
    eigenvalues = diag(eigenvalueMatrix);
    stableMask = real(eigenvalues) < 0;

    if nnz(stableMask) ~= n
        error(['Unable to isolate the stabilising Riccati solution: ', ...
            'expected %d stable Hamiltonian eigenvalues but found %d.'], ...
            n, nnz(stableMask));
    end

    stableVectors = eigenvectors(:,stableMask);
    upperBlock = stableVectors(1:n,:);
    lowerBlock = stableVectors(n+1:end,:);

    if rcond(upperBlock) < 1e-12
        error('The stable Hamiltonian invariant subspace is ill-conditioned.');
    end

    P = real(lowerBlock/upperBlock);
    P = (P + P')/2;
    K = R\(B'*P);

    % Numerical check of the algebraic Riccati equation.
    residual = A'*P + P*A - P*B*(R\B')*P + Q;
    relativeResidual = norm(residual,'fro') / max(1,norm(Q,'fro'));
    if relativeResidual > 1e-7
        warning('LQR:RiccatiResidual', ...
            'Relative Riccati residual is %.3e.', relativeResidual);
    end
end


