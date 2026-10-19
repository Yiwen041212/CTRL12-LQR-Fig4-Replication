% CTRL_12: Inverted pendulum LQR replication (Fig. 4)
clear; clc; close all;

%% Linear model
M = 2.4;
m = 0.23;
l = 0.3;
g = 9.81;
I = 0.0017;
b = 0.001;
q = (I + m*l^2)*(M + m) - (m*l)^2;

A = [0, 1, 0, 0;
     0, -(I + m*l^2)*b/q, -m^2*l^2*g/q, 0;
     0, 0, 0, 1;
     0, m*l*b/q, (M + m)*m*l*g/q, 0];
B = [0; (I + m*l^2)/q; 0; -m*l/q];
C = [1, 0, 0, 0;
     0, 0, 1, 0];
D = zeros(2, 1);
plant = ss(A, B, C, D);

poles = eig(A);
assert(any(real(poles) > 0), 'The open-loop model should be unstable.');
fprintf('Open-loop poles:\n');
disp(poles);

%% LQR settings
Q = diag([800, 200, 10, 100]);
R = 0.1;
xref = [0.05; 0; 0; 0];
x0 = zeros(4, 1);
t = (0:0.005:5)';
nTrials = 5;
qStep = 1000;
xTol = 0.02*xref(1);
thetaTol = 5e-4;  % Fixed angle band; not specified in the paper

Kall = zeros(nTrials, 4);
Qall = zeros(4, 4, nTrials);
xall = zeros(numel(t), nTrials);
thetaall = zeros(numel(t), nTrials);
results = zeros(nTrials, 14);

%% Five tuning trials
for j = 1:nTrials
    i = j - 1;
    Q(1,1) = Q(1,1) + i*qStep;
    Q(3,3) = Q(3,3) + i*qStep;
    K = lqr(A, B, Q, R);

    Acl = A - B*K;
    Bcl = B*(K*xref);
    closedLoop = ss(Acl, Bcl, C, D);
    y = lsim(closedLoop, ones(numel(t), 1), t, x0);

    x = y(:,1);
    theta = y(:,2);

    lastX = find(abs(x - xref(1)) > xTol, 1, 'last');
    lastTheta = find(abs(theta) > thetaTol, 1, 'last');

    xTs = 0;
    thetaTs = 0;
    if ~isempty(lastX)
        if lastX == numel(t)
            xTs = NaN;
        else
            xTs = t(lastX + 1);
        end
    end
    if ~isempty(lastTheta)
        if lastTheta == numel(t)
            thetaTs = NaN;
        else
            thetaTs = t(lastTheta + 1);
        end
    end

    if isnan(xTs) || isnan(thetaTs)
        S = Inf;
    else
        S = xTs^2 + max(x)^2 + min(x)^2 + ...
            thetaTs^2 + max(theta)^2 + min(theta)^2;
    end

    Qall(:,:,j) = Q;
    Kall(j,:) = K;
    xall(:,j) = x;
    thetaall(:,j) = theta;
    results(j,:) = [j, Q(1,1), Q(3,3), K, xTs, thetaTs, ...
                    max(x), min(x), max(theta), min(theta), S];
end

names = {'Trial','Q11','Q33','K1','K2','K3','K4', ...
         'xTs','thetaTs','xMax','xMin','thetaMax','thetaMin','S'};
summary = array2table(results, 'VariableNames', names);
disp(summary);

%% Paper comparison
paperK0 = [-89.4427, -88.5363, -345.9809, -70.0729];
paperK4 = [-328.6335, -216.3436, -694.5218, -116.4222];
errK0 = max(abs(Kall(1,:) - paperK0));
errK4 = max(abs(Kall(5,:) - paperK4));
[bestS, bestTrial] = min(summary.S);
fprintf('Initial gain error: %.6g\n', errK0);
fprintf('Final gain error:   %.6g\n', errK4);
fprintf('Minimum S: trial %d, S = %.6f\n', bestTrial, bestS);
disp('Selected Q:');
disp(Qall(:,:,bestTrial));
if bestTrial ~= 5
    warning('The selected trial differs from the paper. Check settling-time conventions.');
end

%% Plot responses
colors = lines(nTrials);
figure('Color','w','Name','Fig. 4 replication');
hold on;
for j = 1:nTrials
    plot(t, xall(:,j), '-', 'Color', colors(j,:), 'LineWidth', 1.4, ...
         'DisplayName', sprintf('Trial %d: x',j));
    plot(t, thetaall(:,j), '--', 'Color', colors(j,:), 'LineWidth', 1.2, ...
         'DisplayName', sprintf('Trial %d: theta',j));
end
plot(t, xref(1)*ones(size(t)), ':k', 'HandleVisibility', 'off');
hold off; grid on;
xlim([0 5]);
xlabel('Time (s)');
ylabel('Cart position (m) / angle (rad)');
title('LQR Tuning - Fig. 4 Replication');
legend('Location', 'eastoutside');

figure('Color','w','Name','Response comparison');
subplot(2,1,1);
plot(t, xall, 'LineWidth', 1.2);
hold on;
plot(t, xref(1)*ones(size(t)), ':k');
hold off; grid on;
ylabel('Cart position (m)');
legend('Trial 1','Trial 2','Trial 3','Trial 4','Trial 5','Reference');
subplot(2,1,2);
plot(t, thetaall, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Pendulum angle (rad)');

