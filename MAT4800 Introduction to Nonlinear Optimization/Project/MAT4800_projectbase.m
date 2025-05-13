% -------------------------
%Parameters
% -------------------------
p_f = 883;             % Free trade price

b_d = 8;               % Demand intercept
m_d = 0.001716;        % Demand slope
a_d = 1.0;             % Demand curvature

b_s = 1;               % Supply intercept
m_s = 0.000189;        % Supply slope
a_s = 1.0;             % Supply curvature

t_s = 0;                  % Tariff floor
p_s = -0.2*p_f + p_f ;    % Price floor
t_c = 1;                  % Tariff cap
p_c = 0.2*p_f + p_f;      % Price cap

% -------------------------
% Shared Functions
% -------------------------
D = @(p) b_d - m_d * p.^a_d;           % Demand function
S = @(p) b_s + m_s * p.^a_s;           % Supply function
pe = @(t, p_i) p_i ./ (t + 1);         % Export price function

% -------------------------
% Optimization Setup
% -------------------------
x0 = [0.1, 850];     % initial guess
lb = [t_s, p_s];     % lower bounds
ub = [t_c, p_c];     % upper bounds

options = optimoptions('fmincon', ...
    'Display','off', ...
    'Algorithm','interior-point', ...
    'StepTolerance', 1e-10);

% -------------------------
% Without Exporter Loss
% -------------------------

fprintf('--- WITHOUT Exporter Loss ---\n');

W1 = @(x) ...                                                                % objective function without exporter loss
    (p_f - pe(x(1), x(2))) .* (D(x(2)) - S(x(2))) ...                        % G
  - 0.5 * (p_f - pe(x(1), x(2))) .* (S(x(2)) - S(pe(x(1), x(2)))) ...        % B
  - 0.5 * (p_f - pe(x(1), x(2))) .* (D(pe(x(1), x(2))) - D(x(2)));           % D

obj1 = @(x) -W1(x);

[x1_opt, fval1] = fmincon(obj1, x0, [], [], [], [], lb, ub, [], options);

t1 = x1_opt(1);
p1 = x1_opt(2);
W1_val = W1(x1_opt);
imports1 = D(p1) - S(p1);

fprintf('Optimal tariff t: %.4f\n', t1);
fprintf('Optimal price p_i: %.4f\n', p1);
fprintf('Actual welfare: %.4f million USD\n', W1_val);
fprintf('Import quantity (D - S): %.4f million tons\n\n', imports1);

% -------------------------
% With Exporter Loss
% -------------------------

fprintf('--- WITH Exporter Loss ---\n');

W2 = @(x) ...                                                               % objective function with exporter loss
    (p_f - pe(x(1), x(2))) .* (D(x(2)) - S(x(2))) ...                       % G
  - 0.5 * (p_f - pe(x(1), x(2))) .* (S(x(2)) - S(pe(x(1), x(2)))) ...       % B
  - 0.5 * (p_f - pe(x(1), x(2))) .* (D(pe(x(1), x(2))) - D(x(2))) ...       % D
  - 0.5 * (p_f - pe(x(1), x(2))) .* abs( ...                                % Exporter loss
        (S(p_f) - D(p_f)) + ...
        (S(pe(x(1), x(2))) - D(pe(x(1), x(2)))) );

obj2 = @(x) -W2(x);

[x2_opt, fval2] = fmincon(obj2, x0, [], [], [], [], lb, ub, [], options);

t2 = x2_opt(1);
p2 = x2_opt(2);
pe2 = pe(t2, p2);

W2_val = W2(x2_opt);
imports2 = D(p2) - S(p2);

fprintf('Optimal tariff t: %.4f\n', t2);
fprintf('Optimal price p_i: %.4f\n', p2);
fprintf('Actual welfare (including exporter loss): %.4f million USD\n', W2_val);
fprintf('Import quantity (D - S): %.4f million tons\n', imports2);

