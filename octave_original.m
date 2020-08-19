% This implements the solution of a
% multivariate rational expectations model,
% using the example from Problem Set 7, parts (g), (h) and (i).
% The computations performed in the following code therefore
% correspond to pages 41 - 51 of the book
% by Roger E. A. Farmer,
% "Macroeconomics of Self-fulfilling Prophecies", 2nd edition, MIT Press, 1999.
% The comments refer to page numbers and equation numbers from the book.
% Thomas Hintermaier, hinterma@uni-bonn.de, December 16, 2019

% Specify paramter values used in the model
% (Remark: Notation is in line with symbols in the book.)
beta  = 0.95; % discount factor
alpha = 0.3;  % capital share
delta = 0.1;  % depreciation rate
rho   = 0.9;  % autocorrelation of productivity

% Getting non-stochastic steady-state values
% book: p. 46; lecture slides: Part 9, Slide 13
k_bar = ((1/beta - 1 + delta)/alpha)^(1/(alpha - 1));
c_bar = k_bar^alpha  - delta*k_bar;
s_bar = 1;

% Describing relationships underlying linear coefficients
% book: p. 47; lecture slides: Part 9, Slide 15.
a1 = beta*alpha*(alpha - 1)*k_bar^(alpha - 1);
a2 = beta*alpha*k_bar^(alpha - 1);
b1 = 1 - delta + alpha*k_bar^(alpha - 1);
b2 = k_bar^(alpha - 1);
b3 = -c_bar/k_bar;

% Writing dynamic equilibrium conditions into matrix form

% m371 stands for the 1st matrix showing up in equation (3.7),
% i.e. the matrix multiplying period-t variables on the LHS of (3.7).
% book: p. 47; lecture slides: Part 9, Slide 15.

m371 = zeros(3,3);
m371(1,1) = -1;
m371(2,1) = b3;
m371(2,2) = b1;
m371(2,3) = b2;
m371(3,3) = rho;

% m372 stands for the 2nd matrix showing up in equation (3.7)
m372 = zeros(3,3);
m372(1,1) = -1;
m372(1,2) = a1;
m372(1,3) = a2;
m372(2,2) = 1;
m372(3,3) = 1;

% m38 stands for the inverse of m371
% book: equation (3.8), p. 47; lecture slides: Part 9, Slide 16 (top).
m38 = inv(m371);

% Matrix A
% book: equation (3.9), p. 47; lecture slides: Part 9, Slide 16.
% Note that this encodes all the relevant dynamics for the solution.
A = m38*m372;

% Calculate Eigenvectors and Eigenvalues of A
[Q_temp,Lam_temp] = eig(A); % These are labeled "temp",
                            % for being temporary in the sense that we still
                            % need to make sure that the ORDERING of eigenvalues
                            % puts the forward-stable one in the first position

% The next lines are there to identify the appropriate eigenvalue,
% which is forward stable, and therefore can be used to derive
% a restriction to express a free variable as an equilibrium function of predetermined variables
diagonal_Lam_temp = diag(Lam_temp);
select_i_stable = abs(diagonal_Lam_temp) < 1;
mark_unstable = not(select_i_stable);

% Reordering to make sure the stable eigenvalue and the corresponding eigenvector
% are in the first positions in the matrix of eigenvalues and the matrix of eigenvectors.
Lam = diag([diagonal_Lam_temp(select_i_stable);diagonal_Lam_temp(mark_unstable)]);
Q = zeros(3,3);
Q(:,1) =     Q_temp(:,select_i_stable);
Q(:,2:end) = Q_temp(:,mark_unstable);
% The eigenvalue less than one in absolute value is in the first position now
% book: p. 50, equation (3.12); lecture slides: Part 10, Slide 7.

% The inverse of Q, as used in the decomposition A = Q*Lam*inv(Q)
% book, p. 50 (top), lecture slides: Part 10, Slide 4.
Q_inv = inv(Q);

% The key restriction: getting a free variable (consumption)
% as an equilibrium function of predetermined variables (capital, productivity-state).
% book: p. 51, lecture slides: Part 10, Slide 9.
disp('coefficients of rational expectations equilibrium function');
c_coeff_k = -Q_inv(1,2)/Q_inv(1,1)
c_coeff_s = -Q_inv(1,3)/Q_inv(1,1)

% % disp('Q');
% % Q
% %
% % disp('Lam');
% % Lam
% %
% % disp('inverse of Q');
% % Q_inv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part (g). Simulate solution for 50 periods.
T_g = 50;

% Create vector for productivity values.
s_t = NaN*zeros(1,T_g);

% First 10 periods given in part (f)
s_t(1,1:10) = [0, -0.005, -0.009, -0.013, -0.022, -0.021, -0.019, -0.011, -0.012, -0.003];

% Obtain innovations for the productivity process for subsequent periods.
var_v = 0.007^2;
v_t = NaN*zeros(T_g,1);
v_t(11:end)=sqrt(var_v)*randn(length(v_t)-10,1);


% Obtain productivity levels
for t = 11:T_g

    s_t(t) = rho*s_t(t-1) + v_t(t);

end

% Generating the time series of variables, which are equilibria for the given shock realizations
c_t = NaN*zeros(T_g,1);
k_t = NaN*zeros(T_g+1,1);
k_t(1) = 0;

for t = 1:T_g

    % Imposing the key equilibrium restriction from rational expectations:
    % The free variable (consumption) is expressed as a function of
    % the predetermined variables (capital and productivity state)
    % book: p. 51, lecture slides: Part 10, Slide 12.
    % This is the recursive structure described in Part 10, Slide 12.
%
    c_t(t) = -Q_inv(1,2)/Q_inv(1,1)*k_t(t) -Q_inv(1,3)/Q_inv(1,1)*s_t(t);
%     c_t(t) = c_coeff_k*k_t(t) + c_coeff_s*s_t(t);
    k_t(t+1) = b1*k_t(t) + b2*s_t(t) + b3*c_t(t);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part (h). Calculate output deviations for these 50 periods.

% Note that log linearization of production function y_t = s_t*(k^{alpha}_t)
% yields \hat{y}_t = \hat{s}_t + alpha*\hat{k}_t

y_t = s_t' + alpha.*k_t(1:(end-1),1);

% Plot results.
figure(10);

t = 0:1:(T_g-1);
t2 = 0:1:T_g;

plot(t,s_t,'r-x');
hold on;
plot(t,c_t,'b-x');
hold on;
plot(t2,k_t,'g-x');
hold on;
plot(t,y_t,'k-x');

title('productivity state, consumption, capital,output','FontSize',24);
legend('productivity state','consumption','capital','output');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Part (i). Simulate for 1000 periods and compute correlations.

T_i = 1000;

% Create vector for productivity values.
s_t_i = NaN*zeros(1,T_i);

% First 10 periods given in part (f)
s_t_i(1,1:10) = [0, -0.005, -0.009, -0.013, -0.022, -0.021, -0.019, -0.011, -0.012, -0.003];

% Obtain innovations for the productivity process for subsequent periods.
var_v = 0.007^2;
v_t_i = NaN*zeros(T_i,1);
v_t_i(11:end)=sqrt(var_v)*randn(length(v_t_i)-10,1);


% Obtain productivity levels
for t = 11:T_i

    s_t_i(t) = rho*s_t_i(t-1) + v_t_i(t);

end

% Generating the time series of variables, which are equilibria for the given shock realizations
c_t_i = NaN*zeros(T_i,1);
k_t_i = NaN*zeros(T_i+1,1);
k_t_i(1) = 0;

for t = 1:T_i

    % Imposing the key equilibrium restriction from rational expectations:
    % The free variable (consumption) is expressed as a function of
    % the predetermined variables (capital and productivity state)
    % book: p. 51, lecture slides: Part 10, Slide 12.
    % This is the recursive structure described in Part 10, Slide 12.

    c_t_i(t) = -Q_inv(1,2)/Q_inv(1,1)*k_t_i(t) -Q_inv(1,3)/Q_inv(1,1)*s_t_i(t);
%     c_t_i(t) = c_coeff_k*k_t_i(t) + c_coeff_s*s_t_i(t);
    k_t_i(t+1) = b1*k_t_i(t) + b2*s_t_i(t) + b3*c_t_i(t);

end

% Obtain output deviations.
y_t_i = s_t_i' + alpha.*k_t_i(1:(end-1),1);

% Obtain output-capital correlation.

correlation_y_k = corr(y_t_i,k_t_i(1:(end-1),1));
disp('correlation between the deviation of output from its steady-state value and the deviation of capital from its steady-state value');
correlation_y_k

% Plot results.
figure(11);

t = 0:1:(T_i-1);
t2 = 0:1:T_i;

plot(t,s_t_i,'r-x');
hold on;
plot(t,c_t_i,'b-x');
hold on;
plot(t2,k_t_i,'g-x');
hold on;
plot(t,y_t_i,'k-x');

title('productivity state, consumption, capital,output','FontSize',24);
legend('productivity state','consumption','capital','output');
hold off;








