% TTK4135 - Helicopter lab
% Hints/template for problem 2.
% Updated spring 2018, Andreas L. Fl�ten

%% Initialization and model definition
init; % Change this to the init file corresponding to your helicopter

% Continous time state space model
A_c = [0 1 0 0 0 0;0 0 -K_2 0 0 0;0 0 0 1 0 0; 0 0 -K_1*K_pp -K_1*K_pd 0 0; 0 0 0 0 0 1;0 0 0 0 -K_3*K_ep -K_3*K_ed];
B_c = [0 0;0 0;0 0;K_1*K_pp 0;0 0;0 K_3*K_ep];

% Discrete time system model. x = [lambda r p p_dot]'
delta_t	= 0.25; % sampling time
A_d = delta_t*A_c+eye(6);
B_d = delta_t*B_c;

% Number of states and inputs
mx = size(A_d,2); % Number of states (number of columns in A)
mu = size(B_d,2); % Number of inputs(number of columns in B)
nx = mx;
nu = mu;

% Initial values
x1_0 = pi;                              % Lambda
x2_0 = 0;                               % r
x3_0 = 0;                               % p
x4_0 = 0;                               % p_dot
x5_0 = 0;                               % e
x6_0 = 0;                               % e_dot
x0 = [x1_0 x2_0 x3_0 x4_0 x5_0 x6_0]';  % Initial values

% Time horizon and initialization
N  = 40;                                % Time horizon for states
M  = N;                                 % Time horizon for inputs
n = N*nx+M*nu;
z  = zeros(n,1);                        % Initialize z for the whole horizon for optimization           
z0 = [x0 ; zeros(n-6,1)];

%Constraint function
alpha = 0.2;
beta =20;

% Initialization of lqr matrices
q1 = 1;
q2 = 1;
R = [q1 0;0 q2];
Q1 = [1 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0;
      0 0 0 0 0 0];
G = gen_q(Q1,R,N,M);
A_eq = gen_aeq(A_d,B_d,N,nx,nu);
B_eq = zeros(size(A_eq,1),1);
B_eq(1:nx) = A_d*x0;                      %definere x0

% Bounds
pk      = 30*pi/180;                     % boundry on pitch
ul 	    = [-pk; -inf];                   % Lower bound on control
uu 	    = [pk; inf];                    % Upper bound on control

xl(1:nx,1)      = -Inf*ones(mx,1);              % Lower bound on states (no bound)
xu(1:nx,1)      = Inf*ones(mx,1);               % Upper bound on states (no bound)
xl(3)   = ul(1);                           % Lower bound on state x3
xu(3)   = uu(1);                           % Upper bound on state x3

% Generate constraints on measurements and inputs
[vlb,vub]       = gen_constraints(N,M,xl,xu,ul,uu); % hint: gen_constraints
vlb(N*mx+M*mu)  = 0;                    % We want the last input to be zero
vub(N*mx+M*mu)  = 0;                    % We want the last input to be zero

% Objective function and constraint
f = @(z) z'*G*z;
opt = optimoptions('fmincon','Algorithm','sqp','MaxFunEvals',400000);

%% Solve QP problem with SQP solver
tic
[z,lambda] = fmincon(f,z0,[],[],A_eq,B_eq,vlb,vub,@constraint, opt); % hint: fmincon. Type 'doc fmincon' for more info 
toc

%% Extract control inputs and states
u_pc  = [z(N*nx+1:nu:n);z(n-1)];        % Pitch control input from solution
u_ec = [z(N*nx+2:nu:n);z(n)];           % Elevation control input from solution

x1 = [x0(1);z(1:mx:N*mx)];              % State x1 from solution
x2 = [x0(2);z(2:mx:N*mx)];              % State x2 from solution
x3 = [x0(3);z(3:mx:N*mx)];              % State x3 from solution
x4 = [x0(4);z(4:mx:N*mx)];              % State x4 from solution
x5 = [x0(5);z(5:mx:N*mx)];              % State x5 from solution
x6 = [x0(6);z(6:mx:N*mx)];              % State x6 from solution

num_variables = 5/delta_t;
zero_padding = zeros(num_variables,1);
unit_padding  = ones(num_variables,1);

% Optimal control inputs
u_pc   = [zero_padding; u_pc; zero_padding]; % u_pc*
u_ec   = [zero_padding; u_ec; zero_padding]; % u_ec*

% Optimal states
x1  = [pi*unit_padding; x1; zero_padding];   % Lambda*
x2  = [zero_padding; x2; zero_padding];      % r*
x3  = [zero_padding; x3; zero_padding];      % p*
x4  = [zero_padding; x4; zero_padding];      % p_dot*
x5  = [zero_padding; x5; zero_padding];      % e*
x6  = [zero_padding; x6; zero_padding];      % e_dot*

t = 0:delta_t:delta_t*(length(u_pc)-1);

xmodel = [t' x1 x2 x3 x4 x5 x6];             % Optimal states
umodel = [t' u_pc u_ec];                     % Optimal control inputs

%% Plotting optimal trajectory and control inputs
 figure(2)
 subplot(711)
 stairs(t,u_ec),grid
 ylabel('{u_{ec}}^*')
 subplot(712)
 stairs(t,u_pc),grid
 ylabel('{u_{pc}}^*')
 subplot(713)
 plot(t,x1,'m',t,x1,'m'),grid
 ylabel('{lambda}^*')
 subplot(714)
 plot(t,x2,'m',t,x2','m'),grid
 ylabel('{r}^*')
 subplot(715)
 plot(t,x3,'m',t,x3,'m'),grid
 ylabel('{p}^*')
 subplot(716)
 plot(t,x4,'m',t,x4','m'),grid
 xlabel('tid (s)'),ylabel('{pdot}^*')
 subplot(717)
 plot(t,x5,'m',t,x5','m'),grid
 xlabel('tid (s)'),ylabel('{e}^*')
 subplot(717)
 plot(t,x6,'m',t,x6','m'),grid
 xlabel('tid (s)'),ylabel('{edot}^*')