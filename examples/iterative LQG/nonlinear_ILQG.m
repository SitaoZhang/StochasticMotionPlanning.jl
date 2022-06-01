% This piece of code is adapted from: https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization
% Citation: Yuval (2022). iLQG/DDP trajectory optimization (https://www.mathworks.com/matlabcentral/fileexchange/52069-ilqg-ddp-trajectory-optimization), MATLAB Central File Exchange. Retrieved May 23, 2022.


function nonlinear_ILQG

clc;

% make stable nonlinear dynamics
dt = 0.01;          % time step
n = 2;              % state dimension
m = 1;              % control dimension

A = [-1 2; -3 -1];
A = dt*A + eye(length(A));
% A = expm(h*A);  % discrete time
B = dt*[0.5; -2];

% quadratic costs
Q = eye(n);
R = eye(m);

% control limits
Op.lims = ones(m,1)*[-1 1]*.6;

% optimization problem
DYNCST  = @(x,u,i) nonlin_dyn_cst(x,u,A,B,Q,R);
T       = 180;             % horizon
x0      = [-10; 10];        % initial state
h0      = 2*sin(2*(0:h:(T-1)*h));    
u0      = h0;               % initial controls

% run the optimization
iLQG(DYNCST, x0, u0, Op);

%%
function [f,c,fx,fu,fxx,fxu,fuu,cx,cu,cxx,cxu,cuu] = nonlin_dyn_cst(x,u,A,B,Q,R)
    u(isnan(u)) = 0;
    dis = 1.3;
    if nargout == 2
        % Boundary of obstacle
        xleft = -8;
        xright = -3;
        yup = 5.8;
        ylow = 0;
        if(x(1) > xleft-dis && x(1) < xright+dis && x(2) > ylow-dis && x(2) < yup+dis)
            if(x(1)-(xleft-dis) <= xright+dis-x(1))
                x(1) = xleft-dis;
            else
                x(1) = xright+dis;
            end
            if(x(2)-(ylow-dis) <= yup+dis-x(2))
                x(2) = ylow-dis;
            else
                x(2) = yup+dis;
            end
        end
        
        f = A*x + B*u + 0.01*[0; -0.05].*(x(2))^3 + 0.1*B.*0.01*(x(2))^2*randn(1, 1);
        c = sum(x.*(Q*x),1) + sum(u.*(R*u),1);
        
    else 
        N   = size(x,2);
        fx  = repmat(A, [1 1 N]);
        fu  = repmat(B, [1 1 N]);
        cx  = Q*x;
        cu  = R*u;
        cxx = repmat(Q, [1 1 N]);
        cxu = repmat(zeros(size(B)), [1 1 N]);
        cuu = repmat(R, [1 1 N]);
        [f,c,fxx,fxu,fuu] = deal([]);
    end


