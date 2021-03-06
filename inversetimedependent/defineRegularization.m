%% define regularization
%
% defines the regularization operator
%
%

function R = defineRegularization( nwind, dt, bc )

% scale regularization properly
s2day = 1/(60*60*24);

e = ones(nwind,1);
sdt = dt*s2day; % convert time step size

%% Laplacian Dirichlet
%R = speye(nwind) +  (1/(2*dt)).*spdiags([e, e], [-1, 1], nwind, nwind);
R = (1/(sdt^2)).*spdiags([e, -2*e, e], [-1,0, 1], nwind, nwind);

if strcmpi(bc, 'neumann')
    % neumann
    R(1,1) = -(1/sdt^2);
    R(end,end) = -(1/sdt^2);
elseif strcmpi(bc, 'periodic')
    % periodic
    R(1,end) = (1/sdt^2);
    R(end,1) = (1/sdt^2);
end    
%R = (1/(2*dt)).*spdiags([e, -e], [-1, 1], nwind, nwind);

R = (1e-2*speye(nwind) + (R'*R))*sdt;

%R(1,1) = -(1/dt);
%R(end,end) = -(1/dt);
%R(1,end) = (1/dt^2);
%R(end,1) = (1/dt^2);


%R = (1/(dt)).*spdiags([-1*e, e], [ -1, 0], nwind+1, nwind);
%R(1,end) = -1;
%R(end,1) = 1;

end