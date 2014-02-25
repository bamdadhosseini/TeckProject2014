%% contraction and interpolation operator for Multi-Grid

function [I] = getMgOperator( opname, N )

    % N - total number of grid points including boundary
    % n - is the number of grid points on current grid 
    %     without the boundary points.
    
    n = N-2;

    if strcmpi('coarse2fine', opname)
        e = ones(2*n+1 ,1);
        I = speye( 2*N-1, N );
        II = 1/2*spdiags([e, 2*e, e], [-1, 0,1], 2*n+1, 2*n+1);
        I(2:end-1,2:end-1) = II(2:2:end,:)';
        %size(I(2:2:end,:))
        %full(I(2:2:end,:))
%         I = [zeros(1,n); I(2:2:end,:)'; zeros(1,n)];
%         % boundary points
%         I(1,1) =1 ;
%         I(end,end) = 1;
    elseif strcmpi('fine2coarse', opname)
        e = ones(n ,1);
        I = speye(floor(N/2)+1,N);
        II = 1/4*spdiags([e, 2*e, e], [-1, 0,1], n , n);
        %I(2:end-1,2:end-1) = [zeros(1,n); I(2:2:end,:); zeros(1,n)];
        I(2:end-1,2:end-1) = II(2:2:end,:);
        % boundary points
        %I(1,1) = 1;
        %I(end,end) = 1;
    end

end