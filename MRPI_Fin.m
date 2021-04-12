%% (Simsteps-N) step Robust Reachable Set with Policy u = Kx 
% Monimoy Bujarbaruah
% Akhil Shetty 

function P = MRPI_Fin(model, X, W, simsteps, N)

    X0 = X;                                 % initial set constraint
    for j = 1:simsteps-N
        % subtract noise    
        S = X0 - W;
        % backward reachable set
        R = model.reachableSet('X', S, 'direction', 'backward');
        % intersect with the state constraints
        P = R.intersect(X0);
        X0 = P;
    end

end
