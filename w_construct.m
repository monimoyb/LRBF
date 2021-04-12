%% Uncertainty Support Construction in the beginning of an iteration
% Monimoy Bujarbaruah
% Akhil Shetty
%%
% Sets a high confidence value. Uses data. Constructs support. Checks
% feasibility of MPC problem from x_0. If infeasible, lowers confidence and
% repeats. if feasible, this support is rolled out to the main code to
% complete an iteration. Then data is collected and process is repeated

function [w_lb, w_ub, Xn, Pinf, conf_possible, scal_val] = w_construct(w_samples, conf, nx,nu, A, B, C, D, b, Q, R, U, N, x_0, X, simsteps, options)

    w_maxSeen =  max(abs(w_samples),[],2);
    feas_conf = 0;
    count = 0; 

    while feas_conf == 0     
        conf = conf/(2^count);
        alpha = 1/(1-conf)^(1/size(w_samples,2));                                       % for uniform distribution

        l_maxScaled = alpha*w_maxSeen;
        w_lb = -l_maxScaled; w_ub = l_maxScaled;                                    
        W = Polyhedron('lb',w_lb,'ub',w_ub); 

        %%% GETTING THE TERMINAL SET
        [Finf, Pinf] = dlqr(A,B,Q,R);  
        model = LTISystem('A',A-B*Finf,'B',B); 
        Hx = X.A; Hu = U.A; 
        hx = X.b; hu = U.b; 
        Hxaut = [Hx; -Hu*Finf]; hxaut = [hx; hu]; 

        S = Polyhedron('A',Hxaut,'b',hxaut); 
        Xn = MRPI_Fin(model, S, W, simsteps,N);                           

        if isEmptySet(Xn) == 0

            %% Forming matrices appearing in the optimization problem 
            dim_t = size(C,1)*N + size(Xn.A,1);  
            [capA, capE, capB, capC, capD, ~, ~, ~] = obtain_matR(A,B,C,D,Xn,nx,nu,N,dim_t);
            matF = capC*capB + capD; 
            matG = capC*capE; 
            matH = -capC*capA; 
            mat_c = [kron(ones(N,1),b); Xn.b];

            constraints = []; 

           %%  Open Loop W Variables
           Hs=[]; hs =[]; 
            for k = 1:N
                polS = W; 
                Hs_ol = polS.A;  hs_ol = polS.b; 
                Hs = blkdiag(Hs, Hs_ol); hs = [hs; hs_ol];
            end
           dim_a = size(Hs,1);                                 

           %% Creating Open Loop Optimization Variables for MPC 
            M = sdpvar(nu*N,nx*N); 
            for j=1:nu*N
                for k = 2*j-1:nx*N
                        M(j,k) =0;
                end
            end
            v = sdpvar(nu*N,1); 
            Z = sdpvar(dim_a,dim_t);     
            
            %% Checking Feasibility 
            constraints = [constraints; matF*v + Z'*hs <= mat_c + matH*x_0];
            constraints = [constraints; Z>=0];
            constraints = [constraints, matF*M + matG == Z'*Hs];
            diagn=solvesdp(constraints, [], options);
            if  diagn.problem == 0
                feas_conf = 1;
            else
                feas_conf = 0;
            end
        else              
            feas_conf = 0;              
        end

        count = count + 1;

    end

    conf_possible = conf; 
    scal_val = alpha; 

end