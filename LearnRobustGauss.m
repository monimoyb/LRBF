%% Robustness Learning with MPC: Clipped Gaussian Case 
% Monimoy Bujarbaruah
% Akhil Shetty
%%
clear all
close all
clc
yalmip 'clear'
rng(3)                                          

%% Parameters to decide the support estimation method
%%% Make sure the problem is well-posed by checking wellPosedness_Check.m
% Set cvxFlag = 0 for LRBF
% Set cvxFlag = 1 and blow = 0 for cvxhull of seen disturbance samples as support estimate
% Set cvxFlag = 1 and blow = 1 for inflated cvxhull of seen disturbance samples as support estimate
% blow_size decides the inflation factor. 
cvxFlag = 0;                                 % if 1 then ONLY convex hull of seen disturbances is considered 
blow = 0;                                    % one can choose to inflate the cvx hull 
blow_size = 2; 
soft_flg = 1;                                % 1 for slacks and 0 for hard constraints 
iC  = 30;                                    % Number of iterations to run 

%% Loading all system parameters
init_SampleSize = 5;
[A,B,C,D,b,X,U,nx,nu,~,~, x_0,Q,R,N,trueMu,trueStd, x_ref,simsteps] = sys_load(); 
Wslack = 10000;                                                   % constraint softening 
mont_count = 100;                                                 % MC run count 
cost_iter = zeros(mont_count,iC); 
prob_fail = zeros(iC,mont_count);                                 % prob mass missing in \hat{W}^j in an iteration

%% Needed only for LRBF
if cvxFlag == 0
    conf = 0.9875;                                                % desired confidence value 
    conf_array = zeros(mont_count,iC); 
end

%% Monte Carlo Loop Starts Here

for mc = 1: mont_count 
    mc                                                                    
    w_samples = []; 
    gaussPol = Polyhedron('lb',trueMu-3*trueStd,'ub',trueMu+3*trueStd);

    %% CODING REJECTION SAMPLING ALGORITHM  
    while size(w_samples,2) < init_SampleSize
        smple = trueStd.* randn(nx, 1) + trueMu;
        if gaussPol.contains(smple) == 1
            w_samples = [w_samples, smple];
        else
            w_samples = w_samples;
        end
    end

    %% Iterations start here. Run these iterations MC number of times. 

    for iter_count = 1:iC           
        %%% Getting Started with Closed Loop MPC Problem 
         x_cl = x_0; 
         cost_iter(mc,iter_count) = 0; 
         options = sdpsettings('solver','gurobi','verbose',0);

        %% Dealing with only cvx hull of seen samples 
        if cvxFlag == 1
            [W, Xn, Pinf] = w_constructCVX(w_samples,A,B,X,Q,R,U,simsteps,N); 
            if blow == 1
                disp('INSIDE CVX AND BLOW LOOP')
                W = blow_size*W;
                [Xn, Pinf] = w_constructCVXBlow(W,A,B,X,Q,R,U,simsteps,N); 
                if iter_count ~=1
                    intr = intersect(Xn0,Xn); 
                    contXn_flag = Xn0.contains(Xn);                    % should be always 1
                    eqXn_flag = isEmptySet(intr-Xn);                   % 0 is good
                end   
            else
                disp('INSIDE CVX AND NO BLOW LOOP')
                W = W; 
                if iter_count ~=1
                    intr = intersect(Xn0,Xn); 
                    contXn_flag = Xn0.contains(Xn);                     % should be always 1
                    eqXn_flag = isEmptySet(intr-Xn);                    % 0 is good
                end  
            end

        else
            %% Calculating Terminal Set and Best Possible W Bounds 
            disp('INSIDE LRBF LOOP')
            [W, Xn, Pinf, conf_possible] = w_constructGauss(w_samples, conf, nx,nu,A,B,C,D,b,Q,R,U,N,x_0,X,simsteps,options); 
            conf_array(mc,iter_count) = conf_possible; 
        end

        %% Include a piece of code here that would check the probability of failure 
        [prob_fail(iter_count, mc)] = monte_carloSimGauss(W,trueMu,trueStd,nx);

        %% MPC Solving
        dim_t = size(C,1)*N + size(Xn.A,1);          
        polW(iter_count,mc) = W;  
        [capA, capE, capB, capC, capD, Aw_batch, Bu_batch, A_batch] = obtain_matR(A,B,C,D,Xn,nx,nu,N,dim_t);
        matF = capC*capB + capD; 
        matG = capC*capE; 
        matH = -capC*capA; 
        mat_c = [kron(ones(N,1),b); Xn.b];

        %% Starting Closed Loop Computations
        for i=1:simsteps
            constraints = []; 
            Hs=[]; hs =[]; 
            for k = 1:N
                polS = W; 
                Hs_ol = polS.A;  hs_ol = polS.b; 
                Hs = blkdiag(Hs, Hs_ol); hs = [hs; hs_ol];
            end  
           
            %% Creating Open Loop Optimization Variables for MPC 
            dim_a = size(Hs,1); 
            M = sdpvar(nu*N,nx*N); 
            for j=1:nu*N
                for k = 2*j-1:nx*N
                        M(j,k) =0;
                end
            end
            v = sdpvar(nu*N,1); 
            Z = sdpvar(dim_a,dim_t);     

            epsilon = sdpvar(dim_t,1);
            for kk = 1:N
                epsilon(size(C,1)*kk -1: size(C,1)*kk,1) = 0;  % no slacks for input constraints
            end
           %% Solving the Optimization Problem 
            x_pred = sdpvar(nx,N+1);
            x_pred(:,1) = x_cl; 
            cost_state = (x_pred(:,1)-x_ref)'*Q*(x_pred(:,1)-x_ref); 
            
            for k=1:N
                x_pred(:,k+1)  = A*x_pred(:,k) + B*v(1+(k-1)*nu:k*nu,1);
                if k~=N
                    cost_state = cost_state + (x_pred(:,k)-x_ref)'*Q*(x_pred(:,k)-x_ref);
                end
            end
            
            cost_state = cost_state + (x_pred(:,N+1)-x_ref)'*Pinf*(x_pred(:,N+1)-x_ref)...
                         + Wslack*(epsilon')*epsilon*soft_flg;                              

            constraints = [constraints; matF*v + Z'*hs <= mat_c + matH*x_cl + soft_flg*epsilon];
            constraints = [constraints; Z>=0; epsilon >=0];
            constraints = [constraints, matF*M + matG == Z'*Hs];

            obj_ol = v'*kron(R,N)*v + cost_state; 
            diagn=solvesdp(constraints, obj_ol, options);

            %% Obtaining Closed Loop Parameters 
            v_hor = double(v);
            u_cl = v_hor(1:nu,:);                                                     
            eps_cl = double(epsilon);
            eps_cl = eps_cl(1:size(C,1),1);         
            %% Actual System Simulation in Closed Loop 
            flg = 0; 
            %%% REJECTION SAMPLING 
            while flg == 0
                w = trueStd.* randn(nx, 1) + trueMu;
                if gaussPol.contains(w) == 1
                    flg = 1; 
                else
                    flg = 0;
                end
            end           

            w_samples = [w_samples, w]; 
            yalmip 'clear'
            cost_iter(mc,iter_count) = cost_iter(mc,iter_count) + (x_cl-x_ref)'*Q*(x_cl-x_ref)...
                                               + u_cl'*R*u_cl + soft_flg*eps_cl'*Wslack*eps_cl; 
           
            x_cl = A*x_cl + B*u_cl +  w;  

        end

        cost_iter(mc,iter_count) = cost_iter(mc,iter_count) + (x_cl-x_ref)'*Q*(x_cl-x_ref);

        %% Iteration ends. 
        iter_count
        Xn0 = Xn; 

    end

end

%%% Takes time to get to this point
keyboard

%% Getting the results
% Depending on the chosen method, plot the mean cost_iter for each iteration across all monte carlo runs
% Also plot the probability of failures
