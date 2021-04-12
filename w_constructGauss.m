%% Uncertainty Support Construction for Truncated Gaussian Case: Bootsrap
% Monimoy Bujarbaruah
% Akhil Shetty
%%
% Sets a high confidence value. Uses data. Constructs support. Checks
% feasibility of MPC problem from x_0. If infeasible, lowers confidence and
% repeats. if feasible, this support is rolled out to the main code to
% complete an iteration. Then data is collected and process is repeated

function [W, Xn, Pinf, conf_possible] = w_constructGauss(w_samples, conf, nx,nu,A,B,C,D,b,Q,R,U,N,x_0,X,simsteps,options)

    bsSize = 1000;                                                              % Bootstrap copies 
    wbs = zeros(nx,size(w_samples,2),bsSize);                 

    %% Start Bootstrap here
    meanEmp = mean(w_samples,2);
    stdEmp = std(w_samples,0,2); 

    for j = 1: bsSize
        for i = 1: size(w_samples,2)
            entr = randi(size(w_samples,2)); 
            wbs(:,i,j) = w_samples(:,entr);
        end
    end

    meanBatch  = zeros(nx,bsSize);
    stdBatch      = zeros(nx,bsSize); 

    for j = 1:bsSize
        meanBatch(:,j) = mean(wbs(:,:,j),2);
        stdBatch(:,j) =  std(wbs(:,:,j),0,2);
    end

    %%%%%% Take the confidence set %%%%%%%%%
    minMu = zeros(nx,1);
    maxMu = zeros(nx,1);
    maxStd = zeros(nx,1);
    
    %%
    feas_conf = 0;
    count = 0; count2 = 0; 

    while feas_conf == 0     

        if conf >= 0.001
            conf = conf/(2^count);    
            for i = 1:nx
                minMu(i,1) = quantile(meanBatch(i,:),(1-conf)/2);
                maxMu(i,1) = quantile(meanBatch(i,:),1-(1-conf)/2);
                maxStd(i,1) = quantile(stdBatch(i,:),1-(1-conf)/2);
            end            
            w_lb = minMu - 3.08*maxStd;   
            w_ub = maxMu + 3.08*maxStd;                       
            W = Polyhedron('lb',w_lb, 'ub', w_ub);               
        else  
            w_lb = meanEmp - (3.08-0.1*count2)*stdEmp;   
            w_ub = meanEmp + (3.08-0.1*count2)*stdEmp;        
            W = Polyhedron('lb',w_lb, 'ub', w_ub);            
            count2 = count2 + 1; 
        end

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

end