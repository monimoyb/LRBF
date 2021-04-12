%% Defining System and Constraint Matrices 
% Monimoy Bujarbaruah and Akhil Shetty

function [A,B,C,D,b,X,U,nx,nu,wub_true,wlb_true, x_0, Q,R,N, trueMu, trueStd, x_ref,simsteps] = sys_load()

    %% Considering two states and one scalar input 
    A = [1.2, 1.3; 0, 1.5]; 
    B = [0;1];   
    nx = size(A,2); nu = size(B,2); 
    %% Weights and horizon
    Q =  10*eye(nx);
    R =   2*eye(nu);
    N = 4;
    simsteps = 20;                                                                               
    %% Constraints 
    % Considering constraints of the form -a<=x(i)<=a and -ulb<=u<=uub
    % and expressing in Cx+Du <=b format 
    C = [1 0; -1 0; 0 1; 0 -1; 0 0; 0 0]; 
    D = [0; 0; 0; 0; 1; -1]; 
    b = [30;30;30;30;40;40]; 
    X = Polyhedron('A',C(1:4,:),'b',b(1:4,:));
    U = Polyhedron('A',D(5:6,:),'b',b(5:6,:));    
    %% Defining disturbance  bounds 
    wub_true = 3;                                  % Upper bound 
    wlb_true = -3;                                  % Lower bound 
    %% Moments of clipped Gaussian
    trueMu =  zeros(nx,1);
    trueStd = ones(nx,1); 
    %% Starting state and reference
    x_0 = [0; 0]; 
    x_ref = [27; 27];
    
end