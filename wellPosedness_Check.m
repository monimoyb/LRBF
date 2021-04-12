%% Problem's Feasibility Check. If not, then no point of doing anything
% Monimoy Bujarbaruah
% Akhil Shetty
%% 
clear all
close all
clc
yalmip 'clear'

%% Form W based on uniform or clipped Gaussian case
[A,B,C,D,b,X,U,nx,nu,wub_true,wlb_true, x_0,Q,R,N, trueMu, trueStd, x_ref,simsteps] = sys_load(); 
unif_flg = 1;  % 1 for uniform, 0 for clipped gaussian

if unif_flg == 1
    W = Polyhedron('lb',wlb_true*ones(nx,1),'ub',wub_true*ones(nx,1));
else
    W = Polyhedron('lb',trueMu-3*trueStd,'ub',trueMu+3*trueStd);
end

%% Terminal set empty flag
[Finf, Pinf] = dlqr(A,B,Q,R);  
model = LTISystem('A',A-B*Finf,'B',B); 
Hx = X.A; Hu = U.A; 
hx = X.b; hu = U.b; 
Hxaut = [Hx; -Hu*Finf]; hxaut = [hx; hu]; 
S = Polyhedron('A',Hxaut,'b',hxaut); 
Xn = MRPI_Fin(model, S, W, simsteps,N);   
term_flg = isEmptySet(Xn);                                      % Must be 0
keyboard
%%% No point going ahead if terminal set is empty. 

%% Feasibility check 
options = sdpsettings('solver','gurobi','verbose',0);
dim_t = size(C,1)*N + size(Xn.A,1);                          
[capA, capE, capB, capC, capD, Aw_batch, Bu_batch, A_batch] = obtain_matR(A,B,C,D,Xn,nx,nu,N,dim_t);
matF = capC*capB + capD; 
matG = capC*capE; 
matH = -capC*capA; 
mat_c = [kron(ones(N,1),b); Xn.b];
constraints = []; 
Hs=[]; hs =[]; 

for k = 1:N
    polS = W; 
    Hs_ol = polS.A;  hs_ol = polS.b; 
    Hs = blkdiag(Hs, Hs_ol); hs = [hs; hs_ol];
end
dim_a = size(Hs,1);     
      
M = sdpvar(nu*N,nx*N); 
for j=1:nu*N
    for k = 2*j-1:nx*N
        M(j,k) =0;
    end
end
v = sdpvar(nu*N,1); 
Z = sdpvar(dim_a,dim_t);     
      
constraints = [constraints; matF*v + Z'*hs <= mat_c + matH*x_0];
constraints = [constraints; Z>=0];
constraints = [constraints, matF*M + matG == Z'*Hs];
diagn=solvesdp(constraints, [], options);

%% Flag
feas = diagn.problem;   % must be 0
  
