%% Terminal Components with Inflated CVX Hull of Seen Disturbances
% Monimoy Bujarbaruah
% Akhil Shetty

function [Xn, Pinf] = w_constructCVXBlow(W,A,B, X,Q,R,U,simsteps,N)

    [Finf, Pinf] = dlqr(A,B,Q,R);  
    model = LTISystem('A',A-B*Finf,'B',B); 

    Hx = X.A; Hu = U.A; 
    hx = X.b; hu = U.b; 

    Hxaut = [Hx; -Hu*Finf]; hxaut = [hx; hu]; 

    S = Polyhedron('A',Hxaut,'b',hxaut); 
    Xn = MRPI_Fin(model, S, W, simsteps,N);   

 end