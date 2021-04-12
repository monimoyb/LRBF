%% Support Construction using CVX Hull of Seen Disturbance Samples
% Monimoy Bujarbaruah
% Akhil Shetty

function [W, Xn, Pinf] = w_constructCVX(w_samples,A,B,X,Q,R,U,simsteps,N)

    W = Polyhedron('V',w_samples');         
    [Finf, Pinf] = dlqr(A,B,Q,R);  
    model = LTISystem('A',A-B*Finf,'B',B); 

    Hx = X.A; Hu = U.A; 
    hx = X.b; hu = U.b; 
    Hxaut = [Hx; -Hu*Finf]; hxaut = [hx; hu]; 

    S = Polyhedron('A',Hxaut,'b',hxaut); 
    Xn = MRPI_Fin(model, S, W, simsteps,N);   

 end
