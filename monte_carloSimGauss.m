%% Monte Carlo Simulations with in each iteration: Gaussian 
% Monimoy Bujarbaruah
% Akhil Shetty 
%%
% This code will give an idea of probability of failure 

function [prob_fail] = monte_carloSimGauss(W,trueMu,trueStd,nx)

    gaussPol = Polyhedron('lb',trueMu-3*trueStd,'ub',trueMu+3*trueStd);

    flg = 0; mont_count = 1000; flgProb = zeros(mont_count,1);
    
    %%% REJECTION SAMPLING 
    for i = 1:mont_count
        while flg == 0 
            w = trueStd.* randn(nx, 1) + trueMu;
            arr1 = w<=trueMu+3*trueStd;
            arr2 = w>=trueMu-3*trueStd; 
            s1 = sum(arr1); s2 = sum(arr2); 

            if s1==2 && s2==2       
                flg = 1; 
            else
                flg = 0;
            end
        end
         flgProb(i,1) = W.contains(w);                                            % check belonging to our constructed polytope
         flg = 0; 
    end

    prob_fail = (mont_count - sum(flgProb))/mont_count;             % Return empirical failure probability 

end