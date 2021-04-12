%% Get probability of failure: Uniform case

function prob_fail = prob_fail_Uniform(max_w_samples,wub_true,conf,nsamples)

    nx = length(max_w_samples);
    temp = min(max_w_samples/((1-conf)^(1/nsamples)),wub_true);

    prob_fail = 1 - prod(temp)/wub_true^nx;

end

