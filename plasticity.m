
function plast_dot = plasticity(t, var, CI, tj, params)
    n = var(1); p = var(2);
    p0 = CI(2);
    stimuli = length(tj);
    tau_r = params(1);
    tau_f = params(2);
    a_f = params(3);
    %sigma = 0.00001;
    
    % depression
    replenish = (1-n)/tau_r;
    depletion = 0; %unsure
    for i = 1:stimuli 
        %depletion = depletion + ((-1/(sigma*sqrt(2*pi)))*exp(1)^(-(t-tj(i))^2/(2*sigma^2)))*p*n; % una normal muy flaquita aproxima una delta
        %depletion = depletion - normpdf(t,round(tj(i),4),sigma)*p*n; % una normal muy flaquita aproxima una delta
        depletion = depletion + ((t == tj(i)) * p * n);
    end
    n_dot = replenish + depletion;
    
    % facilitation
    recovery = (p0-p)/tau_f;
    calcium = 0;
    for i = 1:stimuli  
        %calcium = calcium + ((1/(sigma*sqrt(2*pi)))*exp(1)^(-(t-tj(i))^2/(2*sigma^2)))*a_f*(1-p); 
        %calcium = calcium + normpdf(t,round(tj(i),4),sigma)*a_f*(1-p); 
        calcium = calcium + ((t == tj(i)) * a_f *(1-p)); 
    end
    p_dot = recovery + calcium;
    
    plast_dot = [n_dot p_dot];
end
