function plast_dot = plasticity5(t, var, CI, params)
    nL = var(2); nT = var(3); % nT are ready to release, nL are not, nE are empty and undocked
    p = var(4); b2 = var(1);
    p0 = CI(4); b2r = CI(1);
    
    
    b1 = params(1); 
    tau_f = params(2);
    k1r = params(3); k2r = params(4); tau_b = params(5); 
    m1 = params(6); m2 = params(7); 
    
    % k1 and k2 are calcium dependent  
    k1 = k1r + m1 * (p-p0);
    k2 = k2r + m2 * (p-p0);
    
    % make b2 change on stimulation
    b2_dot = (b2r-b2)/tau_b;
    
    % n subdivided into 2 different vesicle pools where nT+nL+nE = 1
    nL_dot = -(b1 + k1 + k2) * nL + (b2 - k1) * nT + k1 * 1; %because Ntot = 1;
    nT_dot = k2 * nL - b2 * nT;
    
    % change in p only defined by calcium concentration
    recovery = (p0-p)/tau_f;
    p_dot = recovery;
    
    plast_dot = [b2_dot nL_dot nT_dot p_dot];
end