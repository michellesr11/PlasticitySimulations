function plast_dot = plasticity2(t, var, CI, params)
    nL = var(1); nT = var(2); % nT are ready to release, nL are not, nE are empty and undocked
    p = var(3);
    p0 = CI(3);
    
    b1 = params(1); b2 = params(2); 
    tau_f = params(3); a_f = params(4);
    k1r = params(5); k2r = params(6);
    m1 = params(7); m2 = params(8); 
    
    % k1 and k2 are calcium dependent
    k1 = k1r + m1 * (p-p0);
    k2 = k2r + m2 * (p-p0);
    
    % n subdivided into 2 different vesicle pools where nT+nL+nE = 1
    nL_dot = -(b1 + k1 + k2) * nL + (b2 - k1) * nT + k1 * 1; %because Ntot = 1;
    nT_dot = k2 * nL - b2 * nT;
    
    % change in p only defined by calcium concentration
    recovery = (p0-p)/tau_f;
    p_dot = recovery;
    
    plast_dot = [nL_dot nT_dot p_dot];
end
