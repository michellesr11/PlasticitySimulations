
ISIs = [0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1];
%% CLASSIC: residual calcium and use-dependent vesicle depletion
params = [12 0.24 0.26]; % [tau_r, tau_f, a_f]
CI = [1 0.1]; 
maxtime = 5; % integration period (in seconds)
step_size = 0.0001; % integration step
start_time = 0.1; % time at which to deliver the first impulse

ISI_PIs = NaN(numel(ISIs),1); % store plasticity index for each ISI
for i = 1:numel(ISIs)
    stim_times = [start_time start_time+ISIs(i)]; % stimulus times
    [t_steps,f_steps] = simulatePlast(1,CI,params,stim_times,step_size,maxtime);
    % get the amplitude of the PSP as n*p
    amps = PSPamps(t_steps,f_steps,stim_times);
    % normalise them
    amps_norm = amps / amps(1);
    ISI_PIs(i) = amps_norm(2);
end

%% Lin (Neher): residual calcium and two vesicle pools
params = [3.5 0.4 0.2 0.3 0.5 0.5 1.5 0.01]; %[b1 b2 tau_f a_f k1r k2r m1 m2]
% get steady state
nL0 = (params(2)*params(5))/((params(1)*params(2))+params(5)*(params(2)+params(6)));
nT0 = (params(5)*params(6))/((params(1)*params(2))+params(5)*(params(2)+params(6)));
CI = [nL0 nT0 0.1]; % [nL0 nT0 p0] 
maxtime = 5;
step_size = 0.0001; % integration step
start_time = 0.1; % time at which to deliver the first impulse

ISI_PIs = NaN(numel(ISIs),1); % store plasticity index for each ISI
for i = 1:numel(ISIs)
    stim_times = [start_time start_time+ISIs(i)]; % stimulus times
    [t_steps,f_steps] = simulatePlast(2,CI,params,stim_times,step_size,maxtime);
    % get the amplitude of the PSP as n*p
    amps = PSPamps(t_steps,f_steps,stim_times);
    % normalise them
    amps_norm = amps / amps(1);
    ISI_PIs(i) = amps_norm(2);
end

%% Delayed retrieval: residual calcium and two vesicle pools with delayed retrieval
params = [0.05 0.19 1 1 1.2 35 35 0.30 5.5 0.67]; % [b1 tau_f k1r k2r tau_b m1 m2 a_f y del]
b20 = 0.1;
% get steady state
nT0 = (params(3)*params(4))/((params(1)*b20)+params(3)*(b20+params(4)));
nL0 = (b20*params(3))/((params(1)*b20)+params(3)*(b20+params(4)));
CI = [b20 nL0 nT0 0.1]; % [b2r nL0 nT0 p0]
maxtime = 5;
step_size = 0.0001; % integration step
start_time = 0.1; % time at which to deliver the first impulse

ISI_PIs = NaN(numel(ISIs),1); % store plasticity index for each ISI
for i = 1:numel(ISIs)
    stim_times = [start_time start_time+ISIs(i)]; % stimulus times
    [t_steps,f_steps] = simulatePlast(3,CI,params,stim_times,step_size,maxtime);
    % get the amplitude of the PSP as n*p
    amps = PSPamps(t_steps,f_steps,stim_times);
    % normalise them
    amps_norm = amps / amps(1);
    ISI_PIs(i) = amps_norm(2);
end

%%
function [t_steps,f_steps] = simulatePlast(model,CI,params,stim_times,step_size,int_time)

numVar = numel(CI);
tj = stim_times;
h = step_size;

% consider preallocating size of t_steps and f_steps; the size would be numel(0:h:int_time)
if model == 1 % classic
    for int_period = 1:numel(tj)+1
        if int_period == 1
            t_steps_temp = [0:h:tj(int_period)];
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity(t, var, CI, tj, params), numVar, CI, t_steps_temp);
            p_period = f_steps_period(end,2); n_period = f_steps_period(end,1);
            t_steps = t_steps_period;
            f_steps = f_steps_period;
            
        elseif int_period == numel(tj)+1
            t_steps_temp = [tj(int_period-1)+h:h:int_time];
            newCI = [n_period-p_period*n_period p_period+params(3)*(1-p_period)];
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity(t, var, CI, tj, params), numVar, newCI, t_steps_temp);
            t_steps = [t_steps t_steps_period];
            f_steps = vertcat(f_steps, f_steps_period);
            
        else
            t_steps_temp = [tj(int_period-1)+h:h:tj(int_period)];
            newCI = [n_period-p_period*n_period p_period+params(3)*(1-p_period)];
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity(t, var, CI, tj, params), numVar, newCI, t_steps_temp);
            p_period = f_steps_period(end,2); n_period = f_steps_period(end,1);
            t_steps = [t_steps t_steps_period];
            f_steps = vertcat(f_steps, f_steps_period);
        end
    end
    
elseif model == 2   % lin - neher 
    for int_period = 1:numel(tj)+1
        if int_period == 1
            t_steps_temp = [0:h:tj(int_period)];
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity2(t, var, CI, params), numVar, CI, t_steps_temp);
            p_period = f_steps_period(end,3); nT_period = f_steps_period(end,2);
            t_steps = t_steps_period;
            f_steps = f_steps_period;
        elseif int_period == numel(tj)+1
            t_steps_temp = [tj(int_period-1)+h:h:int_time];
            newCI = [f_steps_period(end,1) nT_period-p_period*nT_period p_period+params(4)*(1-p_period)];
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity2(t, var, CI, params), numVar, newCI, t_steps_temp);
            t_steps = [t_steps t_steps_period];
            f_steps = vertcat(f_steps, f_steps_period);       
        else
            t_steps_temp = [tj(int_period-1)+h:h:tj(int_period)];
            % define state after release
            newCI = [f_steps_period(end,1) nT_period-p_period*nT_period p_period+params(4)*(1-p_period)];
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity2(t, var, CI, params), numVar, newCI, t_steps_temp);
            p_period = f_steps_period(end,3); nT_period = f_steps_period(end,2);
            t_steps = [t_steps t_steps_period];
            f_steps = vertcat(f_steps, f_steps_period);
        end      
    end

elseif model == 3   % retrieval
    tj_orig = tj; ret = params(end);
    tj = sort(horzcat(tj,tj+ret));
    for int_period = 1:numel(tj)+1
        if int_period == 1 % baseline before any stimulation
            t_steps_temp = [0:h:tj(int_period)];
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity5(t, var, CI, params), numVar, CI, t_steps_temp);
            p_period = f_steps_period(end,4); nT_period = f_steps_period(end,3); b2_period = f_steps_period(end,1); p_max = f_steps_period(1,4);
            t_steps = t_steps_period;
            f_steps = f_steps_period;
        elseif int_period == numel(tj)+1
            t_steps_temp = [tj(int_period-1):h:int_time];
            % define state after release
            if ~isempty(find((tj(int_period-1) == tj_orig+ret)))
                newCI = [b2_period+params(9)*(p_max) f_steps_period(end,2) nT_period p_period];
            else
                newCI = [b2_period f_steps_period(end,2) nT_period-p_period*nT_period p_period+params(8)*(1-p_period)];
            end
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity5(t, var, CI, params), numVar, newCI, t_steps_temp);
            t_steps = [t_steps t_steps_period];
            f_steps = vertcat(f_steps, f_steps_period);
        elseif int_period == 2
            t_steps_temp = [tj(int_period-1):h:tj(int_period)];
            % define state after release
            newCI = [b2_period f_steps_period(end,2) nT_period-p_period*nT_period p_period+params(8)*(1-p_period)];
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity5(t, var, CI, params), numVar, newCI, t_steps_temp);
            p_period = f_steps_period(end,4); nT_period = f_steps_period(end,3); b2_period = f_steps_period(end,1); p_max = f_steps_period(1,4);
            t_steps = [t_steps t_steps_period];
            f_steps = vertcat(f_steps, f_steps_period);
        else
            t_steps_temp = [tj(int_period-1):h:tj(int_period)];
            % define state after release
            if ~isempty(find((tj(int_period-1) == tj_orig+ret)))
                newCI = [b2_period+params(9)*(p_max) f_steps_period(end,2) nT_period p_period];
            else
                newCI = [b2_period f_steps_period(end,2) nT_period-p_period*nT_period p_period+params(8)*(1-p_period)];
            end
            [t_steps_period, f_steps_period] = RungeKutta4(@(t,var) plasticity5(t, var, CI, params), numVar, newCI, t_steps_temp);
            p_period = f_steps_period(end,4); nT_period = f_steps_period(end,3); b2_period = f_steps_period(end,1); p_max = f_steps_period(1,4);
            t_steps = [t_steps t_steps_period];
            f_steps = vertcat(f_steps, f_steps_period);
        end
    end
    
end

end

%%
function amps = PSPamps(t_steps,f_steps,tj)

f_steps = f_steps(:,end-1:end); % make sure the last two variables are vesicle pool and probability
amps = zeros(numel(tj),1);
for st = 1:numel(tj)
    tjIdx = find(round(t_steps,4) == round(tj(st),4),1);
    amps(st) = f_steps(tjIdx,1) * f_steps(tjIdx,2);
end

end