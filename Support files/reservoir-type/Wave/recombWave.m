%% Infection phase
function loser = recombWave(winner,loser,config)

% params - input_scaling, leak_rate,
W= winner.input_scaling(:);
L = loser.input_scaling(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
L(pos) = W(pos);
loser.input_scaling = reshape(L,size(loser.input_scaling));

W= winner.leak_rate(:);
L = loser.leak_rate(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
L(pos) = W(pos);
loser.leak_rate = reshape(L,size(loser.leak_rate));

% params - 
W= winner.time_period(:);
L = loser.time_period(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));
L(pos) = W(pos);
loser.time_period = reshape(L,size(loser.time_period));

W= winner.wave_speed(:);
L = loser.wave_speed(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));
L(pos) = W(pos);
loser.wave_speed = reshape(L,size(loser.wave_speed));

W= winner.damping_constant(:);
L = loser.damping_constant(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));
L(pos) = W(pos);
loser.damping_constant = reshape(L,size(loser.damping_constant));

W= winner.boundary_conditions(:);
L = loser.boundary_conditions(:);
pos = randperm(length(L),ceil(config.rec_rate*length(L)));
L(pos) = W(pos);
loser.boundary_conditions = reshape(L,size(loser.boundary_conditions));

% cycle through sub-reservoirs
for i = 1:config.num_reservoirs
    
    % input weights
    W= winner.input_weights{i}(:);
    L = loser.input_weights{i}(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
    L(pos) = W(pos);
    loser.input_weights{i} = reshape(L,size(loser.input_weights{i}));
       
end

% for output weights
if config.evolve_output_weights
    W= winner.output_weights(:);
    L = loser.output_weights(:);
    pos = randperm(length(L),ceil(config.rec_rate*length(L)));         
    L(pos) = W(pos);
    loser.output_weights = reshape(L,size(loser.output_weights));
end
