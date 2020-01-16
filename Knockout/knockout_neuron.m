%%
function [individual_prime, config_prime] = knockout_neuron(individual, neuron_num, config)

    % convert from sparse to full matrix
    W = full(cell2mat(individual.W));
    
    % remove the neuron's internal connections
    W(neuron_num,:) = [];
    W(:,neuron_num) = [];
    
    % remove the neuron's external connections
    Win  = cell2mat(individual.input_weights);
    Wout = individual.output_weights;

    Win(neuron_num,:) = [];
    Wout(neuron_num,:) = [];
    
    config_prime = config;
    
    % keep track of which neurons have been knocked out
    if isfield(config, 'knocked_out_neurons')
        config_prime.knocked_out_neurons = [config.knocked_out_neurons neuron_num];
    else
        config_prime.knocked_out_neurons = [neuron_num];
    end

    
    % write changes to the new individual
    individual_prime                    = individual;
    individual_prime.W{1,1}             = sparse(W);
    individual_prime.input_weights{1,1} = Win;
    individual_prime.output_weights     = Wout;
    individual_prime.nodes              = individual.nodes - 1;
    individual_prime.total_units        = individual.total_units - 1;
    
