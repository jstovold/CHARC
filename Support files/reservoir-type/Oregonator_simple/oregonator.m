%% Oregonator code
% Author: Dr James Stovold <james.stovold@york.ac.uk>
% ======
% 
% This code implements a basic forward Euler integration on the photosensitive (Ru-catalysed)
% variant of the Oregonator model: http://www.scholarpedia.org/article/Oregonator 
%  
% Example call: oregonator([1 0 1 1; 0 1 1 1], [0.05, 0, 0.05, 0.1; 0, 0.07, 0.01, 0.1], config)
%
% ======
%
% TODO:
%  - add GPU support to step function
%  - add trace outputs
%  - add time-varying inputs



function [output_data, output] = oregonator(bitmatrix, phimatrix, input_data, config)
    config = getDefaultParams(config);
    reactor = createReactor(bitmatrix, phimatrix, config);
    state = initReactor(reactor, config);
    updateDisplay(state, config);
    
    maxPhiValue = max(phimatrix, [], 'all');
    
    output_data = runReservoir(state, input_data, maxPhiValue, config);
    
    processed_data = processOutputs(output_data, maxPhiValue, config);
    
    output = processed_data;
%     for t = 1:config.timesteps
%         
%         state = oreg_step(state, config);
%         if mod(t, config.stride) == 0
%             updateDisplay(state, config);
%         end
%     end

end

%% default parameters
function config = getDefaultParams(config)
    config.height      = 110;
    config.width       = 200;
    config.timesteps   = 100000;
    config.stride      = 100;
    config.displayLive = true;
    config.displayZoom = 4.0;
    
    config.epsilon     = 0.0243;
    config.grid_2      = 0.0625; % 0.25^2
    config.diff_coeff  = 0.45;
    config.f           = 1.4;
    config.q           = 0.002;
    config.speed       = 0.0005;
    config.phi_active  = 0.054;    % normal, excitable
    config.phi_passive = 0.0975;   % passive, excitable
    
    config.waveThreshold      = 0.6;
    config.numSamplesPerEpoch = 10;
    
    config.u_idx = 1;
    config.v_idx = 2;
    config.p_idx = 3;
    config.b_idx = 4;

    config.input_x = 75;   % eventually these will be converted to lists for multiple inputs
    config.input_y = 75;   % but for now...
    
    config.trace_x = 180;
    config.trace_y = 30;
    
    
    config.vesicle_radius = 25;
%     config.phi_active  = 0.0758;   % normal, wavelet / sub-ex
%     config.phi_passive = 0.10015;  % passive, wavelet / sub-ex
end

%% display update
function img = updateDisplay(state, config)

    % take state and produce image for display
    % display image in window?
  
    R = state(:,:,config.u_idx);
    G = state(:,:,config.p_idx);
    B = state(:,:,config.v_idx);

    R = min(R + state(:,:,config.b_idx), ones(config.height, config.width, 1));
    
    img = cat(3, R, G, B);
    if config.displayLive
        if config.displayZoom ~= 1.0
            img = imresize(img, [config.height * config.displayZoom, ...
                                 config.width * config.displayZoom]);
        end
        imshow(img);
    end
    

end

%% create reactor from matrix
function reactor = createReactor(bitMatrix, phiMatrix, config)
    
    % takes the bit matrix and uses it to define the corresponding connected vesicle reactor
    % config needs size of vesicle
    
    reactor = zeros(config.height, config.width, 2);
    
    margin_w = 5;
    margin_h = 5;
    
    vesicle_w = config.vesicle_radius * 2;
    vesicle_h = vesicle_w;
    % calculate inter-vesicle distance (IVD)
    
    crossover_dist = 3;
    IVD = vesicle_w - crossover_dist;
    
    % calculate size of vesicle `domain' (i.e the area in the reactor within which we will place the vesicles)    
    
    m_height = size(bitMatrix,1);
    m_width  = size(bitMatrix,2);
    
    d_width  = m_width * IVD + (margin_w * 2);
    d_height = m_height * IVD + (margin_h * 2);
    
    if d_width > config.width
        error(strcat('Vesicle domain: ', string(d_width), ' is too wide for reactor: ', string(config.width)));
    end
    
    if d_height > config.height
        error('Vesicle domain is too high for reactor');
    end
    
    % distribute vesicle centres across vesicle domain, at IVD intervals
    
    % first vesicle centre: (x, y)
    v_centre_x = margin_w + config.vesicle_radius;
    v_centre_y = margin_h + config.vesicle_radius;
    
    % use bitmatrix to draw vesicles around desired vesicle centres 
    for y = 1:m_height
        v_centre_x = margin_w + config.vesicle_radius;
        for x = 1:m_width
           if bitMatrix(y, x) == 1
               % draw vesicle
               
               reactor = drawVesicle(reactor, v_centre_x, v_centre_y, config.vesicle_radius);
               
           end
           v_centre_x = v_centre_x + IVD;
        end
        v_centre_y = v_centre_y + IVD;
    end
    
    % remove vesicle boundary from adjacent vesicles
    
    v_centre_y = margin_h + config.vesicle_radius;
    for y = 1:m_height
        v_centre_x = margin_w + config.vesicle_radius;
        for x = 1:m_width
            left = false;
            top = false;

            % check if vesicle to the left or top
            if x > 1 
                % check left
                if bitMatrix(y, x - 1) == 1 && bitMatrix(y, x) == 1
                    left = true;
                end
            end
            if y > 1
                if bitMatrix(y - 1, x) == 1 && bitMatrix(y, x) == 1
                    top = true;
                end
            end
            % remove overlapping boundary if present
            phi = 0;
            if bitMatrix(y, x) == 1
                phi = phiMatrix(y, x);
            end
            
            if left
                reactor = makeConnection(reactor, v_centre_x, v_centre_y, v_centre_x - IVD, v_centre_y, config.vesicle_radius);
                reactor = makeConnection(reactor, v_centre_x - IVD, v_centre_y, v_centre_x, v_centre_y, config.vesicle_radius);
                
                reactor = setExcitability(reactor, phi, v_centre_x, v_centre_y, v_centre_x - IVD, v_centre_y, config.vesicle_radius);
            end

            if top
                reactor = makeConnection(reactor, v_centre_x, v_centre_y, v_centre_x, v_centre_y - IVD, config.vesicle_radius);
                reactor = makeConnection(reactor, v_centre_x, v_centre_y - IVD, v_centre_x, v_centre_y, config.vesicle_radius);
                
                reactor = setExcitability(reactor, phi, v_centre_x, v_centre_y, v_centre_x, v_centre_y - IVD, config.vesicle_radius);
            end
            
            if ~left && ~top
                reactor = setExcitability(reactor, phi, v_centre_x, v_centre_y, 0, 0, config.vesicle_radius);
            end
            
            
            v_centre_x = v_centre_x + IVD;
        end
        v_centre_y = v_centre_y + IVD;
    end
    
    
    % done.
    
    % later: set phi for each vesicle

end

%% set the phi for the given vesicle
function reactor = setExcitability(reactor, phi, x0, y0, x1, y1, radius)

    % find all pixels within the vesicle
    % if there are connected vesicles, 
    %   find the locus between (x0, y0) and (x1, y1)
    %   this locus forms the boundary between the two vesicles
    %   remove pixels from array that are closer to (x1, y1) than (x0, y0)
    % set phi for all pixels in array
    
    hasConnection = false;
    if x1 > 0 && y1 > 0
        hasConnection = true;
    end
    
    % draw box around the vesicle:
    for y = (y0 - radius):(y0 + radius)
        for x = (x0 - radius):(x0 + radius)
            if euclDist(x, y, x0, y0) <= radius

                if hasConnection
                    if euclDist(x, y, x0, y0) < euclDist(x, y, x1, y1)
                        reactor(y, x, 2) = phi;
                    end
                else
                    reactor(y, x, 2) = phi;
                end
                
            end
        end
    end
    
end

function dist = euclDist(x_0,y_0,x_1,y_1) 

    dx = x_0 - x_1;
    dy = y_0 - y_1;

    dist = sqrt(dx * dx + dy * dy);

end


%% form connection between vesicles
function reactor = makeConnection(reactor, x0, y0, x1, y1, radius)

    % remove all points on the circle which are < radius from the other centre-point
    
    f = 1 - radius;
    ddF_x = 0;
    ddF_y = -2 * radius;
    x = 0;
    y = radius;

  
    
    if euclDist(x0, y0 + radius, x1, y1) < radius
        reactor(y0 + radius, x0, 1) = 0.0;
    end
    
    if euclDist(x0, y0 - radius, x1, y1) < radius
        reactor(y0 - radius, x0, 1) = 0.0;
    end
    
    if euclDist(x0 + radius, y0, x1, y1) < radius
        reactor(y0, x0 + radius, 1) = 0.0;
    end
    
    if euclDist(x0 - radius, y0, x1, y1) < radius
        reactor(y0, x0 - radius, 1) = 0.0;
    end
    
    
    
    while (x < y)
        if f >= 0
            y     = y - 1;
            ddF_y = ddF_y + 2;
            f     = f + ddF_y;
        end
        
        x     = x + 1;
        ddF_x = ddF_x + 2;
        f     = f + ddF_x + 1;
        
        if euclDist(x0 + x, y0 + y, x1, y1) < radius
            reactor(y0 + y, x0 + x, 1) = 0.0;
        end
        
        if euclDist(x0 - x, y0 + y, x1, y1) < radius
            reactor(y0 + y, x0 - x, 1) = 0.0;
        end
        
        if euclDist(x0 + x, y0 - y, x1, y1) < radius
            reactor(y0 - y, x0 + x, 1) = 0.0;
        end
        
        if euclDist(x0 - x, y0 - y, x1, y1) < radius
            reactor(y0 - y, x0 - x, 1) = 0.0;
        end
        
        if euclDist(x0 + y, y0 + x, x1, y1) < radius
            reactor(y0 + x, x0 + y, 1) = 0.0;
        end
        
        if euclDist(x0 - y, y0 + x, x1, y1) < radius
            reactor(y0 + x, x0 - y, 1) = 0.0;
        end
        
        if euclDist(x0 + y, y0 - x, x1, y1) < radius
            reactor(y0 - x, x0 + y, 1) = 0.0;
        end
        
        if euclDist(x0 - y, y0 - x, x1, y1) < radius
            reactor(y0 - x, x0 - y, 1) = 0.0;
        end
    
    end


end

%% draw a vesicle on the reactor
function reactor = drawVesicle(reactor, x0, y0, radius)

    f = 1 - radius;
    ddF_x = 0;
    ddF_y = -2 * radius;
    x = 0;
    y = radius;
    
    reactor(y0 + radius, x0, 1) = 1.0;
    reactor(y0 - radius, x0, 1) = 1.0;
    reactor(y0, x0 + radius, 1) = 1.0;
    reactor(y0, x0 - radius, 1) = 1.0;
    
    
    while (x < y)
        if f >= 0
            y     = y - 1;
            ddF_y = ddF_y + 2;
            f     = f + ddF_y;
        end
        
        x     = x + 1;
        ddF_x = ddF_x + 2;
        f     = f + ddF_x + 1;
        
        reactor(y0 + y, x0 + x, 1) = 1.0;
        reactor(y0 + y, x0 - x, 1) = 1.0;
        reactor(y0 - y, x0 + x, 1) = 1.0;
        reactor(y0 - y, x0 - x, 1) = 1.0;
        
        reactor(y0 + x, x0 + y, 1) = 1.0;
        reactor(y0 + x, x0 - y, 1) = 1.0;
        reactor(y0 - x, x0 + y, 1) = 1.0;
        reactor(y0 - x, x0 - y, 1) = 1.0;

    end

end

%% initialise reactor variables
function state = initReactor(reactor, config)
    
    state = zeros(config.height, config.width, 4);
    
%     state(round(config.height / 2), round(config.width / 2), config.u_idx) = 1.0;
    state(config.input_y, config.input_x, config.u_idx) = 1.0;
    state(:,:,config.p_idx) = reactor(:,:,2); %config.phi_active;
    state(:,:,config.b_idx) = reactor(:,:,1);
    
    
end

function state = stimulateInput(state, amount, config)
    state(config.input_y, config.input_x, config.u_idx) = amount;
end




%% representation
% we need a consistent sampling rate: for how long does each data point last?

function strideLength = getStrideLength(phiValue)
    
    % minimum phi value: 0.03
    % maximum phi value: 0.07 (no waves form)

    % if maximum phi value <= 0.06  ->11000
    % if maximum phi value <= 0.055 -> 9000
    % if maximum phi value <= 0.05  -> 8000
    % if maximum phi value <= 0.04  -> 7000
    
    strideLength = 0;
    
    if phiValue <= 0.04
        strideLength = 7000;
    elseif phiValue <= 0.05
        strideLength = 8000;
    elseif phiValue <= 0.055
        strideLength = 9000;
    else 
        strideLength = 11000;
    end
    
end

function [output_data, state]  = runReservoir(state, input_data, maxPhiValue, config)
    
    trace_y = config.trace_y;
    trace_x = config.trace_x;
    u_idx   = config.u_idx;
    
    % convert magnitude of each data point to frequency of waves?
    % let's assume all data is in the range 0..1
    
    % 1 -> maximum frequency (i.e. minimum time between waves)
    % 0 -> minimum frequency (and so will be dictated by our sampling rate)
    
    % need the mapping from phi value to stride length (i.e. how close together can our waves be?
    
    strideLength = getStrideLength(maxPhiValue);
    numSamples   = config.numSamplesPerEpoch;    % at least 3
    timeLength   = strideLength * numSamples;
    
    output_data  = zeros(1, length(input_data) * timeLength); % currently assuming only 1 output
%     output_data  = zeros(config.numOutputs, length(input_data) * timeLength); % currently assuming only 1 output
    
    t = 1;
    for t_i = 1:length(input_data)
        
        mag = input_data(t_i);
        
        if mag > 0
            [quantity, delay] = convertMagToQuant(mag, timeLength, strideLength);
            if quantity > 0
                for q = 1:quantity
                    % stimulate the inputs
                    % wait \strideLength\ timesteps
                    % repeat q times

                    state = stimulateInput(state, 1.0, config);

                    for i = 1:delay
                        state = oreg_step(state, config);
                        
                        if mod(t, config.stride) == 0
                            updateDisplay(state, config);
                        end

                        output_data(t) = state(trace_y, trace_x, u_idx);
                        
                        t = t + 1;
                    end

                end
            else
                for i = 1:timeLength
                    state = oreg_step(state, config);
                    
                    if mod(t, config.stride) == 0
                        updateDisplay(state, config);
                    end
                    
                    output_data(t) = state(trace_y, trace_x, u_idx);
                    t = t + 1;
                end
            end
        else
            for i = 1:timeLength                
                state = oreg_step(state, config);
                if mod(t, config.stride) == 0
                    updateDisplay(state, config);
                end
                
                output_data(t) = state(trace_y, trace_x, u_idx);
                t = t + 1;

            end
        end
        
    end
    
end

function output_data = processOutputs(outputStream, maxPhiValue, config)

    % mag = freq / sampleRate;
    
    % - if we receive a wave every /sampleRate/ timesteps, we have the maximum frequency (i.e. mag=1)
    % - if we receive no waves in /sampleRate/ timesteps, we have the minimum frequency (i.e. mag=0)
    % - counting the number of waves in multiple rounds of /sampleRate/ timesteps, we can produce a linear mapping
    % between frequency of waves and magnitude in the range 0..1

    % convert frequency of waves to magnitude of data
    strideLength = getStrideLength(maxPhiValue); 
    numSamples   = config.numSamplesPerEpoch;    % at least 3
    timeLength   = strideLength * numSamples;
    threshold    = config.waveThreshold;         % how high does the concentration need to be to register as a wave?
    output_data  = zeros(1, ceil(length(outputStream) / timeLength));
    o = 1;
    for t = 1:timeLength:length(outputStream)
        inWave = false;
        waves  = 0;
        for i = 1:timeLength
            % count number of waves
            u = outputStream(t + i - 1);
            if ~inWave && u > threshold
                % rising edge
                inWave = true;
                
                % increment wave counter
                waves = waves + 1;
            end
            
            if inWave && u <= threshold
                % falling edge
                inWave = false;
            end
        end
        output_data(o) = convertQuantToMag(waves, timeLength, strideLength);
        o = o + 1;
    end
    
    
    
end



function mag = convertQuantToMag(quantity, timeLength, strideLength)

    % having measured \quantity\ waves in \timeLength\ timesteps with \strideLength\ stride, what is
    % the underlying magnitude which was inputted?
    
    %  0 -> 0.0
    %  timeLength / strideLength -> 1.0
    %  0..1 -> mag = quant.strideLength / timeLength
    
    mag = quantity * strideLength / timeLength;


end


function [quantity, delay] = convertMagToQuant(mag, timeLength, strideLength)

    % given \mag\ input (0..1), how many waves should we produce in a given length of time?
    % 1.0 and 0.0 are the easiest to calculate: 
    %    0.0 -> 0 waves
    %    1.0 -> timeLength / strideLength waves
    %    for 0..1 : linear mapping between extreme values ^^
    %       -- mag.timeLength / strideLength
    
    freq = mag * timeLength / strideLength; % frequency of waves
    quantity = floor(freq);                 % number of waves to produce
    
    % how many time steps do I need to wait after triggering the wave before I can trigger the next
    % one?
    
    % 1.0 -> strideLength
    % 0.5 -> strideLength * 2
    % 0.0 -> timeLength
    
    % \quantity\ waves to produce in \timeLength\ steps
    % timeLength / freq = time to wait?
    
    % delay = stride / mag (for mag > 0)
    
    delay = timeLength;
    if mag > 0
        delay = strideLength / mag;
        if delay > timeLength
            delay = timeLength;
        end
    end
end






%% step function
function state = oreg_step(state, config)

    height     = config.height;
    width      = config.width;
    epsilon    = config.epsilon;
    grid_2     = config.grid_2;
    diff_coeff = config.diff_coeff;
    f          = config.f;
    q          = config.q;
    
    u_idx      = config.u_idx;
    v_idx      = config.v_idx;
    b_idx      = config.b_idx;
    p_idx      = config.p_idx;
    
    
    % store changes for t+1 in here
    del_uv     = zeros(height, width, 4);   %need 4 layers here, but we don't change the last 2
    
    for row = 1:height
        cprev    = state(row, 1, u_idx);   % x - 1
        cnext    = state(row, 2, u_idx);   % x + 1
        
        if row == 1, rprev_i = 1; else rprev_i = row - 1; end
%         rprev_i  = iif(row == 1, 1, row - 1);  % boundary condition
        if row == height, rnext_i = height; else rnext_i = row + 1; end
%         rnext_i  = iif(row == height, height, row + 1);  % boundary condition
        
        for col = 1:width
            this_u   = cnext;   % already calculated this value 
            this_v   = state(row, col, v_idx);
            if col == width, cnext_i = width; else cnext_i  = col + 1; end
%             cnext_i  = iif(col == width, width, col + 1);
            cnext    = state(row, cnext_i, u_idx);
            rprev    = state(rprev_i, col, u_idx);
            rnext    = state(rnext_i, col, u_idx);
            phi      = state(row, col, p_idx);       % illumination here
            boundary = state(row, col, b_idx);       % is this a hard boundary / vesicle?
            
            
            if boundary == 1
                del_uv(row, col, u_idx) = -this_u;   % kill the wave
                del_uv(row, col, v_idx) = -this_v;   % 
            else
                laplacian = (rprev + rnext + cprev + cnext - (4 * this_u));
                laplacian = laplacian / grid_2;
                
                del_uv(row, col, u_idx) = ( (...
                                        this_u - (this_u * this_u) - (f * this_v + phi) *  ...
                                        ((this_u - q) / (this_u + q)) ...
                                    ) / epsilon ...
                                    ) + (diff_coeff * laplacian);
                
                del_uv(row, col, v_idx) = this_u - this_v;
                
            end
            
            cprev = this_u;
            
        end
        
    end
    
    state = state + (del_uv * config.speed);

end

