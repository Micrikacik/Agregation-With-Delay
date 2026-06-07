function xInitHist = genInitHist(x,dt,stepDelay,boundConds,dims)

[N, d] = size(x);

xInitHist = zeros([N, d, stepDelay]);

xInitHist(:,:,1) = x - sqrt(dt) * randn(N, d);    

% Apply BCs
switch boundConds
    case "NoBoundary"
        % No BCs
    case "Periodic"
        % Periodic BCs
        xInitHist(:,:,1) = mod(xInitHist(:,:,1),dims);
    case "Reflective"
        % Reflective BCs
        xInitHist(:,:,1) = abs(xInitHist(:,:,1));
        xInitHist(:,:,1) = dims - abs(dims - xInitHist(:,:,1));
    otherwise
        error('Invalid boundary conditions.');
end

for i = 2:stepDelay
    xInitHist(:,:,i) = xInitHist(:,:,i-1) - sqrt(dt) * randn(N, d);   

    % Apply BCs
    switch boundConds
        case "NoBoundary"
            % No BCs
        case "Periodic"
            % Periodic BCs
            xInitHist(:,:,i) = mod(xInitHist(:,:,i),dims);
        case "Reflective"
            % Reflective BCs
            xInitHist(:,:,i) = abs(xInitHist(:,:,i));
            xInitHist(:,:,i) = dims - abs(dims - xInitHist(:,:,i));
        otherwise
            error('Invalid boundary conditions.');
    end
end
