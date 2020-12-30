function [X, rho, nar, ngibbs] = rmvnrnd(mu,sigma,N,A,b,rhoThr)

defaultRhoThr = 2.9e-4;

if nargin<6 || rhoThr<0,rhoThr = defaultRhoThr;end

if nargin<5
    rho = 1; nar = N; ngibbs = 0;
    X = mvnrnd(mu, sigma, N);
    return
end
    
A = A'; b = b';
mu = mu';

p = length(mu);                         % dimension
m = size(A,2);                          % number of constraints
if m==0
    A = zeros(p,1);
end

X = zeros(N,p);                        
nar = 0;                               
ngibbs = 0;
rho = 1; 

if rhoThr<1
    % Approach 1: accept-reject method
    n = 0; % no. accepted
    maxSample = 1e6;
    trials = 0; passes = 0;
    s = N;
    while n<N && ( rho>rhoThr || s<maxSample)
        R  = mvnrnd(mu,sigma,s);
        R = R(sum(R*A<=repmat(b,size(R,1),1),2)==size(A,2),:);
        if size(R,1)>0
            X((n+1):min(N,(n+size(R,1))),:) = R(1:min(N-n,size(R,1)),:);
            nar = nar + min(N,(n+size(R,1))) - n;
        end
        n = n + size(R,1); trials = trials + s;
        rho = n/trials;
        if rho>0
            s = min([maxSample,ceil((N-n)/rho),10*s]);
        else
            s = min([maxSample,10*s]);
        end
        passes = passes + 1;
    end
end

%
% Approach 2: Gibbs sampler of Christian Robert
%
if nar < N
    % choose starting point
    if nar > 0
        x = X(nar,:);
        discard = 0;
    elseif all(A'*mu' <= b')
        x = mu';
        discard = 1;
    else
        % need to find a feasible point
        xf = chebycenter(A',b',1);
        if ~all(A'*xf <= b') 
            error('Failed to find a feasible point')
        end
        % would like a feasible point near mu
        Atheta = -A'*(xf-mu');
        btheta = b' - A'*xf;
        Atheta = [Atheta; 1; -1];
        btheta = [btheta; 1; 0 ];
        ftheta = -1;
        % x = linprog(f,A,b,Aeq,beq,lb,ub,x0,options)
        options = optimset;
        options = optimset(options,'Display', 'off');
        theta = linprog(ftheta,Atheta,btheta,[],[],[],[],[],options);
        x = mu' + (1-theta)*(xf-mu');
        
        discard = 1;
    end
    % set up inverse Sigma
    SigmaInv = inv(sigma);
    n = nar;
    while n<N
        % choose p new components
        for i = 1:p
            Sigma_i_iInv = SigmaInv([1:(i-1) (i+1):p],[1:(i-1) (i+1):p]) - ...
                SigmaInv([1:(i-1) (i+1):p],i)*SigmaInv([1:(i-1) (i+1):p],i)' ...
                / SigmaInv(i,i);
            x_i = x([1:(i-1) (i+1):p]);
            mu_i = mu([1:(i-1) (i+1):p]);
            mui = mu(i) + Sigmai_i' * Sigma_i_iInv * (x_i' - mu_i');
            s2i = sigma(i,i) - Sigmai_i'*Sigma_i_iInv*Sigmai_i; x m matrix derived from A by removing
            % the i-th row.
            A_i = A([1:(i-1) (i+1):p],:);
            Ai = A(i,:);
            c = (b-x_i*A_i)./Ai;
            lb = max(c(Ai<0));
            if isempty(lb), lb=-Inf; end
            ub = min(c(Ai>0));
            if isempty(ub), ub=Inf; end
            x(i) = mui+TruncatedGaussian(-sqrt(s2i),[lb ub]-mui);
        end
        if discard <= 0
            n = n + 1;
            X(n,:) = x; 
            ngibbs = ngibbs+1;
        end
        discard = discard-1;
    end
end
