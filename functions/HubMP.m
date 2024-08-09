function [Ilocs, gamVals] = HubMP(A, Y, K, sig2, q, T)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HubMP: CL Matching Pursuit (MP) using Huber's loss. This greedy covariance 
% learning methid is described in the ICASSP-2025 paper.
%
% INPUT: 
%   A       - Dictionary of size L x N 
%   Y       - Matrix of L x M 
%   K       - number of non-zero sources
%   sig2    - Noise variance
%   q       - quantile that determines c^2 (Default is q=0.8)
%   T       - max # of FP algorithm
%
% OUTPUT:
%   Ilocs      - Support of non-zeros signal powers (K-vector)
%   gamVals  - Power estimates (N-vector)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = size(A, 1);                 % Number of pilots
N = size(A, 2);                 % Number of MTDs
M = size(Y, 2);                 % Number of antennas (=measurement vecs)      
Sigmainv = (1/sig2)*eye(L);    % Initial \Sigma^(0)

Ilocs = zeros(1,K);            % A sparse active device set
gamVals = zeros(1,K);        

ufun = @(t,c) ((t<=c) + (c./t).*(t>c)); % weight function u(t)
rhofun = @(t,c) ( t.*(t<=c) + c*(log(t/c)+1).*(t>c) );

if nargin < 5
    q = 0.8;
    T = 10;
end
csq = chi2inv(q,2*L)/2;  
b = chi2cdf(2*csq,2*(L+1))+(csq/L)*(1-q);

const = (1/(b*M));
onetoN = 1:N;
gam = zeros(1,N);        % Initialize power estimates

% Main Iteration Phase
for k = 1:K

    MD = real(sum(conj(Y).*(Sigmainv*Y))); % y_m'*Sigma^-1*y_m
    B = Sigmainv*A;    % Sigma^-1 a_n , n = 1,..,N
    D = subplus(real(sum(conj(A).*B))); % a_n'Sigma^-1 a_n
    tmp = setdiff(onetoN,Ilocs(1:k-1));
    err = zeros(1,N);

    for i = tmp

        a_i  = A(:,i);
        b_i  = B(:,i); 
        num_i = abs(sum(conj(Y).*b_i)).^2;

        % <-- FP iterations start
        gam0 = gam(i);
        % 10 is the max number of FP iterations!
        for t = 1:T

            c_i = gam0/(1+gam0*D(i));
            eta = MD - c_i*num_i;
            C = const*(Y.*repmat(ufun(eta,csq),L,1))*Y'; 
            gam1 = subplus(real(b_i'*C*b_i)/D(i)^2 -1/D(1));

            if gam1==gam0 || abs(gam1 - gam0)/gam1 < 10^-3
                break;
            end
            gam0 = gam1;
        end
        gam(i) = gam1;
        % <-- FP iteratations end

        %- Compute the error 
        c_i = gam1/(1+gam1*D(i));
        eta = MD -  c_i*num_i;
        invSig_i = Sigmainv -  c_i*b_i*(b_i');
        err(i) = mean(rhofun(eta,csq))/b - real(log(det(invSig_i)));
           
    end

    % Find the index with minimum error 
    [~, indx] = min(err(tmp));
    ix = tmp(indx);
    gamVals(k) = gam(ix);
    Ilocs(k) = ix;

    % Update the Sigma
    d = gamVals(k); 
    c = real(A(:,ix)'*B(:,ix));
    Sigmainv = Sigmainv - (d/(1+d*c))*B(:,ix)*(B(:,ix)');

end
