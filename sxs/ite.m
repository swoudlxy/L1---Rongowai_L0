% This function iteratively solve the positions of specular points
% based on the WGS84 model
% Inputs:
% 1) Tx and Rx coordiates in ECEF XYZ
% Ouputputs
% 1) Sx coordinate in ECEF XYZ

function sx = ite(tx,rx)

%tx = [-22052644.76142904 -13074480.840453148 7530822.467494242];
%rx = [-4913538.73980396	280776.609248639	-4052785.53858799];

% Fibonacci sequence
k = 0:1:59;
term1 = ((1+sqrt(5))/2).^(k+1);
term2 = ((1-sqrt(5))/2).^(k+1);
term3 = term1-term2;
F = term3/sqrt(5);      %Fibonacci sequence

s_t2r = norm(rx-tx);

%determine iteration
N = find(F>s_t2r,1,'first');

%first iteration parameters
a = rx;
b = tx;

for k = 1:N-2

    %testing points
    term1 = F(N-k-1)/F(N-k+1);
    term2 = F(N-k)/F(N-k+1);
        
    m_lambda = a+term1*(b-a);
    m_mu = a+term2*(b-a);
    
    %nadir points
    s_lambda = nadir(m_lambda);
    s_mu = nadir(m_mu);
    
    %propagation distance
    f_lambda = pdis(tx,rx,s_lambda);
    f_mu = pdis(tx,rx,s_mu);
    
    %parameters for next iteration
    if f_lambda > f_mu
       a = m_lambda;
       b = b;
    elseif f_lambda <= f_mu
        a = a;
        b = m_mu;
    end
end

sx = nadir(m_lambda);