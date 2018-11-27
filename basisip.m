function ip = basisip(m,k)
    global N dimvec 
%   |k+m,N-2k-m,k>
    ip = dimvec(N + m + 1) + min(k,k+m) +1;
end