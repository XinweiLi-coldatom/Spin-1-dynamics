clear
clc
global N dimvec H
%% definition
N = 20; % atom number
mvec = -N:N; % values of Lz
mn = length(mvec);
dim = floor((N-abs(mvec))/2)+1; % dimension of sub-Hilbert space
dimvec = cumsum(dim); 
dimvec = [0,dimvec];
hdim = dimvec(end); % dimension of full Hilbert space
%% operators
Lpc1 = [];
Lpc2 = [];
Lpc3 = [];
Lpr1 = [];
Lpr2 = [];
Lpr3 = [];
Lpv1 = [];
Lpv2 = [];
Lpv3 = [];
Lzv = zeros(hdim,1);
N1v = zeros(hdim,1);
N0v = zeros(hdim,1);
Nm1v = zeros(hdim,1);
cnt = 0;

for midx = 1:mn
    m = mvec(midx); % fixed m subspace   
    kvec = max(0,-m):floor((N-m)/2); % range of k for a certain m
  
    for kidx = 1:length(kvec)
        cnt = cnt + 1;
        k = kvec(kidx); % for a certain k        
        % The state is |k+m,N-2k-m,k>
        ip = basisip(m,k); % index of basis
        
        % Find the matrix element of a_1^\dag*a_0
        % <k+m+1,N-2k-m-1,k|a_1^\dag*a_0|k+m,N-2k-m,k> = sqrt(N-2k-m)*sqrt(k+m+1)
        if(m<N) % N-m-2k>0           
            rip = basisip(m+1,k);                
            vau = sqrt((N-2*k-m)*(k+m+1));
            Lpr1 = [Lpr1,rip]; % row index
            Lpc1 = [Lpc1,ip]; % column index
            Lpv1 = [Lpv1,vau]; % element of operator matrix    
        end
        
        % Find the matrix element of a_0^\dag*a_{-1}
        % <k+m,N-2k-m+1,k-1|a_0^\dag*a_{-1}|k+m,N-2k-m,k> = sqrt(k)*sqrt(N-2k-m+1);
        if(k>0)
            rip = basisip(m+1,k-1);                    
            vau = sqrt(k*(N-2*k-m+1));
            Lpr2 = [Lpr2,rip];% row index
            Lpc2 = [Lpc2,ip]; % column index
            Lpv2 = [Lpv2,vau];% element of operator matrix
        end
        
        % Find the matrix element of a_1^\dag*a_{-1}
        % <k+m+1,N-2k-m,k-1|a_1^\dag*a_{-1}|k+m,N-2k-m,k> = sqrt(k)*sqrt(k+m+1);
        if(k>0)
            rip = basisip(m+2,k-1);                    
            vau = sqrt(k*(k+m+1));
            Lpr3 = [Lpr3,rip];% row index
            Lpc3 = [Lpc3,ip]; % column index
            Lpv3 = [Lpv3,vau];% element of operator matrix
        end

        Lzv(cnt) = m;
        N0v(cnt) = N-2*k-m;
        N1v(cnt) = k+m;
        Nm1v(cnt) = k;
        
    end
end
% operators
Lp1 = sparse(Lpr1,Lpc1,Lpv1,hdim,hdim); % a_1^\dag*a_0
Lp2 = sparse(Lpr2,Lpc2,Lpv2,hdim,hdim); % a_0^\dag*a_{-1}
Lp3 = sparse(Lpr3,Lpc3,Lpv3,hdim,hdim); % a_1^\dag*a_{-1}
Lm1 = Lp1';
Lm2 = Lp2';
Lm3 = Lp3';

Lp = Lp1+Lp2; % L_{+}
Lm = Lp'; % L_{-}
Lz = sparse(1:hdim,1:hdim,Lzv,hdim,hdim);
N0 = sparse(1:hdim,1:hdim,N0v,hdim,hdim);
N1 = sparse(1:hdim,1:hdim,N1v,hdim,hdim);
Nm1 = sparse(1:hdim,1:hdim,Nm1v,hdim,hdim);
Lx = (Lp + Lm)/sqrt(2);
Ly = (Lp - Lm)/sqrt(2)/1i;
L2 = Lx^2+Ly^2+Lz^2;
%% initial state |0,N,0> -> m=0, k=0
m = 0;
k = 0;
id = basisip(m,k);
psi = zeros(hdim,1);
psi(id) = 1; 
%% construction of Hamiltonian
c = -1; % spin-mixing rate
q = abs(c);
lambda = c/2/N;
H = lambda*L2-q*N0; % spin-mixing
%% spin-mixing evolution
% Evolution time
ti = 0;
tf = 20;
nt = 101;
dt = (tf-ti)/(nt-1);
tvec = linspace(0,tf,nt);

% Observable(mean atomic number)
n0mean = zeros(1,nt);
n1mean = zeros(1,nt);
n0sqmean = zeros(1,nt);

% evolution
t1 = ti;
for k = 1:nt 

    n0mean(k) = psi'*N0*psi/N;
    n1mean(k) = psi'*N1*psi/N;
    n0sqmean(k) = psi'*N0^2*psi/N/N;
    
%     if mod(k,10) == 0
%         disp(k)
%         disp(psi'*(N0+2*Nm1)*psi/N)
%     end
    
    t2 = k*dt;
    [t,yt] = ode45(@fff,[t1 t2],psi');
    psi = yt(end,:)';
    t1 = t2;
end

figure
plot(tvec,n0mean,'k','linewidth',2)
