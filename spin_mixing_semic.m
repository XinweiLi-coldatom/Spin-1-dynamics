global ntra q c flag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
natom = 10000; % total atom number
ntra = 3000; % number of trajectories
c = -1;
q = abs(c);
ti = 0;
tf = 15; % evolution time
nt = 101;
dt = (tf-ti)/(nt-1);
tvec = linspace(ti,tf,nt);
%% initialization
delta = 1/sqrt(natom);
Lx = random('Normal',0,delta,ntra,1);
Nyz = random('Normal',0,delta,ntra,1);
Ly = random('Normal',0,delta,ntra,1);
Nxz = random('Normal',0,delta,ntra,1);

% Eq(3.29)
chip = atan(-(Ly+Nyz)./(Lx+Nxz));
index = find(Lx+Nxz<0);
chip(index) = atan(-(Ly(index)+Nyz(index))./(Lx(index)+Nxz(index)))+pi;

chim = atan((Ly-Nyz)./(Lx-Nxz));
index = find(Lx-Nxz<0);
chim(index) = atan((Ly(index)-Nyz(index))./(Lx(index)-Nxz(index)))+pi;

rho0 = 1/2 + sqrt(1/4-1/8*(((Lx+Nxz)./cos(chip)).^2+((Lx-Nxz)./cos(chim)).^2));
m = 1/8*(((Lx+Nxz)./cos(chip)).^2-((Lx-Nxz)./cos(chim)).^2)./rho0;
phs = chip + chim;

% initial polar state
yint(1,:) = sqrt((1-rho0+m)/2).*exp(1i*chip);
yint(2,:) = sqrt(rho0);
yint(3,:) = sqrt((1-rho0-m)/2).*exp(1i*chim);
psi = yint;
%% evolution

n0mean = zeros(1,nt);
n0sq = zeros(1,nt);

t1 = ti;
for k = 1:nt
   
    n0mean(k) = mean(abs(psi(2,:)).^2);
    n0sq(k) = mean(abs(psi(2,:)).^4);
    t2 = t1+dt;
    flag = 1;
    [t, yt] = ode45(@semic,[t1 t2],psi);
    psi = reshape(yt(end,:),3,ntra);    
    t1 = t2;
end
%%
figure
hold on
plot(tvec,n0mean,'b','LineWidth',2)