%% Computational Cost

K = 50;
S_max = 3*K;
T = 3;
r = 0.05;
sigma = 0.25;
M = 100:100:1000;
N = 1000:10000:100000;
dt = T./N;
ds = S_max./M;

time = zeros(length(M),4);

for i = 1:length(M)
    
    
    tic,[S,V,tao] = Explicit_B_S(K,S_max,r,T,sigma,ds(i),dt(i),'CALL');
    time(i,1) =toc;
    
    tic,[S,V,tao] = Semi_implicit_B_S(K,S_max,r,T,sigma,ds(i),dt(i),'CALL');
    time(i,2) =toc;
    
    tic,[S,V,tao] = Fully_Implicit_B_S(K,S_max,r,T,sigma,ds(i),dt(i),'CALL');
    time(i,3) =toc;
    
    tic,[S,V,tao] = CN_B_S(K,S_max,r,T,sigma,ds(i),dt(i),'CALL');
    time(i,4) =toc;
    
end
