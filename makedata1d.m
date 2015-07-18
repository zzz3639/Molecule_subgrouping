function [ dataset, mu ] = makedata1d( N, mu0, alpha0, sigma0, D )
%MAKEDATA1D Summary of this function goes here
%   Detailed explanation goes here
rng('shuffle');
Gnumber=length(mu0);
T=length(N);
dataset=cell(T,1);

%generate trajectories
mu=zeros(T,Gnumber);
for i=1:Gnumber
    mu(1,i)=mu0(i);
end
for t=2:T
    for i=1:Gnumber
        mu(t,i)=mu(t-1,i)+D*randn();
    end
end

%generate observations
for t=1:T
    n=N(t);
    datathis=zeros(n,1);
    nthis=mnrnd(n,alpha0);
    k=1;
    for i=1:Gnumber
        datathis(k:k+nthis(i)-1,1)=mu(t,i)+sigma0(i)*randn(nthis(i),1);
        k=k+nthis(i);
    end
    dataset{t}=datathis;
end

[H,X]=hist(dataset{1},20);

for t=1:T
    H=hist(dataset{t},X);
    plot(X,H,'*');
    pause(1);
end

end

