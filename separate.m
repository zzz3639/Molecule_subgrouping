function [ newX,newalpha,likelihood ] = separate( dataset,sigma )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
T=size(dataset,1);
%initialize
X=repmat([-1.5,1.5],T,1);
alpha=rand(T,1);
alpha=[alpha,1-alpha];
oldlikelihood=-realmax;
tc=0;
while (tc<100)
    likelihood=0;
    tc=tc+1;
%E step
    P=cell(T,1);
    for i=1:T
        P{i}=zeros(size(dataset{i},1),2);
        P{i}(:,1)=alpha(i,1)/sqrt(2*pi)/sigma*exp(-(dataset{i}-X(i,1)).^2/2*sigma^2);
        P{i}(:,2)=alpha(i,2)/sqrt(2*pi)/sigma*exp(-(dataset{i}-X(i,2)).^2/2*sigma^2);
        sump=sum(P{i},2);
        likelihood=likelihood+sum(log(sump));
        P{i}=P{i}./repmat(sump,1,2);
    end
%M step
    A=zeros(T,T,2);
    B=zeros(T,2);
    newalpha=zeros(T,2);
    
    newalpha(1,:)=sum(P{1});
    newalpha(T,:)=sum(P{T});
    
    A(1,1,1)=newalpha(1,1);
    
    A(T,T,1)=newalpha(T,1);
    
    A(1,1,2)=newalpha(1,2);
    
    A(T,T,2)=newalpha(T,2);
    
    B(1,:)=sum(P{1}.*repmat(dataset{1},1,2));
    B(T,:)=sum(P{T}.*repmat(dataset{T},1,2));
    
    for i=2:T-1
        newalpha(i,:)=sum(P{i});
        A(i,i,1)=newalpha(i,1);
        A(i,i,2)=newalpha(i,2);
        B(i,:)=sum(P{i}.*repmat(dataset{i},1,2));
    end        
    newX(:,1)=A(:,:,1)^(-1)*B(:,1);
    newX(:,2)=A(:,:,2)^(-1)*B(:,2);
    if likelihood-oldlikelihood<0.00001
        break
    else
        oldlikelihood=likelihood;
        X=newX;
        alpha=newalpha;
    end
end
end

