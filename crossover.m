function [Z_P,Z_Q,Z_R,Z_S]=crossover(P,Q,R,S,n)
% P = population
% n = number of pairs of chromosomes to be crossovered
% pages = number of chromosomes
[term_n,term_m,pages]=size(P);
Z_P=zeros(term_n,term_m,2*n); %Z = term_n,term_m,2*n = 100x100xdouble the number of crossover
Z_Q=zeros(term_n,term_m,2*n); %Z = term_n,term_m,2*n = 100x100xdouble the number of crossover
Z_R=zeros(term_n,term_m,2*n); %Z = term_n,term_m,2*n = 100x100xdouble the number of crossover
Z_S=zeros(term_n,term_m,2*n); %Z = term_n,term_m,2*n = 100x100xdouble the number of crossover
for i = 1:n %for every crossover chromosome pair to be created
        
    %cut off point should be the same for all constants!!

    % Determine which chromosomes to switch
    r1=randi(pages,1,2); %generate ints (chromosomes) from 1 to x1 of array: 1 x 2
    while r1(1)==r1(2) %if r1(1) is equal to r1(2) = if same chromosome
        r1=randi(pages,1,2); %generate another pair
    end
    %cut off point
    r2_n=1+randi(term_n-1); % random cutting point %generate cutting point in chromosome
    r2_m=1+randi(term_m-1); % random cutting point %generate cutting point in chromosome

    %P 
    A1=P(:,:,r1(1)); % parent 1 %get parent 1
    A2=P(:,:,r1(2)); % parent 2 %get parent 2
    B1=A1(r2_n:term_n,r2_m:term_m); %copy parent 1 from cut off point to end
    A1(r2_n:term_n,r2_m:term_m)=A2(r2_n:term_n,r2_m:term_m); %copy parent 2 from cut off point to end to parent 1 from cut pt to end
    A2(r2_n:term_n,r2_m:term_m)=B1; %copy parent 1 from cut off pt to end to parent 2
    Z_P(:,:,2*i-1)=A1; % offspring 1 %save offspring 1
    Z_P(:,:,2*i)=A2; % offspring 2 %save offspring 2

    %Q 
    A1=Q(:,:,r1(1)); % parent 1 %get parent 1
    A2=Q(:,:,r1(2)); % parent 2 %get parent 2
    B1=A1(r2_n:term_n,r2_m:term_m); %copy parent 1 from cut off point to end
    A1(r2_n:term_n,r2_m:term_m)=A2(r2_n:term_n,r2_m:term_m); %copy parent 2 from cut off point to end to parent 1 from cut pt to end
    A2(r2_n:term_n,r2_m:term_m)=B1; %copy parent 1 from cut off pt to end to parent 2
    Z_Q(:,:,2*i-1)=A1; % offspring 1 %save offspring 1
    Z_Q(:,:,2*i)=A2; % offspring 2 %save offspring 2

    %R 
    A1=R(:,:,r1(1)); % parent 1 %get parent 1
    A2=R(:,:,r1(2)); % parent 2 %get parent 2
    B1=A1(r2_n:term_n,r2_m:term_m); %copy parent 1 from cut off point to end
    A1(r2_n:term_n,r2_m:term_m)=A2(r2_n:term_n,r2_m:term_m); %copy parent 2 from cut off point to end to parent 1 from cut pt to end
    A2(r2_n:term_n,r2_m:term_m)=B1; %copy parent 1 from cut off pt to end to parent 2
    Z_R(:,:,2*i-1)=A1; % offspring 1 %save offspring 1
    Z_R(:,:,2*i)=A2; % offspring 2 %save offspring 2

    %S 
    A1=S(:,:,r1(1)); % parent 1 %get parent 1
    A2=S(:,:,r1(2)); % parent 2 %get parent 2
    B1=A1(r2_n:term_n,r2_m:term_m); %copy parent 1 from cut off point to end
    A1(r2_n:term_n,r2_m:term_m)=A2(r2_n:term_n,r2_m:term_m); %copy parent 2 from cut off point to end to parent 1 from cut pt to end
    A2(r2_n:term_n,r2_m:term_m)=B1; %copy parent 1 from cut off pt to end to parent 2
    Z_S(:,:,2*i-1)=A1; % offspring 1 %save offspring 1
    Z_S(:,:,2*i)=A2; % offspring 2 %save offspring 2
end