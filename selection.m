function [YY1_P,YY1_Q,YY1_R,YY1_S,YY2] = selection(P,Q,R,S,F,p)
% P - population, F - fitness value, p - population size
[term_n,term_m,pages]=size(P);
Y1_P = zeros(term_n,term_m,p); %chrom dim 1 x chrom dim 2 x p
Y1_Q = zeros(term_n,term_m,p); %chrom dim 1 x chrom dim 2 x p
Y1_R = zeros(term_n,term_m,p); %chrom dim 1 x chrom dim 2 x p
Y1_S = zeros(term_n,term_m,p); %chrom dim 1 x chrom dim 2 x p

% elite selection
e=3;
%Find the top three best chromosomes based on fitness score
for i = 1:e
    [r1 c1]=find(F==max(F)); %find which the row (r1) col (c1) ind of max fit chromosome r1 always 1 and c1 is actl

    Y1_P(:,:,i)=P(:,:,max(c1)); %store max indice of chromosome w highest fitness score
    Y1_Q(:,:,i)=Q(:,:,max(c1)); %store max indice of chromosome w highest fitness score
    Y1_R(:,:,i)=R(:,:,max(c1)); %store max indice of chromosome w highest fitness score
    Y1_S(:,:,i)=S(:,:,max(c1)); %store max indice of chromosome w highest fitness score

    P(:,:,max(c1))=[]; %empty chromosome of previously mentioned chromosome
    Q(:,:,max(c1))=[]; %empty chromosome of previously mentioned chromosome
    R(:,:,max(c1))=[]; %empty chromosome of previously mentioned chromosome
    S(:,:,max(c1))=[]; %empty chromosome of previously mentioned chromosome

    Fn(i)=F(max(c1)); %save score function of previously mentioned chromosome
    F(:,max(c1))=[]; %empty score func of previously mentioned chromosome
end
% https://en.wikipedia.org/wiki/Selection_(genetic_algorithm)
D=F/sum(F); % Determine selection probability %this is like a normalization
E= cumsum(D); % Determine cumulative probability  %cumsum is very weird function
N=rand(1); % Generate a vector constaining normalised random numbers %generate single random number between 0 and 1
d1=1; 
d2=e;
%assign chromosomes 4:100
while d2 <= p-e %d2 <= 97
    %find first accum norm score that is greater than N
    if N <= E(d1) % if rand num is <= first cumsum probs. E is accumulated normalized value

       Y1_P(:,:,d2+1)=P(:,:,d1); %store 4:100 = 1:97
       Y1_Q(:,:,d2+1)=Q(:,:,d1); %store 4:100 = 1:97
       Y1_R(:,:,d2+1)=R(:,:,d1); %store 4:100 = 1:97
       Y1_S(:,:,d2+1)=S(:,:,d1); %store 4:100 = 1:97

       Fn(d2+1)=F(d1); %store fitness functions
       N=rand(1);
       d2 = d2 +1; % 3 ->4 ...->98
       d1=1;
     else
       d1 = d1 + 1; %increment d1
    end
end
YY1_P = Y1_P; %store chromosomes
YY1_Q = Y1_Q; %store chromosomes
YY1_R = Y1_R; %store chromosomes
YY1_S = Y1_S; %store chromosomes
YY2 = Fn; %store scores
end