function [YY1_P,YY1_Q,YY1_R,YY1_S,YY2] = selection_elitism(P,Q,R,S,F,p)
% P - population, F - fitness value, p - population size
%elitism selection:
% Often to get better parameters, strategies with partial reproduction are used. 
% One of them is elitism, in which a small portion of the best individuals 
% from the last generation is carried over (without any changes) to the next one.


[sortScore,I] = sort(F,'descend');
Y1_P = P(:,:,I);
Y1_Q = Q(:,:,I);
Y1_R = R(:,:,I);
Y1_S = S(:,:,I);

Y1_P = Y1_P(:,:,1:p);
Y1_Q = Y1_Q(:,:,1:p);
Y1_R = Y1_R(:,:,1:p);
Y1_S = Y1_S(:,:,1:p);
sortScore = sortScore(1:p);


YY1_P = Y1_P; %store chromosomes
YY1_Q = Y1_Q; %store chromosomes
YY1_R = Y1_R; %store chromosomes
YY1_S = Y1_S; %store chromosomes
YY2 = sortScore; %store scores
end