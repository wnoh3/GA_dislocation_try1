term_n = 10;
term_m = 10;
vary_const = 0.1;
binmat_vary = round(rand(term_n,term_m)); %which consts to vary or not
binmat_inc_dec = round(rand(term_n,term_m)); %of the consts to vary, which to increase or decrease
inc = binmat_inc_dec.*(1+rand(term_n,term_m)*vary_const); %inc = binmat _inc_dec*[1.01,1.1] => only 1 (increase) given 1.01-1.1 value
dec = (1-binmat_inc_dec).*(1-rand(term_n,term_m)*vary_const); %dec = (1- binmat_inc_dec)*[0.9,0.99] => make sure 1- binmat _inc_dec. only 0 (decrease which are now 1 values) given 0.9-0.99 value
inc_dec = inc+dec; %Add inc and dec together to get a matrix of increase and decrease percentages
vary = inc_dec.*binmat_vary %Inc_dec*binmat_vary should assign percent variations to respective 1 spots
final_vary = vary+(1-binmat_vary) %change all 0s to 1 to maintain same const and not 0