function [Z_P,Z_Q,Z_R,Z_S]=mutation(P,Q,R,S,n,gene_mutation_rate,mut_vary)
% PQRS = population
% n = chromosomes to be mutated
% pages = number of chromosomes
[term_n,term_m,pages]=size(P);
Z=zeros(term_n,term_m,n); % term_n x term_m x number of mutation
for i = 1:n %for every mutated chromosome
    
    %select a chromosome to mutate
    r1=randi(pages); %generate random int from 1 to pages
    Pnm=P(:,:,r1); % random parent %call the random parents
    Qnm=Q(:,:,r1); % random parent %call the random parents
    Rnm=R(:,:,r1); % random parent %call the random parents
    Snm=S(:,:,r1); % random parent %call the random parents

    %select which genes (consts) to mutate
    r2_n=randi(term_n); %generate rand int from 1 to 100
    r2_m=randi(term_m); %generate rand int from 1 to 100

    % mutation
    %determine mutation rate
    mut_rate_n = rand(1)*gene_mutation_rate; %generate rand num [0,gene_mutation_rate]
    mut_rate_m = rand(1)*gene_mutation_rate; %generate rand num [0,gene_mutation_rate]
    
    const_mut_num_n = round(mut_rate_n*term_n);  %determine numbers of constants in term_n will be mutated
    const_mut_num_m = round(mut_rate_m*term_m);  %determine numbers of constants in term_m will be mutated

    %generate which constants indices will be mutated
    const_mut_n = randi(term_n,const_mut_num_n,1); %term_n x 1
    const_mut_m = randi(term_m,1,const_mut_num_m); %1 x term_m
%     const_mut_nm = [const_mut_n,const_mut_m']; %term_n x 2

    bin_mut = zeros(term_n,term_m);
    bin_mut(const_mut_n,const_mut_m') = 1; 

    % WHICH CONSTS INDICES ARE MUTATED ARE THE SAME FOR ALL PQRS
    %The consts varied are diff for each PQRS

    % P 
    %mutate by changing increase or decreasing value by mut_vary

    binmat_inc_dec = round(rand(term_n,term_m)); %of the consts to vary, which to increase or decrease
    inc = binmat_inc_dec.*(1+rand(term_n,term_m)*mut_vary); %inc = binmat _inc_dec*[1.01,1.1] => only 1 (increase) given 1.01-1.1 value
    dec = (1-binmat_inc_dec).*(1-rand(term_n,term_m)*mut_vary); %dec = (1- binmat_inc_dec)*[0.9,0.99] => make sure 1- binmat _inc_dec. only 0 (decrease which are now 1 values) given 0.9-0.99 value
    inc_dec = inc+dec; %Add inc and dec together to get a matrix of increase and decrease percentages 
    vary = inc_dec.*bin_mut; %Inc_dec*binmat_vary should assign percent variations to respective 1 spots
    final_vary = vary+(1-bin_mut); %change all 0s to 1 to maintain same const and not 0
    Z_P(:,:,i) = Pnm.*final_vary;

    % Q
    %mutate by changing increase or decreasing value by mut_vary
    binmat_inc_dec = round(rand(term_n,term_m)); %of the consts to vary, which to increase or decrease
    inc = binmat_inc_dec.*(1+rand(term_n,term_m)*mut_vary); %inc = binmat _inc_dec*[1.01,1.1] => only 1 (increase) given 1.01-1.1 value
    dec = (1-binmat_inc_dec).*(1-rand(term_n,term_m)*mut_vary); %dec = (1- binmat_inc_dec)*[0.9,0.99] => make sure 1- binmat _inc_dec. only 0 (decrease which are now 1 values) given 0.9-0.99 value
    inc_dec = inc+dec; %Add inc and dec together to get a matrix of increase and decrease percentages 
    vary = inc_dec.*bin_mut; %Inc_dec*binmat_vary should assign percent variations to respective 1 spots
    final_vary = vary+(1-bin_mut); %change all 0s to 1 to maintain same const and not 0
    Z_Q(:,:,i) = Qnm.*final_vary;

    % R
    %mutate by changing increase or decreasing value by mut_vary
    binmat_inc_dec = round(rand(term_n,term_m)); %of the consts to vary, which to increase or decrease
    inc = binmat_inc_dec.*(1+rand(term_n,term_m)*mut_vary); %inc = binmat _inc_dec*[1.01,1.1] => only 1 (increase) given 1.01-1.1 value
    dec = (1-binmat_inc_dec).*(1-rand(term_n,term_m)*mut_vary); %dec = (1- binmat_inc_dec)*[0.9,0.99] => make sure 1- binmat _inc_dec. only 0 (decrease which are now 1 values) given 0.9-0.99 value
    inc_dec = inc+dec; %Add inc and dec together to get a matrix of increase and decrease percentages 
    vary = inc_dec.*bin_mut; %Inc_dec*binmat_vary should assign percent variations to respective 1 spots
    final_vary = vary+(1-bin_mut); %change all 0s to 1 to maintain same const and not 0
    Z_R(:,:,i) = Rnm.*final_vary;

    % S
    %mutate by changing increase or decreasing value by mut_vary
    binmat_inc_dec = round(rand(term_n,term_m)); %of the consts to vary, which to increase or decrease
    inc = binmat_inc_dec.*(1+rand(term_n,term_m)*mut_vary); %inc = binmat _inc_dec*[1.01,1.1] => only 1 (increase) given 1.01-1.1 value
    dec = (1-binmat_inc_dec).*(1-rand(term_n,term_m)*mut_vary); %dec = (1- binmat_inc_dec)*[0.9,0.99] => make sure 1- binmat _inc_dec. only 0 (decrease which are now 1 values) given 0.9-0.99 value
    inc_dec = inc+dec; %Add inc and dec together to get a matrix of increase and decrease percentages 
    vary = inc_dec.*bin_mut; %Inc_dec*binmat_vary should assign percent variations to respective 1 spots
    final_vary = vary+(1-bin_mut); %change all 0s to 1 to maintain same const and not 0
    Z_S(:,:,i) = Rnm.*final_vary;

end