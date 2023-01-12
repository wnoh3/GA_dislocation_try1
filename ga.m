clear all
close all
clc
%Admin Info ---------------------------------------------------------------
%{
-Derivative code from:
https://learnwithpanda.com/2020/07/22/genetic-algorithm-general-concept-matlab-code-and-example/
-Code adapted for Genetic Algorithm optimize constants of isolated
dislocation
-Author: William Noh
-UIUC 
-Date: January 6, 2022
%}
%Additional Information----------------------------------------------------
%Chromosome: Constants for Single Isolated Dislocation
%--------------------------------------------------------------------------
p=10; % Population size: 100 Chromosomes with (n x m) matrix
c=3; % number of pairs of chromosomes to be crossovered %30 percent crossover
m=3; % number chromosomes to be mutated %30 mutation
tg=2; % Total number of generations 
vary_const = 0.1; %how much the single isolated dislocation consts will vary
gene_mutation_rate = 0.3; %how much of the chromosome should be mutated
mut_vary = 0.1; %up to much a mutation will vary const whether increase or decrease
%--------------------------------------------------------------------------
%set following from const_fp down 
const_fp = 'D:\Research\FinalAFOSRReview\SpectralAnalysis\AnalyticalTry2\Dec6Lpimplementation\spectral\6Jan_singledisloconst\n100_m100_sub288x288_LxLy345.6A_PBC_i0_optnonpiece_nonpieceplot_V2\';
subsetx = 288; %notneeded
subsety = 288; %notneeded

%set this to the number of dislocations and spacing
d = 130.2903-67.013; %A %35_1_0

%--------------------------------------------------------------------------
%MD info ------------------------------------------------------------------
MD_path = 'D:\Research\FinalAFOSRReview\GeneticAlgorithm\FirstTry\MD_GBs\Jan8_V2\GB35_1_0\';
GB = 'GB35_1_0';
MD_path = strcat(MD_path,GB);
%Destination Directory-----------------------------------------------------
dest = 'D:\Research\FinalAFOSRReview\GeneticAlgorithm\FirstTry\dislocation_try1\';
fpname = 'Jan9_test';
destdir = strcat(dest,fpname,'\');
mkdir(destdir)

%--------------------------------------------------------------------------
%fourier stuff
p = ceil(345.6/d); %# of d needed to exceed single dislo influence
Lx = p*d;%set as function of d AND > 345.6A (dislo influence)
Ly = Lx;
PBC_i = 2*p-1; %recalculate this if you need to prove formula is true

%import constants from spectral method of analytical single dislocation
Pnm_S11 = readmatrix(strcat(const_fp,'Pnm_S11.csv'));
Pnm_S12 = readmatrix(strcat(const_fp,'Pnm_S12.csv'));
Pnm_S22 = readmatrix(strcat(const_fp,'Pnm_S22.csv'));

Qnm_S11 = readmatrix(strcat(const_fp,'Qnm_S11.csv'));
Qnm_S12 = readmatrix(strcat(const_fp,'Qnm_S12.csv'));
Qnm_S22 = readmatrix(strcat(const_fp,'Qnm_S22.csv'));

Rnm_S11 = readmatrix(strcat(const_fp,'Rnm_S11.csv'));
Rnm_S12 = readmatrix(strcat(const_fp,'Rnm_S12.csv'));
Rnm_S22 = readmatrix(strcat(const_fp,'Rnm_S22.csv'));

Snm_S11 = readmatrix(strcat(const_fp,'Snm_S11.csv'));
Snm_S12 = readmatrix(strcat(const_fp,'Snm_S12.csv'));
Snm_S22 = readmatrix(strcat(const_fp,'Snm_S22.csv'));

fig_11 = figure(1);
title('Blue - Average         Red - Maximum S11');
xlabel('Generation')
ylabel('Objective Function Value')
hold on

fig_12 = figure(2);
title('Blue - Average         Red - Maximum S12');
xlabel('Generation')
ylabel('Objective Function Value')
hold on

fig_22 = figure(3);
title('Blue - Average         Red - Maximum S22');
xlabel('Generation')
ylabel('Objective Function Value')
hold on

%generate initial population by randomly varying constants by 10%
[P_11,Q_11,R_11,S_11]=population(p,Pnm_S11,Qnm_S11,Rnm_S11,Snm_S11,vary_const); % (n x m) x p = row page sheet
[P_12,Q_12,R_12,S_12]=population(p,Pnm_S12,Qnm_S12,Rnm_S12,Snm_S12,vary_const); % (n x m) x p = row page sheet
[P_22,Q_22,R_22,S_22]=population(p,Pnm_S22,Qnm_S22,Rnm_S22,Snm_S22,vary_const); % (n x m) x p = row page sheet

K_11=0;
K_12=0;
K_22=0;
[term_n,term_m,n]=size(P_22);
P1 = 0;
for i=1:tg   
%crossover-----------------------------------------------------------------
    [Cr_P_11,Cr_Q_11,Cr_R_11,Cr_S_11]=crossover(P_11,Q_11,R_11,S_11,c); %create crossovers
    [Cr_P_12,Cr_Q_12,Cr_R_12,Cr_S_12]=crossover(P_12,Q_12,R_12,S_12,c); %create crossovers
    [Cr_P_22,Cr_Q_22,Cr_R_22,Cr_S_22]=crossover(P_22,Q_22,R_22,S_22,c); %create crossovers
    fprintf("Crossover for gen %d finished.\n",i);

%mutation------------------------------------------------------------------
    [Mu_P_11,Mu_Q_11,Mu_R_11,Mu_S_11]=mutation(P_11,Q_11,R_11,S_11,m,gene_mutation_rate,mut_vary); %create mutations
    [Mu_P_12,Mu_Q_12,Mu_R_12,Mu_S_12]=mutation(P_12,Q_12,R_12,S_12,m,gene_mutation_rate,mut_vary); %create mutations
    [Mu_P_22,Mu_Q_22,Mu_R_22,Mu_S_22]=mutation(P_22,Q_22,R_22,S_22,m,gene_mutation_rate,mut_vary); %create mutations
    fprintf("Mutation for gen %d finished.\n",i);

%assign crossover and mutations to PQRS------------------------------------
    %11
    P_11(:,:,p+1:p+2*c)      =Cr_P_11; %append population array with crossovers
    P_11(:,:,p+2*c+1:p+2*c+m)=Mu_P_11; %append population array with mutations
    Q_11(:,:,p+1:p+2*c)      =Cr_Q_11; %append population array with crossovers
    Q_11(:,:,p+2*c+1:p+2*c+m)=Mu_Q_11; %append population array with mutations
    R_11(:,:,p+1:p+2*c)      =Cr_R_11; %append population array with crossovers
    R_11(:,:,p+2*c+1:p+2*c+m)=Mu_R_11; %append population array with mutations
    S_11(:,:,p+1:p+2*c)      =Cr_S_11; %append population array with crossovers
    S_11(:,:,p+2*c+1:p+2*c+m)=Mu_S_11; %append population array with mutations
	
    %12
    P_12(:,:,p+1:p+2*c)      =Cr_P_12; %append population array with crossovers
    P_12(:,:,p+2*c+1:p+2*c+m)=Mu_P_12; %append population array with mutations
    Q_12(:,:,p+1:p+2*c)      =Cr_Q_12; %append population array with crossovers
    Q_12(:,:,p+2*c+1:p+2*c+m)=Mu_Q_12; %append population array with mutations
    R_12(:,:,p+1:p+2*c)      =Cr_R_12; %append population array with crossovers
    R_12(:,:,p+2*c+1:p+2*c+m)=Mu_R_12; %append population array with mutations
    S_12(:,:,p+1:p+2*c)      =Cr_S_12; %append population array with crossovers
    S_12(:,:,p+2*c+1:p+2*c+m)=Mu_S_12; %append population array with mutations

    %22
    P_22(:,:,p+1:p+2*c)      =Cr_P_22; %append population array with crossovers
    P_22(:,:,p+2*c+1:p+2*c+m)=Mu_P_22; %append population array with mutations
    Q_22(:,:,p+1:p+2*c)      =Cr_Q_22; %append population array with crossovers
    Q_22(:,:,p+2*c+1:p+2*c+m)=Mu_Q_22; %append population array with mutations
    R_22(:,:,p+1:p+2*c)      =Cr_R_22; %append population array with crossovers
    R_22(:,:,p+2*c+1:p+2*c+m)=Mu_R_22; %append population array with mutations
    S_22(:,:,p+1:p+2*c)      =Cr_S_22; %append population array with crossovers
    S_22(:,:,p+2*c+1:p+2*c+m)=Mu_S_22; %append population array with mutations
    fprintf("Appending for gen %d finished.\n",i);

%Evaluation----------------------------------------------------------------
    E=evaluation(P_11,Q_11,R_11,S_11,P_12,Q_12,R_12,S_12,P_22,Q_22,R_22,S_22...
                    ,Lx,Ly,PBC_i,d,...
                    MD_path); %evaluate performance of each chromosome 3 x chromosomes array
    fprintf("Evaluation for gen %d finished.\n",i);

%Selection-----------------------------------------------------------------
    %11
    [P_11,Q_11,R_11,S_11,scores_11]=selection(P_11,Q_11,R_11,S_11,E(1,:),p); %returns 100 best selected chromosomes and fit scores
    %11
    [P_12,Q_12,R_12,S_12,scores_12]=selection(P_12,Q_12,R_12,S_12,E(2,:),p); %returns 100 best selected chromosomes and fit scores
    %11
    [P_22,Q_22,R_22,S_22,scores_22]=selection(P_22,Q_22,R_22,S_22,E(3,:),p); %returns 100 best selected chromosomes and fit scores
    fprintf("Selection for gen %d finished.\n",i);

%scores--------------------------------------------------------------------
    %11
    K_11(i,1)=sum(scores_11)/p; %find avg of fitness scores
    K_11(i,2)=scores_11(1); % %save best fit score
    %12
    K_12(i,1)=sum(scores_12)/p; %find avg of fitness scores
    K_12(i,2)=scores_12(1); % %save best fit score
    %22
    K_22(i,1)=sum(scores_22)/p; %find avg of fitness scores
    K_22(i,2)=scores_22(1); % %save best fit score
%plot scores---------------------------------------------------------------
    figure(1);
    plot(K_11(:,1),'b.'); drawnow
    hold on
    plot(K_11(:,2),'r.'); drawnow

    figure(2);
    plot(K_12(:,1),'b.'); drawnow
    hold on
    plot(K_12(:,2),'r.'); drawnow

    figure(3);
    plot(K_12(:,1),'b.'); drawnow
    hold on
    plot(K_12(:,2),'r.'); drawnow


end
fprintf("GA Finished.\n");

Max_fitness_value_S11=max(K_11(:,2)) %find best of best
Max_fitness_value_S12=max(K_12(:,2)) %find best of best
Max_fitness_value_S22=max(K_22(:,2)) %find best of best

Optimal_solution_P_11 = P_11(:,:,1); % Best chromosome
Optimal_solution_Q_11 = Q_11(:,:,1); % Best chromosome
Optimal_solution_R_11 = R_11(:,:,1); % Best chromosome
Optimal_solution_S_11 = S_11(:,:,1); % Best chromosome

Optimal_solution_P_12 = P_12(:,:,1); % Best chromosome
Optimal_solution_Q_12 = Q_12(:,:,1); % Best chromosome
Optimal_solution_R_12 = R_12(:,:,1); % Best chromosome
Optimal_solution_S_12 = S_12(:,:,1); % Best chromosome

Optimal_solution_P_22 = P_22(:,:,1); % Best chromosome
Optimal_solution_Q_22 = Q_22(:,:,1); % Best chromosome
Optimal_solution_R_22 = R_22(:,:,1); % Best chromosome
Optimal_solution_S_22 = S_22(:,:,1); % Best chromosome

%plot using optimal solution
plot_fourier(Optimal_solution_P_11,Optimal_solution_Q_11,Optimal_solution_R_11,Optimal_solution_S_11,...
            Optimal_solution_P_12,Optimal_solution_Q_12,Optimal_solution_R_12,Optimal_solution_S_12,...
            Optimal_solution_P_22,Optimal_solution_Q_22,Optimal_solution_R_22,Optimal_solution_S_22,...
            destdir)



