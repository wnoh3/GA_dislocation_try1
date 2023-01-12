%LD Lx is subject to change which means you must adjust spec_plot_xlength
%and ylength
%the 2.4 in subsetx and subsety could also be subject to change
function H=evaluation_parfor(Pnm_2_11_full,Qnm_2_11_full,Rnm_2_11_full,Snm_2_11_full,...
                        Pnm_2_12_full,Qnm_2_12_full,Rnm_2_12_full,Snm_2_12_full,...
                        Pnm_2_22_full,Qnm_2_22_full,Rnm_2_22_full,Snm_2_22_full,...
                        Lx,Ly,PBC_i,d,...
                        MDpath)

%import weighted area MD 
%MD Lx = 50.4;
MD_S11 = readmatrix(strcat(MDpath,'_S11'));
MD_S12 = readmatrix(strcat(MDpath,'_S12'));
MD_S22 = readmatrix(strcat(MDpath,'_S22'));

% %for test only
% [testx,testy] = size(MD_S11);
% MD_S11 = MD_S11(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);
% MD_S12 = MD_S12(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);
% MD_S22 = MD_S22(testx/2-21+1:testx/2+21,testx/2-21+1:testx/2+21);


N = length(MD_S11(:)); 

% PQRS = population
% n = chromosomes to be mutated
% pages = number of chromosomes
[term_n,term_m,pages]=size(Pnm_2_11_full);
H=zeros(3,pages); % S11 S12 S22 x number of chromosomes

%create n x m matrix
m_mat = zeros(term_m,term_n); %y=row x=col
n_mat = zeros(term_m,term_n); %y=row x=col
m_row = 1:term_m;
n_row = 1:term_n;
for len = 1:term_n %for every column
    m_mat(:,len) = m_row';
end %len = 1:n
for len = 1:term_m %for every row
    n_mat(len,:) = n_row;
end %len = 1:m


%Lx and Ly of reconstructed GB is set manually

%the following is defined above--------------------------------------------
%define new spectral plot values-------------------------------------------
spec_plot_xlength = 50.4;
spec_plot_ylength = 50.4;
spec_plot_subsetx = 2*spec_plot_xlength/2.4;
spec_plot_subsety = 2*spec_plot_xlength/2.4;
subx = 2*spec_plot_xlength/spec_plot_subsetx; % xdist value of each subsetx
suby = 2*spec_plot_ylength/spec_plot_subsety; % ydist value of each subsety
subArea = subx*suby;

topleftxycoord = [-spec_plot_xlength ,spec_plot_ylength];

%reconstruct stress field of GB using constants and fourier series
for i = 1:pages

    %assign pages
    Pnm_2_11 = Pnm_2_11_full(:,:,i);
    Qnm_2_11 = Qnm_2_11_full(:,:,i);
    Rnm_2_11 = Rnm_2_11_full(:,:,i);
    Snm_2_11 = Snm_2_11_full(:,:,i);

    Pnm_2_12 = Pnm_2_12_full(:,:,i);
    Qnm_2_12 = Qnm_2_12_full(:,:,i);
    Rnm_2_12 = Rnm_2_12_full(:,:,i);
    Snm_2_12 = Snm_2_12_full(:,:,i);

    Pnm_2_22 = Pnm_2_22_full(:,:,i);
    Qnm_2_22 = Qnm_2_22_full(:,:,i);
    Rnm_2_22 = Rnm_2_22_full(:,:,i);
    Snm_2_22 = Snm_2_22_full(:,:,i);
    
    
    %% plot 
    
    sigxx_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    sigxy_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    sigyy_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    
    %fprintf('start parfor\n');
    % %traverse through subset Voronoi Spectral
%     tic
    parfor j = 1:spec_plot_subsetx %4workers
        %fprintf("blah1\n")
        for k = 1:spec_plot_subsety
            
            % get sigma (vir) value at point
            % NOTE: For now we assume, subset is small enough to be
            % represented as a point and not weighted area average
    
            %find midpoint of subset
            %NOTE: coord sys is now dislo which is center of subset
            midpoint = topleftxycoord + [j*subx-0.5*subx,-k*suby+0.5*suby]; %formula verified Research Log12/20 pg 249
            xvalue = midpoint(1);
            yvalue = midpoint(2);
    
            %nonpiecewise--------------------------------------------------------------
            %
            dx = xvalue;
            dy = yvalue;
            
            sigxx_2(k,j) = sigxx_2(k,j) +  sum(Pnm_2_11.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
					                       +  Qnm_2_11.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
					                       +  Rnm_2_11.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
					                       +  Snm_2_11.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
					                       ,'all');
            sigxy_2(k,j) = sigxy_2(k,j) +  sum(Pnm_2_12.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
					                       +  Qnm_2_12.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
					                       +  Rnm_2_12.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
					                       +  Snm_2_12.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
					                       ,'all');
            sigyy_2(k,j) = sigyy_2(k,j) + sum(Pnm_2_22.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                   +  Qnm_2_22.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                   +  Rnm_2_22.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                   +  Snm_2_22.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                   ,'all');
            
            %find summation terms
            for add_i = 1:PBC_i
    
			    dxR = xvalue-( add_i*d);
                dxL = xvalue-(-add_i*d);
			    dy = yvalue;
    
                %dislocations contributions for left dislocations
                dx = xvalue+add_i*d;
                dy = yvalue;
                sigxx_2(k,j) = sigxx_2(k,j) + sum(Pnm_2_11.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Qnm_2_11.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           +  Rnm_2_11.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Snm_2_11.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           ,'all');
                sigxy_2(k,j) = sigxy_2(k,j) + sum(Pnm_2_12.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Qnm_2_12.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           +  Rnm_2_12.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Snm_2_12.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           ,'all');
                sigyy_2(k,j) = sigyy_2(k,j) + sum(Pnm_2_22.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Qnm_2_22.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           +  Rnm_2_22.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						           +  Snm_2_22.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						           ,'all');
                
    
    
                %dislocations contributions for right dislocations
                dx = xvalue-add_i*d;
                dy = yvalue;
    
		        sigxx_2(k,j) = sigxx_2(k,j) + sum(Pnm_2_11.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Qnm_2_11.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       +  Rnm_2_11.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Snm_2_11.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       ,'all');
		        sigxy_2(k,j) = sigxy_2(k,j) + sum(Pnm_2_12.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Qnm_2_12.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       +  Rnm_2_12.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Snm_2_12.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       ,'all');
		        sigyy_2(k,j) = sigyy_2(k,j) + sum(Pnm_2_22.*sin(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Qnm_2_22.*sin(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       +  Rnm_2_22.*cos(pi*n_mat*dx/Lx).*sin(pi*m_mat*dy/Ly)...
						                       +  Snm_2_22.*cos(pi*n_mat*dx/Lx).*cos(pi*m_mat*dy/Ly)...
						                       ,'all');
            end %add_i = 1:PBC_i
            %end non piecewise-------------------------------------------------------------
            
        end %j = 1:subsety
        %fprintf("blah\n")
    end %k = 1:subsetx
%     toc
%     fprintf("page\n");
%     tic
    %calculate RMSE between MD and reconstructed
    RMSE_S11 = sqrt(sum((sigxx_2-MD_S11).^2,'all')/N);
    RMSE_S12 = sqrt(sum((sigxy_2-MD_S12).^2,'all')/N);
    RMSE_S22 = sqrt(sum((sigyy_2-MD_S22).^2,'all')/N);
%     toc
    %bc RMSE is for minimization but the code is more maximization,
    % we will maximize the negative fitness which is effectively
    % minimization
    %RMSE will always produce a 0 or positive value, so the max of negative
    %should work
    H(1,i) = -RMSE_S11;
    H(2,i) = -RMSE_S12;
    H(3,i) = -RMSE_S22;


end %for i = 1:pages
