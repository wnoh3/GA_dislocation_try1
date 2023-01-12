%objective: plot SINGLE dislocation using best constants

%LD Lx is subject to change which means you must adjust spec_plot_xlength
%and ylength
%the 2.4 in subsetx and subsety could also be subject to change
function [sigxx_2,sigxy_2,sigyy_2] = plot_fourier(Pnm_2_11,Qnm_2_11,Rnm_2_11,Snm_2_11,...
                        Pnm_2_12,Qnm_2_12,Rnm_2_12,Snm_2_12,...
                        Pnm_2_22,Qnm_2_22,Rnm_2_22,Snm_2_22,...
                        Lx,Ly,...
                        destDirectory)

    [term_n,term_m]=size(Pnm_2_11);
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


    %% plot 
    
    %Lx and Ly of reconstructed GB is set manually

    %the following is defined above--------------------------------------------
    %define new spectral plot values-------------------------------------------
    spec_plot_xlength = 345.6;
    spec_plot_ylength = 345.6;
    spec_plot_subsetx = round(2*spec_plot_xlength/2.4);%make sure this is an actual whole number!!!! 
    spec_plot_subsety = round(2*spec_plot_ylength/2.4);
    subx = 2*spec_plot_xlength/spec_plot_subsetx; % xdist value of each subsetx
    suby = 2*spec_plot_ylength/spec_plot_subsety; % ydist value of each subsety
    subArea = subx*suby;
    
    topleftxycoord = [-spec_plot_xlength ,spec_plot_ylength];

    sigxx_2          = zeros(spec_plot_subsety,spec_plot_subsetx); %make sure this is an actual whole number!!!! 
    sigxy_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    sigyy_2          = zeros(spec_plot_subsety,spec_plot_subsetx);
    
%     sigxx_2          = zeros(round(spec_plot_subsety),round(spec_plot_subsetx)); %make sure this is an actual whole number!!!! 
%     sigxy_2          = zeros(round(spec_plot_subsety),round(spec_plot_subsetx));
%     sigyy_2          = zeros(round(spec_plot_subsety),round(spec_plot_subsetx));
    
    % %traverse through subset Voronoi Spectral
    parfor j = 1:spec_plot_subsetx
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
    
            %end non piecewise-------------------------------------------------------------
    
        end %j = 1:subsety
    end %k = 1:subsetx

    %save plots
    %%
    figName= 'Optimal S11 Scaled Contours';
    f = figure('Name', figName ,'visible','on');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),sigxx_2,'CDataMapping','scaled');
    caxis([min(sigxx_2(:)),max(sigxx_2(:))]*0.1)
    colorbar
    axis equal
    saveas(f, [destDirectory,figName,'.png']);
    figName= 'Optimal S12 Scaled Contours';
    f = figure('Name', figName ,'visible','on');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),sigxy_2,'CDataMapping','scaled');
    caxis([min(sigxy_2(:)),max(sigxy_2(:))]*0.1)
    colorbar
    axis equal
    saveas(f, [destDirectory,figName,'.png']);
    figName= 'Optimal S22 Scaled Contours';
    f = figure('Name', figName ,'visible','on');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),sigyy_2,'CDataMapping','scaled');
    caxis([min(sigyy_2(:)),max(sigyy_2(:))]*0.1)
    colorbar
    axis equal
    saveas(f, [destDirectory,figName,'.png']);
    
    %%
    figName= 'Optimal S11';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),sigxx_2,'CDataMapping','scaled');
    caxis([min(sigxx_2(:)),max(sigxx_2(:))])
    colorbar
    axis equal
    saveas(f, [destDirectory,figName,'.png']);
    figName= 'Optimal S12';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),sigxy_2,'CDataMapping','scaled');
    caxis([min(sigxy_2(:)),max(sigxy_2(:))])
    colorbar
    axis equal
    saveas(f, [destDirectory,figName,'.png']);
    figName= 'Optimal S22';
    f = figure('Name', figName ,'visible','off');
    image(linspace(-spec_plot_xlength,spec_plot_xlength,spec_plot_subsetx),linspace(-spec_plot_ylength,spec_plot_ylength,spec_plot_subsety),sigyy_2,'CDataMapping','scaled');
    caxis([min(sigyy_2(:)),max(sigyy_2(:))])
    colorbar
    axis equal
    saveas(f, [destDirectory,figName,'.png']);

end %end function