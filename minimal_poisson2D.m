% 2D Staggering example using a 2D Mimetic laplacian

clc
clear all
close all force

kkk = [2 4 6 8];
%kkk = [2];
mxmxmx = [64];
%mxmxmx = [512];
%omomom = [0.2 0.4 0.6 0.8 1.2 1.4 1.6 1.8];
%omomom = [0.2 0.4 0.6 0.8 1.2];
omomom = zeros(1, 15);

for iii = 1:size(kkk, 2)
    for jjj = 1:size(mxmxmx, 2)
        for ppp = 1:size(omomom, 2)
            % You may want to change these 2 lines with 
            % cd(FOLDER_TO_SCRIPT)
            folder_to_cover = fullfile(getenv('MSCPATH'), ...
                'csrc-mole/examples/matlab/');
            cd(folder_to_cover)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % GLOBAL PARAMETERS - DISCRETIZATION
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            k = kkk(1, iii); % Order of accuracy
            %h = 0.02; % Discretization size
            %m = 1/h; % Vertical resolution
            %n = 1/h; % Horizontal resolution
            mx = mxmxmx(1, jjj);
            ny = mx;
            n = mx + 2;  % Grid size (n x n), meaning (n-2)*(n-2) interior points
            disp('k');
            disp(k);
            disp('mx');
            disp(mx);
            h = 1/mx;
            hy = 1/ny;
            
            % GLOBAL PARAMETERS - GMRES
            tol = 1e-12; %1e-15: sensible to round-off errors
            maxit = 1000;
            
            % Robin parameters, alpha = 1 beta = 0 for Dirichlet;
            alpha = 1; 
            beta = 0;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TURN THESE ON/OFF FOR SOME FEATURES (PLOT FIGURES, SAVE FILES, ETC.)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            problema = 2; % 1: Laplace, analytical solution is a Fourier series, MOLE
                          % 2: Laplace, synthetic solution, Ford
                          % 3: Poisson, Robin conditions, 2022
                          % 4: Laplace, u(x,y)=sin(pi*x)*sin(pi*y)
            
            is_discretized = 0; % Always in 1, except if we need to change
                                % the linear system
            is_preconditioned = 1;
            compute_autovals = 0; compute_condest = 0;
            use_iterative_methods = 1;
            lgmress = ~use_iterative_methods;
            plot_figures = use_iterative_methods;
            plot_figures = 0;
            save_figures = 0;
            save_experiment = 1;
            save_linear_systems = 0;
            finite_diff = 0;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % GLOBAL PARAMETERS - PRECONDITIONING
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Uncomment only once at a time
            preconditioner = 'none';
            % preconditioner = 'Jacobi';
            % preconditioner = 'Gauss_Seidel';
            omega = omomom(1, ppp); % Only useful for SOR
            % preconditioner = strcat('SOR_omega_', num2str(omega));
            % preconditioner = 'ILU';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            % if(strcmp(preconditioner, 'none')) == 1
            %     % Will load existing non-preconditioned system
            %     linear_system = strcat('minimal_poisson2D_k_',num2str(k),...
            %         '_m_', num2str(mx),'.mat')';
            % else
            %     % Will load existing preconditioned system
            %     % First experiments will need to save this system before 
            %     % its recycling
            %     linear_system = strcat('minimal_poisson2D_k_',num2str(k),...
            %         '_m_', num2str(mx), '_', preconditioner, '.mat');
            % end
            
            % Eigenvalues of non-preconditioned matrix
            % autovals = eig(full(L));
            % autovName = strcat('eig_k_', num2str(k), '_m_', num2str(mx),...
            %         '_none.mat');
            % save(autovName, "autovals");
            
            
            % fig8 = figure(8);
            % plot(autovals, '+');
            % hold on;
            % figName = strcat('eig_k_', num2str(k), '_m_', num2str(mx),...
            %         '_none.fig');
            % savefig(fig8, figName);
            
            
            %cd('../');
            
            
            
            % fig7 = figure(7);
            % spy(L);
            % hold on;
            % figName = strcat('sparsity_k_', num2str(k), '_m_', num2str(mx),...
            %         '_none.fig');
            % savefig(fig7, figName);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % DISCRETIZATION - MUST CHANGE GLOBAL PARAMETERS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if is_discretized == 0
                folder_to_mole = fullfile(getenv('MSCPATH'), ...
                                            'csrc-mole/src/matlab/');
                cd(folder_to_mole)
                
                RHS2 = zeros(mx+2, ny+2);
                
                % Right-hand side
                if problema == 1
                    RHS2(1, :) = 100; % Known value at the bottom boundary
                elseif problema == 2
                    % New source term
                    v = [0 h/2:h:1-h/2 1];
                    w = sin(pi.*v);
                    w(1, end) = 0;
                    RHS2(1, :) = w;
                elseif problema == 3
                    lambda = +1;
                    alpha = -exp(lambda);
                    beta = (exp(lambda) - 1)/lambda;
                    F  = @(x) lambda*lambda*exp(lambda*x)/(exp(lambda)-1);
                    v = linspace(0,1,mx+2);
                    w = F(v);
                    for i = 1:n
                        RHS2(i,:) = w;
                    end
                else
                    F = @(x, y) -2.*pi.*pi.*sin(pi.*x).*sin(pi.*y);
                    x=[0 h/2:h:1-h/2 1]; % Staggered grid
                    y=[0 h/2:h:1-h/2 1];
                    [X, Y]=meshgrid(x, y);
                    RHS2 = h*h*F(X,Y);
                    RHS2(1,:) = 0;
                    RHS2(end, :) = 0;
                    RHS2(:, 1) = 0;
                    RHS2(:, end) = 0;
                end
            
                L = lap2D(k, mx, 1, ny, 1); % 2D Mimetic laplacian operator
                L = L + robinBC2D(k, mx, 1, ny, 1, alpha, beta);
            
                RHS = reshape(RHS2, [], 1);
            
            else
                % LOADING ALREADY DISCRETIZED LINEAR SYSTEMS
                % You may want to change these 2 lines with 
                % cd(FOLDER_TO_LIN_SYS)
                lin_sys_folder = fullfile(getenv('MSCPATH'), ...
                    'csrc-mole/examples/matlab/linear_systems');
                cd(lin_sys_folder)
                
                linear_system = strcat('minimal_poisson2D_k_',num2str(k),...
                        '_m_', num2str(mx),'.mat')';
            
                load(linear_system);
            end

            if compute_condest == 1
                CondNumber = condest(L);
	    else
		CondNumber = 0;
	    end

            disp('Condition Number before prec.');
            disp(CondNumber);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ANALYTICAL SOLUTION
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % TO BE CHANGED in terms of whether we solve by finite differences or by 
            % mimetic operators
            
            if use_iterative_methods == 1 || lgmress == 1
                x=[0 h/2:h:1-h/2 1]; % Staggered grid
                y=[0 h/2:h:1-h/2 1];
            end
            
            if finite_diff == 1
                x=[0:h:1]; % Finite difference grid
                y=[0:h:1];
            end
            
            [X, Y]=meshgrid(x, y);
            uxy=zeros(size(X));
            
            zzz = zeros(size(X,1),1);
            
            if problema == 1
                nFuncs = 220;
                for nn=1:nFuncs
                    cnn = (1-(-1)^nn)/(nn*sinh(nn*pi));
                    uxy = uxy + cnn*sin(nn*pi*X).*sinh(nn*pi*Y);
                end
                uxy=(200/pi).*uxy;
                uxy(1,:) = 0;
                uxy(:,1) = 0;
                uxy(:, n) = 0;
                uxy(n,:) = 100;
            elseif problema == 2
                uxy= (1/sinh(pi)).*sin(pi.*X).*sinh(pi.*Y);
            elseif problema == 3
                uxy= (exp(lambda.*X)-1)/(exp(lambda)-1);
            else
                uxy = sin(pi.*X).*sin(pi.*Y);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LEFT-PRECONDITIONING
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            timePREC = 0;
            
            disp(strcat('PRECONDICIONADOR = ', preconditioner));
            
            % If there is preconditioning
            if is_preconditioned == 1
            
                % Jacobi
                if strcmp(preconditioner, 'Jacobi') == 1
                    tic
                    M_jacobi = (diag(diag(L)))^-1;
                    L = M_jacobi * L;
                    RHS = M_jacobi * RHS;
                    timePREC = toc();
                    disp(timePREC);
                end
                
                % Gauss-Seidel
                if strcmp(preconditioner, 'Gauss_Seidel') == 1
                    tic
                    M_GS = (tril(L))^-1;
                    L= M_GS*L;
                    RHS = M_GS*RHS;
                    timePREC = toc();
                    disp(timePREC);
                end
                
                % SOR
                if strcmp(preconditioner, strcat('SOR_omega_', num2str(omega))) == 1
                    %disp('omega');
                    %disp(omega);
                    tic
                    D=diag(diag(L));
                    L2=-(tril(L)-D);
                    M_sorr=(1/omega)*D-L2;
                    M_sor=M_sorr^-1;
                    L= M_sor*L;
                    RHS = M_sor*RHS;
                    %L=M_sor\L;
                    %RHS=M_sor\RHS;
                    timePREC = toc();
                    disp(timePREC);
                end
                
                % ILU
                if strcmp(preconditioner, 'ILU') == 1
                    tic
                    setup.type = 'ilutp';
                    setup.droptol = 1e-6;
                    setup.udiag=1;
                    
                    [LT,UT] = ilu(sparse(L), setup);
                    M=LT * UT;
                    L=M\L;
                    RHS=M\RHS;
                    timePREC = toc();
                    disp(timePREC);
                end
                
                % Estimating condition number of coefficient matrix AFTER preconditioning
                if compute_condest == 1
		    CondNumber = condest(L);
		else
		    CondNumber = 0;
		end
                disp('Condition Number after prec.');
                disp(CondNumber);
                
                if plot_figures == 1
                    fig9 = figure(9);
                    spy(L);
                    hold on;
                    if save_figures == 1
                        figName = strcat('sparsity_k_', num2str(k), '_m_', num2str(mx),...
                                '_', preconditioner, '.fig');
                        savefig(fig9, figName);
                    end
                end
            
            else
                % NO Preconditioning
                disp(timePREC);
            end
            
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SPECTRAL ANALYSYS - EIGENVALUES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if compute_autovals == 1
                nEigs = 6;
            
                %n largest-magnitude eigenvalues of coefficient matrix
                max_autovals = eigs(L, nEigs, 'largestabs');
                
                %n smallest-magnitude eigenvalues of coefficient matrix
                min_autovals = eigs(L, nEigs, 'smallestabs');
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % NUMERICAL SOLUTION - BACKSLASH
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('BACKSLASH');
            tic;
            SOLB = L\RHS;
            timeB = toc();
            SOLBRS = reshape(SOLB, mx+2, []);
            disp(timeB);
            relresB = norm(RHS - L*SOLB)/norm(RHS);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % KRYLOV SUBSPACE METHODS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if use_iterative_methods == 1
            
                % You may want to change these 2 lines with 
                % cd(FOLDER_TO_SOLVERS)
                solver_folder = fullfile(getenv('MSCPATH'), ...
                    'krysbas-dev/test-prec-gmres/src/solvers');
                cd(solver_folder)
                
                % BUILT-IN GMRES
                % m = 30;
                % [SOL0, flag0, relres0, iter0, resvec0] = ...
                %     gmres(L, RHS, m, tol, maxit);
                
                % BUILT-IN GMRES + JACOBI
                % [SOL00, flag00, relres00, iter00, resvec00] = ...
                %     gmres(L, RHS, m, tol, maxit, diag(diag(L)));
                
                % PREC_GMRES
                % disp('MI-PREC-GMRES');
                % m = 30;
                % print = 1;
                % color = 'c-';
                % flagprec = 1; 
                % itermax = 1000; 
                % %itermax = maxit;
                % Name_Matrix = 'poisson';
                % tic
                % [vec11, vny1, xx1, relresPG]= mi_prec_GMRES_m(L,RHS,m,itermax,tol, color,print, Name_Matrix, flagprec);
                % timePG = toc();
                % disp(timePG);
                
                % PD_GMRES
                mInitial = 30;
                disp('PD-GMRES');
                [SOLPD, flagPD, relresvecPD, mvec, timePD] = ...
                    pd_gmres(L, RHS, mInitial, [], [], tol, maxit, [], [1, 1]);
                iternPD = sum(mvec);
                disp(timePD);
                SOLPDRS = reshape(SOLPD, mx + 2, ny + 2);
                
                % LGMRES
                m = 27;
                kk = 3;
                disp('LGMRES');
                [SOLL, flagL, relresvecL, kdvecL, timeL] = ...
                    lgmres(L, RHS, m, kk, tol, maxit);
                iternL = sum(kdvecL);
                disp(timeL);
                SOLLRS = reshape(SOLL, mx + 2, ny + 2);
                
                % GMRES-E
                m = 27;
                d = 3;
                disp('GMRES-E');
                [SOLE, flagE, relresvecE, kdvecE, timeE] = ...
                    gmres_e(L, RHS, m, d, tol, maxit);
                iternE = sum(kdvecE);
                disp(timeE);
                SOLERS = reshape(SOLE, mx + 2, ny + 2);
                
                % Parameters for built-in MATLAB iterative methods
                m = 30;
                maxcycles = maxit; % maxcycles: max number of outer iterations
                maxit = m * maxcycles; % maxit: max number of inner iterations
                                        % necessary as built-in gmres() input parameter 
                
                % GMRES(m)
                disp('GMRES(m)');
                tic;
                [SOLG, flagG, resG, iterG, resvecG] = gmres(L, RHS, m, tol, maxcycles);
                iternG = (iterG(1,1)-1)*m + iterG(1,2);
                relresvecG = resvecG ./ norm(RHS);
                relresvecGc = relresvecG(1:m:size(relresvecG, 1), 1);
                timeG = toc();
                disp(timeG);
                SOLGRS = reshape(SOLG, mx + 2, ny + 2);
                
                % GMRES(m) + PRECONDITIONER
                % disp('PREC-GMRES(m), GAUSS-SEIDEL');
                % tic;
                % [SOLGP, flagGP, resGP, iterGP, resvecGP] = gmres(L, RHS, m, tol, maxcycles, tril(L));
                %                                                 % Jacobi: M = diag(diag(L));
                %                                                 % GS: M = tril(L);
                % relresvecGP = resvecGP ./ norm(RHS);
                % timeGP = toc();
                % disp(timeGP);
                % SOLGPRS = reshape(SOLGP, mx + 2, ny + 2);
                
                % % BICGSTAB
                % disp('BICG')
                % tic
                % [SOL4, flag4, res4 , ITER4, resvec4] = bicg(L, RHS, tol, maxit);
                % relresvec4 = resvec4 ./ norm(RHS);
                % time4 = toc();
                % disp(time4)
                % SOL4RS = reshape(SOL4, mx+2, ny+2);
            
            end
            
            if lgmress == 1
                
                % You may want to change these 2 lines with 
                % cd(FOLDER_TO_SOLVERS)
                solver_folder = fullfile(getenv('MSCPATH'), ...
                    'krysbas-dev/test-prec-gmres/src/solvers');
                cd(solver_folder)
            
                % LGMRES
                m = 27;
                kk = 3;
                disp('LGMRES');
                [SOLL, flagL, relresvecL, kdvecL, timeL] = ...
                    lgmres(L, RHS, m, kk, tol, maxit);
                disp(timeL);
                SOLLRS = reshape(SOLL, mx + 2, ny + 2);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FINITE-DIFFERENCE METHOD
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            if finite_diff == 1
                
                % Parameters
                n = mx + 1;  % mx cells, meaning n = mx + 1 1d points, 
                             % meaning (n-2)*(n-2) interior points
                hh = 1/mx;  % Grid spacing
                
                % Initialize the coefficient matrix A and right-hand side vector b
                N = (n-2)^2;  % Number of interior points
                AA = sparse(N, N);  % Use sparse matrix for efficiency
                bb = zeros(N, 1);
            
                x = linspace(0,1,mx+1);
                
                % Assemble the system Ax=b
                for i = 1:(n-2)
                    for j = 1:(n-2)
                        % Current point index in the linear system
                        idx = (i-1)*(n-2) + j;
                        
                        % Diagonal term
                        AA(idx,idx) = -4;
                        
                        % Connect to left neighbor
                        if j > 1
                            AA(idx, idx-1) = 1;
                        end
                        
                        % Connect to right neighbor
                        if j < n-2
                            AA(idx, idx+1) = 1;
                        end
                        
                        % Connect to bottom neighbor
                        if i < n-2
                            AA(idx, idx+n-2) = 1;
                        end
                        
                        % Connect to top neighbor
                        if i > 1
                            AA(idx, idx-(n-2)) = 1;
                        end
                        
                        % Modify b for boundary conditions
                        bb(idx) = 0;
                        
                        % Top boundary condition (u(x,1) = 100)
                %         if i == 1
                %             bb(idx) = bb(idx) - 100;
                %         end
                
                        % Bottom boundary condition (u(x,0) = 100)
                        if i == (n-2)
                            %bb(idx) = bb(idx) - 100;
                            bb(idx) = bb(idx) - sin(pi*x(j+1));
                        end
                        
                        % Other boundary conditions (u(x,0) = u(0,y) = u(1,y) = 0)
                        % are already handled by the zero entries in b
                    end
                end
                
                % Print the full system for small grids
                % if n <= 6
                %     fprintf('Coefficient Matrix A:\n');
                %     full(AA)
                %     fprintf('\nRight-hand side vector b:\n');
                %     bb
                % end
                
                % Solve the system
            
            
                % AA = sparse(AA);
                
                disp('FINITE-DIFF METHOD')
                tic
                % uu_interior = AA\bb;
                solver_folder = fullfile(getenv('MSCPATH'), ...
                    'krysbas-dev/test-prec-gmres/src/solvers');
                cd(solver_folder)
                uu_interior = lgmres(AA, bb, 27, 3, 1e-9, 100);
                timeBFD = toc();
                disp(timeBFD);
            
                % Reconstruct the full solution
                uu = zeros(n, n);
                % Fill interior points
                for i = 1:(n-2)
                    for j = 1:(n-2)
                        uu(i+1,j+1) = uu_interior((i-1)*(n-2) + j);
                    end
                end
                % Fill boundary points
                % uu(1,:) = 0;        % Bottom boundary
                uu(n,:) = sin(pi.*x);      % Top boundary
                % uu(:,1) = 0;        % Left boundary
                % uu(:,n) = 0;        % Right boundary
                
                % Create grid for plotting
                x=[0:h:1]; % Finite difference grid
                y=[0:h:1];
                [X, Y]=meshgrid(x, y);
                
                % Create 3D surface plot
                if plot_figures == 1
                    fig10 = figure(10);
                    %figure('Position', [100, 100, 800, 600]);
                    %surf(X, Y, u);
                    imagesc(flipud(uu));
                    shading interp;
                    %colormap jet;
                    %set(gca, 'CLim', [0 1]);
                    set(gca, 'YDir', 'Normal')
                    colorbar;
                    
                    
                    % Customize the plot
                    title('2D Laplace Equation Solution (Linear System)', 'FontSize', 14);
                    xlabel('x', 'FontSize', 12);
                    ylabel('y', 'FontSize', 12);
                    % zlabel('u(x,y)', 'FontSize', 12);
                    % view(45, 45);
                    
                    % Create contour plot
                    % figure('Position', [100, 100, 800, 600]);
                    % contourf(X, Y, u, 20);
                    % colormap jet;
                    % colorbar;
                    % title('Contour Plot of Solution', 'FontSize', 14);
                    % xlabel('x', 'FontSize', 12);
                    % ylabel('y', 'FontSize', 12);
                    
                    if save_figures == 1
                        figName = strcat('fdm_k_', num2str(k), '_m_', num2str(mx),...
                                '_prec_', preconditioner, '_', datestr(now,30), '.fig');
                        savefig(fig10, figName);
                    end
                
                end
                
                % Print solution statistics
                fprintf('\nSolution Statistics:\n');
                fprintf('Maximum value: %f\n', max(uu(:)));
                fprintf('Minimum value: %f\n', min(uu(:)));
                fprintf('System size: %d x %d\n', N, N);
                fprintf('Condition number of A: %e\n', condest(AA));
            
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % FIGURES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % if plot_figures == 1      
            %     figure(1)
            %     imagesc(flipud(uxy))
            %     title('2D Poisson''s equation - analytical solution')
            %     xlabel('x grid point')
            %     ylabel('y grid point')
            %     colorbar
            %     %set(gca, 'CLim', [0 100])
            %     set(gca, 'YDir', 'Normal')
            %     hold on;
            % end
            if plot_figures == 1      
                figure(1)
                imagesc(flipud(uxy))
                title('2D Poisson''s equation - analytical solution')
                xlabel('x grid point')
                ylabel('y grid point')
                colorbar
                %set(gca, 'CLim', [0 100])
                set(gca, 'YDir', 'Normal')
                hold on;
                
                figure(2)
                % CHANGE PARAMETER 'SOL' HERE
                % RS: rectangular grid
                % SOLBRS: backslash
                % SOLPDRS: PD-GMRES
                % SOLLRS: LGMRES
                % SOLERS: GMRES-E
                imagesc(SOLLRS)
                title('2D Poisson''s equation - numerical solution')
                xlabel('x grid point')
                ylabel('y grid point')
                set(gca, 'YDir', 'Normal')
                colorbar
                hold on;
                
                fig3 = figure(3);
                plot(mvec);
                hold on;
                titleMVEC = strcat('m(j)\_k\_', num2str(k), '\_m\_', num2str(mx),...
                    '\_prec\_', preconditioner);
                title(titleMVEC, 'FontSize', 16);
                if save_figures == 1
                    figName = strcat('m(j)_k_', num2str(k), '_m_', num2str(mx),...
                        '_prec_', preconditioner, '_', datestr(now,30), '.fig');
                    savefig(fig3, figName);
                end
            
                fig4 = figure(4);
                % GMRES(m)
                semilogy(relresvecGc, 'r', 'LineWidth', 2);
                hold on;
                % PD-GMRES
                semilogy(relresvecPD, 'm', 'LineWidth', 2);
                hold on;
                % LGMRES
                semilogy(relresvecL, 'g', 'LineWidth', 2);
                hold on;
                % GMRES-E
                semilogy(relresvecE, 'b', 'LineWidth', 2);
                titleRRN = strcat('k = ', num2str(k), ', m = ', num2str(mx),...
                    ', prec = ', preconditioner);
                title(titleRRN, 'FontSize', 16);
                
                % If maxit == 100
                % xlim([0 100]);
                % xticks([0 20 40 60 80 100]);
                % xticklabels({'0','20','40','60','80','100'})
                
                % If maxit == 1000
                xlim([0 1000]);
                xticks([0 200 400 600 800 1000]);
                xticklabels({'0','200','400','600','800','100'})
                
                % Tol: 1e-9
                % ylim([1e-9 1]);
                % yticks([1e-9 1e-6 1e-3 1]);
                % yticklabels({'10^{-9}', '10^{-6}', '10^{-3}', '1'});
                
                % Tol: 1e-15
                ylim([1e-15 1]);
                yticks([1e-15 1e-12 1e-9 1e-6 1e-3 1]);
                yticklabels({'10^{-15}' '10^{-12}' '10^{-9}', '10^{-6}', '10^{-3}', '1'});
                
            
                xlabel('Number of restarts');
                ylabel('||r_j||  / ||r_0 ||');
                legendG = strcat('GMRES(30), cycles = ', num2str(iterG(1,1)), ...
                    ', time = ', num2str(timeG), 's');
                legendL = strcat('LGMRES(27,3), cycles = ', num2str(size(relresvecL,1)), ...
                    ', time = ', num2str(timeL), 's');
                legendE = strcat('GMRES-E(27,3), cycles = ', num2str(size(relresvecE,1)), ...
                    ', time = ', num2str(timeE), 's');
                legendPD = strcat('PD-GMRES(30, \alpha_P=1, \alpha_D=1), cycles = ', num2str(size(mvec,1)), ...
                    ', time = ', num2str(timePD), 's');
                legend(legendG, legendL, legendE, legendPD);
                pbaspect([4 3 1]); % or ([16 9 1])
                hold on;
                
                if save_figures == 1
                    figName = strcat('relresvec_k', num2str(k), '_m_', num2str(mx),...
                        '_prec_', preconditioner, '_', datestr(now,30), '.fig');
                    savefig(fig4, figName);
                end
                
                    figure(5)
                    subplot(1,2,1)
                    imagesc(uxy)
                    title('2D Poisson''s equation - analytical solution')
                    xlabel('x grid point')
                    ylabel('y grid point')
                    colorbar
                    %set(gca, 'Clim', [0 100])
                    set(gca, 'YDir', 'Normal')
                    subplot(1,2,2)
                    imagesc(flipud(SOLLRS))
                    title('2D Poisson''s equation - numerical solution')
                    xlabel('x grid point')
                    ylabel('y grid point')
                    colorbar
                    %set(gca, 'Clim', [0 100])
                    set(gca, 'YDir', 'Normal')
                    hold on;
                
                if strcmp(preconditioner, "none") ~= 1
                    figure(6)
                    spy(L);
                    hold on;
                end
            
                fig8 = figure(8);
                % GMRES(m)
                semilogy(relresvecGc, 'r', 'LineWidth', 2);
                hold on;
                titleRRN = strcat('k = ', num2str(k), ', m = ', num2str(mx),...
                    ', prec = ', preconditioner);
                title(titleRRN, 'FontSize', 16);
                xlim([0 100]);
                xticks([0 20 40 60 80 100]);
                xticklabels({'0','20','40','60','80','100'})
                
                % Tol: 1e-9
                % ylim([1e-9 1]);
                % yticks([1e-9 1e-6 1e-3 1]);
                % yticklabels({'10^{-9}', '10^{-6}', '10^{-3}', '1'});
                
                % Tol: 1e-15
                ylim([1e-15 1]);
                yticks([1e-15 1e-12 1e-9 1e-6 1e-3 1]);
                yticklabels({'10^{-15}' '10^{-12}' '10^{-9}', '10^{-6}', '10^{-3}', '1'});
            
                xlabel('Number of restarts');
                ylabel('||r_j||  / ||r_0 ||');
                legendG = strcat('GMRES(30), cycles = ', num2str(iterG(1,1)), ...
                    ', tiempo = ', num2str(timeG), 's');
                hold on;
            
                if save_figures == 1
                    figName = strcat('relresvec_GMRESm_k', num2str(k), '_m_', num2str(mx),...
                        '_prec_', preconditioner, '_', datestr(now,30), '.fig');
                    savefig(fig8, figName);
                end
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SAVE NEW DISCRETIZED LINEAR SYSTEMS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if save_linear_systems == 1 && strcmp(preconditioner, 'none') == 0
                new_linear_system = strcat('minimal_poisson2D_k_',num2str(k),...
                    '_m_', num2str(mx), '_', preconditioner, '.mat');
                save(new_linear_system, "L", "RHS");
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % ERROR ANALYSIS
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % For mimetic operators
            
            if mx < 512
                if use_iterative_methods == 1 || lgmress == 1
                    %SOLLRS = SOLLRS(2:mx+1, 2:mx+1);
                    %uxy = uxy(2:mx+1, 2:mx+1);
                
                    %SOLLRS = SOLLRS(mx/2+1-4:mx/2+1-4+9, mx/2+1-4:mx/2+1-4+9);
                    %uxy = uxy(mx/2+1-4:mx/2+1-4+9, mx/2+1-4:mx/2+1-4+9);
                    % 
                    uxy = flipud(uxy);
                    diff = SOLLRS - uxy;
                    diff = reshape(diff, 1, []);
                    eMInf = norm(diff, Inf);
                    eM2 = norm(diff);
                    disp('eMInf');
                    disp(eMInf);
                end
                
                % For finite difference method
                % uu = uu(5:mx-3, 5:mx-3);
                % uxy = uxy(5:mx-3, 5:mx-3);
                if finite_diff == 1
                    diff = uu - uxy;
                    diff = reshape(diff, 1, []);
                    eMInf = norm(diff, Inf);
                    eM2 = norm(diff);
                end
            
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % SAVE EXPERIMENT VARIABLES
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if save_experiment == 1
                folder_for_exps = fullfile(getenv('MSCPATH'), ...
                    'high-order-mimetic-operators/experiments');
                cd(folder_for_exps)
                
                experiment_filename = strcat('experiment_k_',num2str(k),'_m_',...
                    num2str(mx),'_',preconditioner,'_', datestr(now,30),'.mat');
                save(experiment_filename, "problema", ...
                    "k", "mx", "preconditioner", "omega", "CondNumber", ... %"eMInf", "eM2", ...
                    "tol", "maxit", "relresB", "relresvecGc", "relresvecL", "relresvecE", "relresvecPD", "iterG",...
                    "iternG", "iternPD", "iternL", "iternE", "mvec",...
                    "flagG", "flagE", "flagL", "flagPD", ...
                    "timeB", "timeG", "timeL", "timeE", "timePD", "timePREC"); %,...
                    %"max_autovals", "min_autovals");
		pause(1);
            end
        end
    end
end

load handel
sound(y, Fs)