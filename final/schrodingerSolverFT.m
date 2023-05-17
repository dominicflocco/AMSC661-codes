
%% AMSC 661 Final Exam Problem 2
% Method of Fourier Transforms
% Author: Dominic Flocco 
% Date: May 17, 2023


function schrodingerSolverFT
    %% Initialization
    % Input parameters 
    k0 = 10; s0 = 0.1;
    Tmax = 0.4; L = 20; nPts = 256; 
    
    dx = 2*L/nPts;
    dt = 10^(-4);
    tSteps = floor(Tmax/dt);
    % Initialize space and time intervals
    tSteps = floor(Tmax/dt);
    t = linspace(0,Tmax,tSteps);
    x = linspace(-L,L,nPts);

    % Initialize solution array 
    u = zeros(nPts, tSteps); 
    
    % Compute solution at time 0
    u(:,1) = initial(x,s0,k0)';
    % Fourier Transform at time 0
    f0 = fftshift(fft(u(:,1)'));

    % Initialize frequencies
    freqs = pi*[-nPts/2 : nPts/2 - 1]/L;
    
    %% Time March
    for k =1:tSteps
        % Compute Fourier mode at time t(k)
        ft = f0.*exp(-1j*freqs.^2*t(k)/2);
        % Perform IFT on Fourier mode to obtain solution 
        u(:,k) = ifft(ifftshift(ft));
    end
    
    % Compute integral of probability density over interval at time 0
    density = trapz(abs(u(:,1)).^2,x);
    fprintf('Density integral at time t = %d: %f\n', 0,abs(density));
    
    %% Plot Solution 
    figure;
    subplot(3,2,1);
    hold on;
    ex = exact(x,0,s0,k0);
    plot(x, abs(u(:,1)),'DisplayName','Numerical','LineWidth', 1,'Color','r');
    plot(x,abs(ex), 'DisplayName','Exact','LineWidth', 1,'Color','b');
    
    title('$$t=0$$', 'interpreter','latex','Fontsize',18);
    legend('interpreter','latex', 'FontSize',12)
    axis([-20,20,0,2]);
    grid;
    set(gca,'FontName','Times','fontsize',12);
    xlabel('$$x$$','FontSize',14,'interpreter','latex'); 
    ylabel('$$|\psi(x,t)|$$','FontSize',14,'interpreter','latex')
    hold off;
    for i = 1:5
        subplot(3,2,i+1)
        if i ==5
               title('$$t= T_{max}$$','interpreter','latex','Fontsize',18);
        else
               title(sprintf('$$t= (%d/5)T_{max}$$',i),'interpreter','latex','Fontsize',18);
        end
        t = Tmax*(i/5);
        n = floor(t/dt);
        ex = exact(x,t,s0,k0);

        % Compute integral of probability density over interval
        density = trapz(abs(u(:,n)).^2,x);
        fprintf('Density integral at time t = %.2f: %f\n', t,abs(density));
        
        % Plot solution 
        hold on;
        plot(x, abs(u(:,n)),'DisplayName','Numerical','LineWidth', 1,'Color','r');
        plot(x,abs(ex), 'DisplayName','Exact','LineWidth', 1,'Color','b');
        legend('interpreter','latex', 'FontSize',12)
        axis([-20,20,0,2]);
        grid;
        set(gca,'FontName','Times','fontsize',12);
        xlabel('$$x$$','FontSize',14,'interpreter','latex'); 
        ylabel('$$|\psi(x,t)|$$','FontSize',14,'interpreter','latex')
        hold off;
    end
    sgtitle('Solution to Schr\"{o}dinger Eq. using FT method with $$N_x = 256$$','interpreter','latex','Fontsize',22,'FontWeight','Bold');
end
%% Compute Exact Solution
function sol = exact(x,t,s0,k0)
    A = (2*pi*s0^2)^(-1/4);
    aux = 1 + 1j*t/(2*s0^2);
    num = x.^2 - 4j*s0^2*k0*x +2j*s0^2*k0^2*t; 
    denom = 4*s0^2*aux;
    sol = A*exp(-num/denom)/sqrt(aux);
end
%% Compute Initial Condition
function packet = initial(x,s0,k0)
    A = (2*pi*s0^2)^(-1/4);
    packet = A*exp(-x.^2/(4*s0^2) +1j*k0*x);
end
