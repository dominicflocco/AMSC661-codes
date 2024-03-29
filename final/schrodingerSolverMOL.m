%% AMSC 661 Final Exam Problem 2
% Method of Lines (D0^2+RK4)
% Author: Dominic Flocco 
% Date: May 17, 2023

function schrodingerSolverMOL
    %% Initialization 
    % Input parameters
    k0 = 10; s0 = 0.1;
    Tmax = 0.4; nPts = 256; 
    x = linspace(-20,20,nPts);
    dx = 40/nPts;
    dt = 10^(-4);
    tSteps = floor(Tmax/dt);
    
    %% Define Central Difference Scheme
    e = ones(nPts,1);
    A = spdiags([e -2*e e],-1:1, nPts,nPts);
    A(nPts,1) = 1; A(1,nPts) = 1;
    A = -(1j/(2*dx^2))*A;
    % Initialize Solution
    u = zeros(nPts, tSteps);
    u(:,1) = initial(x,s0,k0)';
    %% Time March
    for k = 1:tSteps-1
        % Runge Kutta 4th order stages
        k1 = A*(u(:,k)); 
        k2 = A*(u(:,k) + 0.5*dt*k1);
        k3 = A*(u(:,k) + 0.5*dt*k2); 
        k4 = A*(u(:,k) + dt*k3); 
        % Compute Solution 
        u(:,k+1) = u(:,k) + dt*(k1+2*k2+2*k3+k4)/6;
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
    sgtitle('Solution to Schr\"{o}dinger Eq. using $$D_0^2$$-RK4 with $$N_x = 256$$','interpreter','latex','Fontsize',22,'FontWeight','Bold');
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