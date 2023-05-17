
function schrodingerSolverFT
    % define parameters 
    k0 = 10; s0 = 0.1;
    Tmax = 0.4;
    nPts = 4096; 
    x = linspace(-20,20,nPts);
    
    dx = 40/nPts;
    dt = 10^(-4);
    tSteps = floor(Tmax/dt);
    t = linspace(0,Tmax,tSteps);
    L = 40;
    u = zeros(nPts, tSteps);
    u(:,1) = initial(x,s0,k0)';
    f0 = fftshift(fft(u(:,1)'));
    freqs = 2*pi*[-nPts/2 : nPts/2 - 1]/L;


    for k =1:tSteps
        ft = f0.*exp(-1j*freqs.^2*t(k)/2);
        u(:,k) = ifft(ifftshift(ft));

    end
    ex = exact(x,0,s0,k0);
    figure;
    subplot(3,2,1);
    hold on;
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

    sgtitle('Solution to Schr\"{o}dinger Eq. using FT method with $$N_x = 4096$$','interpreter','latex','Fontsize',22,'FontWeight','Bold');
end

function sol = exact(x,t,s0,k0)
    A = (2*pi*s0^2)^(-1/4);
    aux = 1 + 1j*t/(2*s0^2);
    num = x.^2 - 4j*s0^2*k0*x +2j*s0^2*k0^2*t; 
    denom = 4*s0^2*aux;
    sol = A*exp(-num/denom)/sqrt(aux);
end
function packet = initial(x,s0,k0)
    A = (2*pi*s0^2)^(-1/4);
    packet = A*exp(-x.^2/(4*s0^2) +1j*k0*x);
end
