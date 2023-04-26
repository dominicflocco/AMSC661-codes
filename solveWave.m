% AMSC Homework 11 
% Dominic Flocco; April 26, 2023
function solveWave()
    method = 3;

    h = 0.05;
    a = sqrt(2);
    xx = [-6 : h : 6];
    n = length(xx) -1;
    x = xx(1 : n);
    tmax = 4/a;
    N = ceil(a*tmax/h);
    k = tmax/N;
     % solve for 0 <= t <= tmax
    
    xi = advection(-a,k,h,tmax, 4); %right
    eta = advection(a,k,h,tmax, 3); %left
    ut = a*xi - a*eta; 
    ux = xi + eta;
    sol = cumtrapz(h,ux(2:n,:));
    
    
    ts = [1/(2*a),1/a,2/a,4/a]; 
    figure
    for i = 1:4
        subplot(2,2,i)
        t = ts(i);
        j = floor(t/k);
        ex = exact(x,t,a)';
        exx = (1/2)*(dphi(x+a*t) + dphi(x-a*t));
        hold on
        plot(x(1:n-1), sol(:,j+1)', 'Linewidth', 2,'DisplayName', 'Numerical u');
        %plot(x, ux(:,j)', 'DisplayName', 'Numerical ux');
        title(sprintf('Time = %.1f\n',t),'Fontsize',14)
        xlabel('x','Fontsize',12)
        ylabel('u','Fontsize',12)
        plot(x, ex', 'Linewidth', 2,'DisplayName', 'Exact');
        hold off
        legend;
        grid;
    
    end

    if method == 1 % central difference
        sgtitle('Central difference','fontweight','bold','Fontsize',20)
    end
    if method == 2 % Lax - Friedrichs
        sgtitle('Lax - Friedrichs','fontweight','bold','Fontsize',20)
    end
    if method == 3 % upwind, left
         sgtitle('Upwind','fontweight','bold','Fontsize',20)
    end
    if method == 5 % Lax-Wendroff
         sgtitle('Lax-Wendroff','fontweight','bold','Fontsize',20)
    end
    if method == 6 % beam-warming, left
        sgtitle('Beam - Warming','fontweight','bold','Fontsize',20)
    end
    if method == 8 % leap-frog
        st = 'Leap - frog';
        sgtitle('Leap - frog','fontweight','bold','Fontsize',20)
    end

end
function U = advection(a,k,h,tmax,method)
% Solves the advection equation u_t + a*u_x = 0
% on the interval [-6,6] with periodic BC u(0,t) = u(25,t)
% using different finite difference methods
if method == 1 % central difference
    st = 'Central difference';
end
if method == 2 % Lax - Friedrichs
    st = 'Lax - Friedrichs';
end
if method == 3 % upwind, left
     st = 'Upwind, left';
end
if method == 4 % upwind, right
     st = 'Upwind, right';
end
if method == 5 % Lax-Wendroff
     st = 'Lax-Wendroff';
end
if method == 6 % beam-warming, left
    st = 'Beam - Warming, left';
end
if method == 7 % beam-warming, rights
    st = 'Beam - Warming, right';
end
if method == 8 % leap-frog
    st = 'Leap - frog';
end
nu = (a*k)/h; 
xx = [-6 : h : 6];
n = length(xx) -1;
x = xx(1 : n);
n = length(x);
N = ceil(abs(a)*tmax/h);
k = tmax/N;

U = zeros(n, N);
U(:,1) = (1/2)*dphi(x);
% clf; hold on; grid;
% hplot1 = plot(x,u,'linewidth',2,'color','r');
% hplot2 = plot(x,u,'linewidth',2,'color','k');
% axis([-6,6,-0.5,1.5]);
% set(gca,'Fontsize',20);
% title(st,'Fontsize',20);
t = k;
for j = 2:N+1
    u = U(:,j-1)';
    ujp1 = circshift(u,[0, -1]);
    ujm1 = circshift(u,[0, 1]);
    if method == 1 % central difference
        unew = u - 0.5*nu*(ujp1 - ujm1);
    end
    if method == 2 % Lax - Friedrichs
        unew = 0.5*(ujm1 + ujp1 - nu*(ujp1 - ujm1));
    end
    if method == 3 % upwind, left
        unew = u - nu*(u - ujm1);
    end
    if method == 4 % upwind, right neg
        unew = u - nu*(ujp1 - u);
    end
    if method == 5 % Lax-Wendroff
        unew = u - 0.5*nu*(ujp1 - ujm1 - nu*(ujp1 - 2*u + ujm1));
    end
    if method == 6 % beam-warming, left
        ujm2 = circshift(u,[0, 2]);
        unew = u - 0.5*nu*(3*u - 4*ujm1 + ujm2 - nu*(u - 2*ujm1 + ujm2));
    end
    if method == 7 % beam-warming, right
        ujp2 = circshift(u,[0, -2]);
        unew = u - 0.5*nu*(-3*u + 4*ujp1 - ujp2 - nu*(u - 2*ujp1 + ujp2));
    end
    if method == 8 % leap-frog
        if t < k
            uold = u;
            unew = u - 0.5*nu*(ujp1 - ujm1 + nu*(ujp1 - 2*u + ujm1));
        else
            unew = uold - nu*(ujp1 - ujm1);
            uold = u;
        end
    end
    t = t + k;
    U(:,j) = unew';
end

end

%%
function u = phi(x)
 % the initial condition
    u = max(1-abs(x), 0);
end
function ic = dphi(x)
    ic = zeros(length(x),1);
    for i = 1:length(x)
        if x(i) < 0 && x(i) > -1
            ic(i,1) = 1;
        elseif x(i) > 0 && x(i) < 1 
            ic(i,1) = -1;
        end
    end
end

function sol = exact(x,t,a)
    sol = (1/2)*(phi(x+a*t) + phi(x-a*t));
end
            


