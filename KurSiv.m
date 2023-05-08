function KurSiv()
close all
fsz = 20; % fontsize
% solves u_t + u_{xxxx} + u_{xx} + (0.5u^2)_x = 0, i.e.,
% u_t = -u_{xxx} - u_{xx} - (0.5u^2)_x

N = 256;
L = 32*pi;
x = linspace(0,L,N+1);
x(N + 1) = [];
k = -N/2 : (N/2 - 1); % wave numbers
%
tmax = 200;
t = 0;
dt = 0.1; % time step
u = zeros(tmax/dt,N);
% initial data
u0 = cos(x/16).*(ones(1,N)+sin(x/16));

figure; clf; 
hpic = plot(x,u0,'LineWidth',2,'color','r'); % plot of the numerical solution
hold on;
grid
xlim([0,L]);
set(gca,'Fontsize',fsz);
xlabel('x','FontSize',fsz);
ylabel('u','FontSize',fsz);
drawnow

freq = k.*(2*pi/L); % frequencies
freq2 = (freq).^2;
freq4 = (freq).^4;
% in the Fourier space, uhat = e3.*vhat
e24 = exp((freq2-freq4)*dt);
u(1,:) = u0;
i=1;
while (t<tmax) 
    t=t+dt; i=i+1;
    vhat=fftshift(fft(u(i-1,:))); % v in the Fourier space
    % RK4 step in the Fourier space
    k1=rhs(0,vhat);
    k2=rhs(0.5*dt,vhat+0.5*dt*k1);
    k3=rhs(0.5*dt,vhat+0.5*dt*k2);
    k4=rhs(dt,vhat+dt*k3);
    vhat_new=vhat+dt*(k1+2*k2+2*k3+k4)/6;
    % return to the original space and the original variable u
    unew=ifft(ifftshift(e24.*vhat_new)); % return to u in the x-space
    set(hpic,'xdata',x,'ydata',real(unew));
    
    axis([0 L -0.01 2]);

    u(i,:)=unew;
    % drawnow
end
close all;
imagesc(real(u));
set(gca,'YDir','normal');
xlabel('x')

end
%%
function RHSvhat=rhs(dt,vhat)
% v should be a row vector
% RHSvhat = - e^{-tL}(1i*k*hat{(e^{tL}v)^2/2} 
N=size(vhat,2);
L = 32*pi;
k=-N/2 : (N/2 - 1);
freq = k.*(2*pi/L); % frequencies
freq2 = (freq).^2;
freq4 = (freq).^4;
e24 = exp((freq2-freq4)*dt);
em24=exp(-(freq2-freq4)*dt);
vhat1=vhat.*e24;          % e^{tL}v in the Fourier space 
v1=ifft(ifftshift(vhat1));      % exp(tL)v in the x-space
v2=0.5*v1.^2;          % [exp(tL)v]^2 in the x-space
RHSvhat=-em24.*(1i*freq).*fftshift(fft(v2)); % exp(-tL)[[(exp(tL)v)]_x] in the Fourier space
end
