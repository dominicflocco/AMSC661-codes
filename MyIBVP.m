% AMSC 661, Homework 10
% Dominic Flocco, April 19, 2023

function MyIBVP()
close all
% load mesh 
msh = load('ring_mesh.mat');
pts = msh.pts;
tri = double(msh.tri);
tri = tri + ones(size(tri));

% get dirichlet boundary indices
dirichletOuter = [];
dirichletInner = [];
tol = 1e-2;
for k = 1:size(pts,1)
    x = pts(k,1);
    y = pts(k,2);
    if abs((x^2 + y^2 - 4)) <= tol
        % pt on outer ring
        dirichletOuter = [dirichletOuter, k];
    elseif abs((x^2 + y^2 - 1)) <= tol
        % pt on inner ring
        dirichletInner = [dirichletInner, k];
    end
end
dirichlet = [dirichletInner, dirichletOuter]';
% call FEM
fem2d_heat(pts,tri,dirichlet);
end 

function fem2d_heat(pts,tri,dirichlet)
% Initialization
Npts = size(pts,1);
Ntri = size(tri,1);
FreeNodes=setdiff(1:Npts,unique(dirichlet));
A = sparse(Npts,Npts);
B = sparse(Npts,Npts); 
T = 1; dt = 0.01; N = T/dt;
U = zeros(Npts,N+1);

% Assembly
for j = 1:Ntri
	A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:)) + stima3(pts(tri(j,:),:));
end
for j = 1:Ntri
	B(tri(j,:),tri(j,:)) = B(tri(j,:),tri(j,:)) + ...
        det([1,1,1;pts(tri(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
end
% Initial Condition
U(:,1) = IC(pts); 
% time steps
for n = 2:N+1
    b = sparse(Npts,1);
    % Volume Forces
    for j = 1:Ntri
        b(tri(j,:)) = b(tri(j,:)) + ... 
            det([1,1,1; pts(tri(j,:),:)']) * ... 
            dt*myf(sum(pts(tri(j,:),:))/3,n*dt)/6;
    end
    
    % previous timestep
    b=b+(B-dt/2*A)*U(:,n-1);
    % Dirichlet conditions
    u = sparse(Npts,1);
    u(unique(dirichlet)) = myu_d(pts(unique(dirichlet),:),n*dt);
    b=b-(dt/2*A+B)*u;
    % Computation of the solution
    u(FreeNodes) = (dt/2*A(FreeNodes,FreeNodes)+ ...
            B(FreeNodes,FreeNodes))\b(FreeNodes);
    U(:,n) = u;
end
save('myIBVP_sol.mat')
figure 
subplot(1,2,1)
t = 0.1;
p = ceil(t/dt) + 1;
u = U(:,p);
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp')
title(sprintf('Solution at Time  t = %.1f\n',t),'Fontsize',14);
xlabel('x');
ylabel('y');
axis ij
colorbar
view(2)
set(gca,'Fontsize',14);

subplot(1,2,2)
t = 1;
p = ceil(t/dt) + 1;
u = U(:,p);
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp')
title(sprintf('Solution at Time t = %.1f\n',t),'Fontsize',14);
xlabel('x');
ylabel('y');
axis ij
colorbar
view(2)
set(gca,'Fontsize',14);


u = U(:,N+1); % N+1 corresponds to t=1.
r = sqrt(pts(:,1).^2 + pts(:,2).^2);
[rsort,isort] = sort(r,'ascend');
exact = (1 - rsort.^2)/4 + (3*log(rsort))/(4*log(2));
usort = u(isort);

hold on 
plot(rsort,exact,'Linewidth',2, 'DisplayName', 'Exact');
plot(rsort,usort,'Linewidth',2, 'DisplayName', 'Numerical');
title('Solution as a Function of r');
xlabel('r');
ylabel('u(r)');
legend;
grid;
set(gca,'Fontsize',14);
hold off
end
%%
function M = stima3(vertices)
d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
end
%%
function heatsource = myf(x,t)
heatsource = ones(size(x,1),1);
end
%%
function DirichletBoundaryValue = myu_d(x,t)
DirichletBoundaryValue =  zeros(size(x,1),1);
end

%%
function u = IC(pts)
    r = sqrt(pts(:,1).^2 + pts(:,2).^2);
    cosPhi = pts(:,1)./r;
    u = r + cosPhi;
end
