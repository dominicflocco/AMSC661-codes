%% AMSC 661 Final Exam Problem 3
% Author: Dominic Flocco
% Modified from Appendix C of Alberty et al. 
% Date: May 17, 2023 

function mypb3
    %% Load in mesh 
    % Mesh generated in accompanying file generateMesh.ipynb
    msh = load('mesh.mat');
    pts = msh.pts;
    tri = double(msh.tri);
    tri = tri + ones(size(tri));

    %% Determine Dirichlet Points in Mesh
    dirichlet = [];
    tol = 10^(-8);
    for k = 1:size(pts,1)
        x = pts(k, 1);
        y = pts(k, 2);
        % Omega = [-1,1]^2
        if abs(x-1) <= tol || abs(y-1) <= tol || abs(x+1) <= tol || abs(y+1) <= tol
            dirichlet = [dirichlet; k];
        end
    end

    % Compute free nodes as mesh pts not on dirichlet boundary
    FreeNodes=setdiff(1:size(pts,1),unique(dirichlet));
    
    %% Initialize Solution u
    % u = -ones(size(pts,1),1);
    u = sign(pts(:,1));
    % Set dirichlet BC
    u(unique(dirichlet)) = uD(pts(unique(dirichlet),:));

    %% Newton-Raphson iteration
     for i=1:50
        % Assembly of DF(u)
        A = sparse(size(pts,1),size(pts,1));
        for j = 1:size(tri,1)
           A(tri(j,:),tri(j,:)) = A(tri(j,:), tri(j,:)) ...
               + localdj(pts(tri(j,:),:),u(tri(j,:)));
        end
        % Assembly of F(u)
        b = sparse(size(pts,1),1);
        for j = 1:size(tri,1)
            b(tri(j,:)) = b(tri(j,:)) + localj(pts(tri(j,:),:),u(tri(j,:)));
        end
        % Set Dirichlet BC for Newton iterate
        w = zeros(size(pts,1),1);
        % test funcion = 0 on Dirichlet Boundary
        W(unique(dirichlet)) = 0; 

        % Solve Newton step
        w(FreeNodes) = A(FreeNodes,FreeNodes)\b(FreeNodes);
        u = u - w;

        % Check if solution converged
        if norm(w) < 10^(-10)
            break
        end
     end


     % Plot solution
     trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp','LineWidth', 0.25);
     colormap(jet);
     set(gca,'FontName','Times','fontsize',14);
     xlabel('$$x$$', 'Fontsize', 16,'interpreter','latex')
     ylabel('$$y$$', 'Fontsize', 16,'interpreter','latex')
     zlabel('$$u(x,y)$$', 'Fontsize', 16,'interpreter','latex')
     ylim([-1 1]); xlim([-1 1]);
     title('Solution to Ginzburg-Landau Equation with $$u_D = $$sign$$(x)$$', 'Fontsize', 22,'interpreter','latex')
end

%% Set Dirichlet Boundary Conditions
function initial = uD(vrts)
    initial = zeros(size(vrts,1),1);
    % "Frustrated" BCs in part (e):
    tol =10^(-8);
    for i=1:size(vrts,1)
        % left and right boundaries 
        if abs(vrts(i,1) - 1) <= tol ||abs(vrts(i,1) +1) <= tol
            initial(i) = 1;
        end 
        % top and bottom booundaries
        if abs(vrts(i,2) - 1) <= tol || abs(vrts(i,2)+1) <= tol 
            initial(i) = -1;
        end 
    end
end
%% Compute local value of function F being minized
function b = localj(vertices,u)
    Eps = 1/100;
    G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
    Area = det([ones(1,3);vertices']) / 2;
    b=Area*((Eps*G*G'-[2,1,1;1,2,1;1,1,2]/12)*u+ ...
    [4*u(1)^3+ u(2)^3+u(3)^3+3*u(1)^2*(u(2)+u(3))+2*u(1) ...
    *(u(2)^2+u(3)^2)+u(2)*u(3)*(u(2)+u(3))+2*u(1)*u(2)*u(3);
    4*u(2)^3+ u(1)^3+u(3)^3+3*u(2)^2*(u(1)+u(3))+2*u(2) ...
    *(u(1)^2+u(3)^2)+u(1)*u(3)*(u(1)+u(3))+2*u(1)*u(2)*u(3);
    4*u(3)^3+ u(2)^3+u(1)^3+3*u(3)^2*(u(2)+u(1))+2*u(3) ...
    *(u(2)^2+u(1)^2)+u(2)*u(1)*(u(2)+u(1))+2*u(1)*u(2)*u(3)]/60);
end
%% Compute Local Jacobian DF matrix over triangle  
function M = localdj(vertices,u)
    Eps = 1/100;
    G = [ones(1,3);vertices'] \ [zeros(1,2);eye(2)];
    Area = det([ones(1,3);vertices']) / 2;
    M = Area*(Eps*G*G'-[2,1,1;1,2,1;1,1,2]/12 + ...
    [12*u(1)^2+2*(u(2)^2+u(3)^2+u(2)*u(3))+6*u(1)*(u(2)+u(3)),...
    3*(u(1)^2+u(2)^2)+u(3)^2+4*u(1)*u(2)+2*u(3)*(u(1)+u(2)),...
    3*(u(1)^2+u(3)^2)+u(2)^2+4*u(1)*u(3)+2*u(2)*(u(1)+u(3));
    3*(u(1)^2+u(2)^2)+u(3)^2+4*u(1)*u(2)+2*u(3)*(u(1)+u(2)),...
    12*u(2)^2+2*(u(1)^2+u(3)^2+u(1)*u(3))+6*u(2)*(u(1)+u(3)),...
    3*(u(2)^2+u(3)^2)+u(1)^2+4*u(2)*u(3)+2*u(1)*(u(2)+u(3));
    3*(u(1)^2+u(3)^2)+u(2)^2+4*u(1)*u(3)+2*u(2)*(u(1)+u(3)),...
    3*(u(2)^2+u(3)^2)+u(1)^2+4*u(2)*u(3)+2*u(1)*(u(2)+u(3)),...
    12*u(3)^2+2*(u(1)^2+u(2)^2+u(1)*u(2))+6*u(3)*(u(1)+u(2))]/60);
end