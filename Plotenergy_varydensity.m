clc;
clear;
close all;

%% Define parameters
%% membrane properties
kappa = 25;  %% bending rigidity, Unit: k_BT
sigma = 0.012;  %% membrane tension, Unit: k_BT/nm^2
lambda = 0.0;   %% Line tension, Unit: k_BT/nm
c_0 = 0.06;   %% spontaneous curveture, Unit: nm^(-1)
%% Polymer properties
a = 10;      %% monomer size,  Unit: nm
Chi = 0;   %% Flory parameter
v = (1-2*Chi)*a^3;     %% Excluded volume parameter, Unit: nm^3

R_0 = 50;        %%   membrane patch radius,  Unit: nm
A_s = pi*R_0^2;   %%   membrane patch area,  Unit: nm^2

% N_t = linspace(15,50,36);
% N_t = [15 20 25 30 35 55];           %% Total number of monomers
N_t = 20;
t = 0;
u = 0;
for N = N_t
    N
    u = u+1;
    R_G = (1/sqrt(6))*a*N^(3/5);  %% radius of gyration, Unit: nm
    rho_c1 = 1/(pi*R_G^2);     %% Critical density, if rho_g < rho_c1, then polymer is in mushroom regime. else it is in brush regime, Unit: nm^(-2)
    R_F= a*N^(3/5);  %% Flory radius, Unit: nm
    rho_c2 = 1/(R_F^2);  %% Critical density, if rho_g < rho_c2, then polymer is in mushroom regime. Unit: nm^(-2)
    
%     xi_j = [10 11.1 12.9 15.8];       %% Here we let the blob size on the surface = grafting diatance, Unit: nm
%     xi_j = 10;
    xi_j = linspace(10,30,201);
    j = 0;
    for xi = xi_j
        xi
        j = j+1;
        S_R = xi^2;     %% local area per chain, Unit: nm^2
        rho_g = 1/xi^2;  %% grafting density, Unit: nm^(-2)
        rho_anchor = 1/xi^2;  %% grafting density, Unit: nm^(-2)
        N_p = A_s/S_R;   %% Number of polymer chain
        if xi<2*(1/sqrt(6))*a*N^(3/5)
            brush = 1;
        else
            brush = 0;
        end

        eta_span = 0:0.0005:4; %% fraction of which characterize the area ratio between a spherical cap and a sphere.
        i = 0;
        F_tot = 0;
        for eta = eta_span
            i = i+1;
            n_b = fix(eta);
            H_flat = N*S_R^(-1/3)*(v*a^2/3)^(1/3);
            F_flat1 = N_p*N*(9/2)*(v/3)^(2/3)*(S_R*a)^(-2/3);
            F_flat2 = (9/2)*(pi*R_0^2/xi^2)*N*(v/(3*a*xi^2))^(2/3);
            if eta==0
                H_brush(i) = N*S_R^(-1/3)*(v*a^2/3)^(1/3);
                Z_max(i) = 0;
                F_poly(i) = (9*N/2/kappa)*(R_0^2/xi^2)*(v/(3*a*xi^2))^(2/3);
                F_bend(i) = 8*(sqrt(eta)-R_0*c_0/4)^2;
                F_tension(i) = sigma/kappa*eta*R_0^2;
                F_line(i) = 2*lambda/kappa*R_0*sqrt(1-eta);
                F(i) = (9*N/2/kappa)*(R_0^2/xi^2)*(v/(3*a*xi^2))^(2/3)+8*(sqrt(eta)-R_0*c_0/4)^2+sigma/kappa*eta*R_0^2+2*lambda/kappa*R_0*sqrt(1-eta);
            else
                R = R_0/(2*sqrt(eta));
                H_spherical = R*(1+5/3*H_flat/R)^(3/5)-R;
                H_brush(i) = R*(1+5/3*H_flat/R)^(3/5)-R;
                Z_max(i) = R_0*sqrt(eta);
                F_poly(i) = (9/2/kappa)*(R_0^2/xi^2)*(R_0/(2*sqrt(eta)))*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+10*N*sqrt(eta)/3/R_0*(v*a^2/3/xi^2)^(1/3))^(1/5)-1);
                F_bend(i) = 8*(sqrt(eta)-R_0*c_0/4)^2;
                F_tension(i) = sigma/kappa*R_0^2*(n_b+(eta-n_b)^2)/eta;
                F_line(i) = 2*lambda/kappa*R_0*sqrt((eta-n_b)*(1-eta+n_b)/eta);
                F(i) = (9/2/kappa)*(R_0^2/xi^2)*(R_0/(2*sqrt(eta)))*(3*v^(1/2)/(xi*a^2))^(2/3)*((1+10*N*sqrt(eta)/3/R_0*(v*a^2/3/xi^2)^(1/3))^(1/5)-1)...
                    +8*(sqrt(eta)-R_0*c_0/4)^2+sigma/kappa*R_0^2*(n_b+(eta-n_b)^2)/eta+2*lambda/kappa*R_0*sqrt((eta-n_b)*(1-eta+n_b)/eta);
                F_tot = F_tot+exp(-F(i));
            end
        end
        %% Energy minimum
        [y,index_min] = min(F);
        eta_min(j) = eta_span(index_min);
        % eta_min(t) = eta_span(index_min);
        if eta_span(index_min) == 0
            Budding = 0;
        end
        if 0<eta_span(index_min) && eta_span(index_min)<=1/2
            Budding = 0.4;
        end
        if 1/2<eta_span(index_min) && eta_span(index_min)<1
            Budding = 0.6;
        end
        if eta_span(index_min)==1
            Budding = 1;
        end
        if 1<eta_span(index_min) && eta_span(index_min)<2
            Budding = 1.5;
        end
        if eta_span(index_min)>=2
            Budding = 2;
        end
        F = F';
        F_poly = F_poly';
        F_bend = F_bend';
        F_tension = F_tension';
        F_line = F_line';
        F_force = F_force';
        E{j} = F;
        leg{j}=strcat('\xi=',num2str(xi),'nm');
        t = t+1;
        % p=phase;
        Results(t,:) = [xi,rho_anchor,N_p,N,eta_span(index_min),y,Budding,H_brush(index_min),Z_max(index_min),...
            kappa,sigma,lambda,c_0,a,R_0,brush];
    end
end

eta_span = eta_span';
subplot(2,2,1);
for j=1:length(xi_j)
    plot(eta_span,E{j},'linewidth',2)
    hold on
end
hold off
legend(leg)
xlabel('\eta');
ylabel('\itF_{tot}/k_BT')
str = {'\kappa=25 k_BT','\sigma=0.012 k_BT/nm^2','\lambda=0.5 k_BT/nm','c_0=0.04 nm^{-1}','a=10 nm','N=20','R_0=50 nm'};
text(0.7,3000,str)

subplot(2,2,2);
% plot(N_j,eta_min,'--o',xi_j,Results(:,5),'r',xi_j,Results(:,6),'b','linewidth',2)
hl1 = line(xi_j,eta_min,'Color','b','linewidth',2);
ax1 = gca;
xlabel('\it\xi');
ylabel('\it\eta_{min}');
ax1.XColor = 'b';
ax1.YColor = 'b';

subplot(2,2,3);
%% plot the spherical cap
% memColor = [0.9139    0.7258    0.3063];
memColor = [1.0    0.5    0.0];
for u=1:length(xi_j)
    if eta_min(u) == 0
        alpha = 0;
        Rb = 0;
        x = Rb:Rb+40;
        y = 0*x;
        plot(x,y,'color','b','linewidth',3)
        hold on;
        plot(-x,y,'color','b','linewidth',3)
        xlim([-20 20])
        ylim([-20 20])
        axis equal
    else
        n_b = fix(eta_min(u));
        % if n_b == eta_min(end)
        %     alpha = acos(-1-2*eta_min(end)+2*n_b);
        % else
        %     alpha = acos(1-2*eta_min(end)+2*n_b);
        % end
        alpha = acos(1-2*eta_min(u)+2*n_b);
        Rp = R_0/(2*sqrt(eta_min(u)));
        Rb = R_0*sqrt((eta_min(u)-n_b)*(1-eta_min(u)+n_b)/eta_min(u));
        Rb2 = Rp*sin(alpha);
        x = Rb:Rb+40;
        y = 0*x;
        for bead = 1:1:n_b
            h_b = 2*bead*Rp-Rp+Rp*(1-cos(alpha));
            z_b = 2*bead*Rp+Rp*(1-cos(alpha));
            theta1 = linspace(0,2*pi,1001);
            cir_x1 = Rp * cos(theta1);
            % cir_y1 = Rp * sin(theta1) - Rp * cos(alpha) + h_b;
            cir_y1 = Rp * sin(theta1) + h_b;
            % fill(cir_x1, cir_y1, [0.4940 0.1840 0.5560])
            plot(cir_x1, cir_y1,'Color', memColor, 'LineWidth',4)
            plot(cir_x1, cir_y1, 'LineWidth',4)
            hold on;
        end
        theta2 = linspace(pi/2-alpha,pi/2+alpha,101);
        cir_x2 = Rp * cos(theta2);
        cir_y2 = Rp * sin(theta2) - Rp * cos(alpha);
        plot(x,y,'color','b','linewidth',3)
        hold on;
        plot(-x,y,'color','b','linewidth',3)
        hold on
        % plot(cir_x2, cir_y2, 'Color', memColor, 'LineWidth',4)
        plot(cir_x2, cir_y2,  'LineWidth',4)
        xlim([-60 60])
        ylim([-20 60])
        axis equal
        axis off;
        % set(gcf,'color','none');
        % set(gca,'color','none');
    end
end
hold off
legend(leg)

subplot(2,2,4);
scatter(Results(:,3),Results(:,4),[],Results(:,5),'filled')
u=colorbar;
set(u,'FontName','Times New Roman','FontSize',15,'linewidth',2.5,'FontWeight','bold');
set(get(u,'title'),'string','\eta');
% u.Label.String = '\eta';
set(gca,'xscale','log');
xlim([8.5 81]);
xlabel('\itN_p');
ylabel('\itN')
box on