% geophysic_final.m
%
% This routine : 
% - Model heat conduction by solving the one-dimensional heat-conduction equation.
%
% Requires :
% - Liao_BC.m
% - Liao_suit.m
%
%   Primitive equation:
%   dT/dt = k/(rho*cp)*(d^2T/dz^2) + A/(rho*cp),
%   where d/dt, d/dz are partial deriviative.
%
%   # Simulation 1: 
%     B.C. (i) T=T0 at z=0 (ii) Q=-k*dT/dz=-Q0 at z=0
%     I.C. (i) T=0 from z=1 to z=end
%   # Simulation 2: 
%     B.C. (i) T=T0 at z=0 (ii) Q=-k*dT/dz=-Qd at z=d
%     I.C. (i) T=0 from z=1 to z=d-1
%   # Simulation 3:
%     B.C. (i) T(t,z) = T(t,z-1) at z=w (ii) T(t,z) = T(t,z+1) at z=-w
%     I.C. (i) T=T0 at t=0 for -w<=x<=w (ii) T=0 at t=0 for |x|>w
%
%----------------------------------------------------------
% Numerical Method:
% Forward-Difference + Leap-Frog for time & Central-Difference for sapce derivative.
%       ____________________________
%  zj+1 |__|__|__|__|__|__|__|__|__|
%   zj  |__|__|__|__|__|__|__|__|__|
%    .  |__|__|__|__|__|__|__|__|__|
%    .  |__|__|__|__|__|__|__|__|__|
%    .  |__|__|__|__|__|__|__|__|__|
%    .  |__|__|__|__|__|__|__|__|__|
%   z1  |__|__|__|__|__|__|__|__|__|
%   z0  |__|__|__|__|__|__|__|__|__|
%        t0 t1 t2  .  .  . . tm-1 tm
%
% (T(t+1,j)-T(t-1,j))/2/dt 
% = k/(rho*Cp)*(T(t,j+1)-2*T(t,j)+T(t,j-1))/(dz^2) + A/(rho*Cp)
%
%----------------------------------------------------------
%
clear
clc
close all
%% User Define Section
%%%%%%%%%% Start of User Input %%%%%%%%%%
% domain
zd = 10*10^0;                                   % Depth (m)
z0 = 0;                                         % Surface position (m)
dz = 0.01;                                      % Grid length in z direction (m)
dt = 20;                                        % Time step (s)
tottime = 3*24*3600;                            % Total model integration time (s)

% BC at top
BC_top = 'Mixed';                               % B.C. for top grid
                                                % = 'Dirichlet', Dirichlet B.C.
                                                % =   'Neumann', Neumann B.C.
                                                % =     'Mixed', Mixed B.C. (i.e.,Dirichlet+Neumann)
                                                % =     'Close', Close B.C. (i.e., zero-gradient)
                                                
BC_top_var = 'Mix';                             % Material at BC_top 
                                                % = 'Temperature', used if BC_top = 'Dirichlet'
                                                % =        'Flux', used if BC_top = 'Neumann'
                                                % =         'Mix', used if BC_top = 'Mixed'
                                                % =        'none', used if BC_top = 'Close'
                                                
BC_top_val = 'Q';                               % Value of BC_top_var
                                                % =     0,  or any value you want
                                                % =   'Q',  also can in character, but need designing it.
                                                %           See how I design it for simulation 1 at line 144. 
                                                % = 'none', used if BC_top_var = 'Close'
                                                
BC_top_type = 'Inflow';                         % The direction of material pass top grid.
                                                % = 'Inflow'
                                                % = 'Outflow'

% BC at bottom
BC_bot = 'Close';                               % B.C. for bottom grid
                                                % = 'Dirichlet'
                                                % = 'Neumann'
                                                % = 'Mixed'
                                                % = 'Close'

BC_bot_var = 'none';                            % Material at BC_bot
                                                % = 'Temperature'
                                                % = 'Flux'
                                                % = 'Mix'
                                                % = 'none'

BC_bot_val = 'none';                            % Value of BC_top_var
                                                % = '0'
                                                % = 'Qd' 
                                                % = 'none'

BC_bot_type = 'Outflow';                        % The direction of material pass bottom grid.
                                                % = 'Inflow'
                                                % = 'Outflow'
                                                
% I.C.
IC = ["T(1,:)=T0;"];                            % Initial Condition, T = T(t,z)

% physic
T0 = 0;                                         % Surface temperature (degC)
Td = 0;                                         % Bottom surface temperature (degC)
A  = 1.25*10^-6;                                % Radioactive heat generation (W/m^3)
k  = 1.7;                                       % Conductivity (W/m/degC)
Qd = 21*10^-3;                                  % Basal heat flow (W/m^2) 
rho = 2.3*10^3;                                 % Density (kg/m^3)
Cp = 10^3;                                      % Specific Heat (J/(kg*degC))

% output
outtime      = 3600;                            % Output interval (s)

expcase      = 1;                               % Simulation case, = 1 for case 1
                                                %                  = 2 for case 2
                                                %                  = 3 for case 3
                                                
if_animation = 1;                               % If make animation, = 0 for no
                                                %                  , = 1 for yes
                                                
ani_name = 'final_proj.gif';                    % Output animation name

% suit
simcase = 1;                                    % The "suit" I used for my simulation HW. It would override 
                                                % above options if user open it.
                                                % = 0, no
                                                % = 1, #Simulation 1
                                                % = 2, #Simulation 2
                                                % = 3, #Simulation 3
             
%%%%%%%%%% End of User Input %%%%%%%%%%
%% Program begin
if (simcase==1 || simcase==2 || simcase==3)
    [zd,z0,dz,dt,tottime,...
     BC_top,BC_top_var,BC_top_val,BC_top_type,BC_bot,BC_bot_var,BC_bot_val,BC_bot_type,IC,...
     T0,Td,A,k,Qd,rho,Cp,...
     outtime,expcase,if_animation,ani_name] = Liao_suit(simcase);
end

e_vert = (zd-z0)/dz +1;
z = linspace(z0,zd,e_vert);
nt = tottime/dt +1;

tt = linspace(0,tottime,nt);

%% Build Initial field
% Define Q for simulation 1
tt1 = linspace(0,345,7*3600/dt+1); tt2 = linspace(345,0,6*3600/dt+1); tt3 = linspace(0,0,11*3600/dt+1);
Q = [tt1(2:end),tt2(2:end),tt3(2:end)...
     tt1(2:end),tt2(2:end),tt3(2:end)...
     tt1(2:end),tt2(2:end),tt3(2:end)];

% Apply I.C.
T = zeros(nt,e_vert);
for ii = 1:length(IC)
    eval(IC(ii));
end

% Apply B.C.
if (ischar(BC_top_val)==1)
    if (strcmp(BC_top_val,'none')==0)
        BTV = eval(BC_top_val); % BTV = BC Top Value
    else
        BTV = 999;  % Arbitrary value if BC_top_val = 'none'
    end
else
    BTV = BC_top_val;
end
if (ischar(BC_bot_val)==1)
    if (strcmp(BC_bot_val,'none')==0)
        BBV = eval(BC_bot_val);
    else
        BBV = 999;
    end
else
    BBV = BC_bot_val;
end

if (max(size(BTV))>1); BTV2 = BTV(1);else; BTV2 = BTV;end       % if boundary condition is not scolor                                                      
T(1,1) = Liao_BC(BC_top,BC_top_var,BTV2,T0,dz,k,T(1,2),T(1,1)); % i.e., #simulation 1
if (max(size(BBV))>1); BBV2 = BBV(1);else; BBV2 = BBV;end
T(1,end) = Liao_BC(BC_bot,BC_bot_var,BBV2,Td,dz,k,T(1,end-1),T(1,end));
     
%% Solve time-stepping
for t = 1:nt-1
    for j = 1:e_vert

        if (j==1)
            if (length(BTV)>1); BTV2 = BTV(t);else; BTV2 = BTV;end
            T(t+1,j) = Liao_BC(BC_top,BC_top_var,BTV2,T0,dz,k,T(t,j+1),T(t,j));
        elseif (j==e_vert)
            if (length(BBV)>1); BBV2 = BBV(t);else; BBV2 = BBV;end
            T(t+1,j) = Liao_BC(BC_bot,BC_bot_var,BBV2,Td,dz,k,T(t,j-1),T(t,j));
        else
            if (t==1)
                s = dt*k/rho/Cp/(dz^2); c = A*dt/rho/Cp;
                T(t+1,j) = s*(T(t,j+1)+T(t,j-1)) + (1-2*s)*T(t,j) + c;
            else
                s = 2*dt*k/rho/Cp/(dz^2); c = 2*A*dt/rho/Cp;
                T(t+1,j) = s*(T(t,j+1)+T(t,j-1)) -2*s*T(t,j) + T(t-1,j) + c;
                T(t,j) = T(t,j) + 0.5*1.0*(T(t-1,j)-2.*T(t,j)+T(t+1,j)); % Apply filter to remove time splitting
            end
        end
        
    end
end
 
%}
 %% Plot and Print
 %
 if (expcase==1 || expcase==2)
  dd = 1;
    time_tick = ['05:00 AM'; '06:00 AM'; '07:00 AM'; '08:00 AM'; '09:00 AM'; '10:00 AM'; '11:00 AM'; '12:00 PM';...
                '13:00 PM'; '14:00 PM'; '15:00 PM'; '16:00 PM'; '17:00 PM'; '18:00 PM'; '19:00 PM'; '20:00 PM';...
                '21:00 PM'; '22:00 PM'; '23:00 PM'; '24:00 PM'; '01:00 AM'; '02:00 AM'; '03:00 AM'; '04:00 AM'];
    for t = 1:nt
        if (mod((t-1)*dt,outtime)==0)
            dd
            %[(t-1)*dt,dd];
            figure('visible','off')
            plot(T(t,:),z,'LineWidth',2)
            set(gca,'xaxislocation','top');
            set(gca,'YDir',' reverse ');
            xlim([0 ceil(max(max(T)))]); ylim([0 zd])
            %xlim([0 1700]); ylim([0 55*10^3]);
            %set(gca,'xtick',[0:500:1500,1700]);
            %set(gca,'ytick',[0:10*10^3:50*10^3]);
            %set(gca,'yticklabel',num2cell([0:10:50]));
            if (expcase == 1)
                dd2 = mod(dd,24);
                if (dd2 == 0)
                    dd2 = 24;
                end
                xlabel('Temperature (\circC)','FontSize',18); ylabel('Depth (m)','FontSize',18);
                title([time_tick(dd2,1:end)]);
            elseif (expcase==2)
                xlabel('Temperature (\circC)','FontSize',18); ylabel('Depth (m)','FontSize',18);
                title([num2str(dd-1),' Ma']);
            end

            grid on
            
            frames(dd) = getframe(gcf);
            dd = dd + 1;
        end
    end
 end
if (expcase==3)
    dd = 1;
    for t = 1:nt
        if (mod(tt(t),outtime)==0)
            dd
            figure('visible','off')
            plot(z(1:end),T(t,:),'LineWidth',2)
            set(gca,'xaxislocation','top');
            xlim([0 z(end)]); ylim([0 T0]);
            xlabel('Width (m)','FontSize',12); ylabel('Temperature (\circC)','FontSize',18);
            title([num2str(tt(t)/24/3600), 'day'],'FontSize',12)
    
            grid on
            frames(dd) = getframe(gcf);
            dd = dd + 1;
        end
    end
end
if (if_animation==1)
    for i=1:dd-1
        image=frame2im(frames(i));
        [im,cm]=rgb2ind(image,256);
        if i==1
            imwrite(im,cm,ani_name,'gif','LoopCount',Inf,'DelayTime',0.1);
        else
            imwrite(im,cm,ani_name,'gif','WriteMode','append','DelayTime',0.1);
        end
    end
end
 %}
     
