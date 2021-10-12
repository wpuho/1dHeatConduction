function [zd,z0,dz,dt,tottime,...
          BC_top,BC_top_var,BC_top_val,BC_top_type,BC_bot,BC_bot_var,BC_bot_val,BC_bot_type,IC,...
          T0,Td,A,k,Qd,rho,Cp,...
          outtime,expcase,if_animation,ani_name] = Liao_suit(simcase)
    
    if (simcase==1)
        % domain
        zd = 10*10^0;              
        z0 = 0;               
        dz = 0.01;          
        dt = 20; 
        tottime = 3*24*3600;

        BC_top = 'Mixed';
        BC_top_var = 'Mix';
        BC_top_val = 'Q';
        BC_top_type = 'Inflow';
        %
        BC_bot = 'Close'; 
        BC_bot_var = 'none'; 
        BC_bot_val = 'none';
        BC_bot_type = 'Outflow';

        IC = ["T(1,:)=T0;"];

        % physic
        T0 = 0;     
        Td = 1700;     
        A  = 1.25*10^-6;  
        k  = 2.5;         
        Qd = 10.5*10^-3;    
        rho = 2.3*10^3; 
        Cp = 10^3;

        % output
        outtime      = 3600; 
        expcase      = simcase;
        if_animation = 1; 
        ani_name = 'Simulation_1.gif';
        
    end
    if (simcase==2)
        % domain
        zd = 50*10^3;              
        z0 = 0;               
        dz = 5*10^3;          
        dt = 0.1*10^6*365*24*3600; 
        tottime = 100*10^6*365*24*3600;

        BC_top = 'Dirichlet';
        BC_top_var = 'Temperature';
        BC_top_val = 0;
        BC_top_type = 'Outflow';
        %
        BC_bot = 'Neumann'; 
        BC_bot_var = 'Flux'; 
        BC_bot_val = 'Qd';
        BC_bot_type = 'Inflow';

        IC = ["T(1,:)=T0;"];

        % physic
        T0 = 0;     
        Td = 0;     
        A  = 1.25*10^-6;  
        k  = 2.5;         
        Qd = 42*10^-3;    
        rho = 2.3*10^3; 
        Cp = 10^3;

        % output
        outtime      = 1*10^6*365*24*3600; 
        expcase      = simcase;
        if_animation = 1; 
        ani_name = 'Simulation_2.gif';
        
    end
    if (simcase==3)
        % domain
        zd = 10*10^0;             
        z0 = 0;                   
        dz = 0.02;                
        dt = 60;                  
        tottime = 365*24*3600;    

        BC_top = 'Dirichlet';
        BC_top_var = 'Temperature';
        BC_top_val = 0;
        BC_top_type = 'Outflow';
        %
        BC_bot = 'Dirichlet'; 
        BC_bot_var = 'Temperature'; 
        BC_bot_val = 0;
        BC_bot_type = 'Outflow';

        IC = ["T(1,1:200)=0;", "T(1,201:301)=T0;", "T(1,302:501)=0;"];

        % physic
        T0 = 1000;       
        Td = 0;           
        A  = 1.25*10^-6;  
        k  = 1.7;         
        Qd = 21*10^-3;    
        rho = 2.3*10^3;    
        Cp = 10^3;       

        % output
        outtime      = 24*3600;
        expcase      = simcase;
        if_animation = 1;
        ani_name = 'Simulation_3.gif';
    end
    
return

