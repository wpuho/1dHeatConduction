function T = Liao_BC(BC,BC_var,BC_val,Ts,dz,k,Tjm1,Tj,c)

    if (strcmp(BC,'Dirichlet')==1) && (strcmp(BC_var,'Temperature')==1)
        T = BC_val;
    end
    
    if (strcmp(BC,'Neumann')==1) && (strcmp(BC_var,'Flux')==1)
        T = 1*dz*BC_val/k + Tjm1;
    end
    
    if (strcmp(BC,'Mixed')==1) && (strcmp(BC_var,'Mix')==1)
        T = 1*dz*BC_val/k + Ts;
    end
 
    if (strcmp(BC,'Close')==1)
        T = Tjm1;
        %T = s*Tjm1 + (1-s)*Tj + c;
        %T = Ts;
    end
    
end