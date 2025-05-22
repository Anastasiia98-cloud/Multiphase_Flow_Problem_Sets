%% Hyperbolic flow solver for computation of new state vector at time step n+1

function [W] = solveEulerHLLC(W,F,p,c,Nx,hx,dt)

    W_buf = zeros(3,Nx-1);

    for k=1:3
        
        for i=1:Nx-1

            if(i==1)
                F_L = F(k,i);
                [F_R] = flux_HLLC_F(F,W,p,c,k,i,1);

            elseif(i==Nx-1)
                [F_L] = flux_HLLC_F(F,W,p,c,k,i,-1);
                F_R = F(k,i);

            else
                [F_L] = flux_HLLC_F(F,W,p,c,k,i,-1);
                [F_R] = flux_HLLC_F(F,W,p,c,k,i,1);
            end

            %% COMPUTE STATE VECTOR AT NEXT TIME STEP IN W_buf
            % Your Code Here 
            W_buf(k,i) = W(k,i) - (dt/hx)*(F_R-F_L);
            
        end
    end

    W = W_buf;
end

function [flux] = flux_HLLC_F(F,W,p,c,ind,Row,sign)

%% Definition of states, fluxes, and natural variables on both sides of cell boundaries
    if (sign>0) % F_R from solveEulerHLLC
		rho_L = W(1,Row); 
		u_L = W(2,Row)/rho_L;
		p_L = p(Row,1);
        c_L = c(Row,1);
        
        rho_R = W(1,Row+sign);
		u_R = W(2,Row+sign)/rho_R;
		p_R = p(Row+sign,1);
        c_R = c(Row+sign,1);

		W_L = W(ind,Row); % state at i
		W_R = W(ind,Row+sign); % state at i+1

		F_L = F(ind,Row); % flux at i
		F_R = F(ind,Row+sign); % flux at i+1
    else % F_L from solveEulerHLLC
        rho_L = W(1,Row+sign);
		u_L = W(2,Row+sign)/rho_L;
		p_L = p(Row+sign,1);
        c_L = c(Row+sign,1);
        
        rho_R = W(1,Row);
		u_R = W(2,Row)/rho_R;
		p_R = p(Row,1);
        c_R = c(Row,1);

		W_L = W(ind,Row+sign); % state at i-1
		W_R = W(ind,Row); % state at i

		F_L = F(ind,Row+sign); % flux at i-1
		F_R = F(ind,Row); % flux at i
    end

%% Calculation of wave speeds S_l, S_R, S_star

    S_L = min(u_L-c_L,u_R-c_R); % Toro (10.38)
    S_R = max(u_L+c_L,u_R+c_R);

    S_star = ( p_R - p_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R) )/( rho_L*(S_L-u_L) - rho_R*(S_R-u_R) ); % Toro (10.58)
  
%% Calculation of states W_star_L and W_star_R at the cell boundaries (Toro (10.33))

    if(ind==1) % Mass
		W_star_L = rho_L*(S_L-u_L)/(S_L-S_star);
		W_star_R = rho_R*(S_R-u_R)/(S_R-S_star);
    elseif(ind==2) % Momentum
		W_star_L = rho_L*(S_L-u_L)/(S_L-S_star)*S_star;
		W_star_R = rho_R*(S_R-u_R)/(S_R-S_star)*S_star;
    else  % Energy
		W_star_L = rho_L*(S_L-u_L)/(S_L-S_star)*(W_L/rho_L + (S_star-u_L)*(S_star+p_L/rho_L/(S_L-u_L)));
		W_star_R = rho_R*(S_R-u_R)/(S_R-S_star)*(W_R/rho_R + (S_star-u_R)*(S_star+p_R/rho_R/(S_R-u_R)));
    end

%% Calculation of HLLC flux (Toro (10.34))   
    
    F_star_L = F_L + S_L*(W_star_L-W_L);
    F_star_R = F_R + S_R*(W_star_R-W_R);

    if(S_L >= 0.0)					
        flux = F_L;
    elseif( (S_L <= 0.0) && (S_star >= 0.0) ) 
        flux = F_star_L;
    elseif( (S_star <= 0.0) && (S_R >= 0.0) ) 	
        flux = F_star_R;
    elseif(S_R <= 0.0) 				
        flux = F_R;
    end
   
end