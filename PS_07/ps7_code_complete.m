%========================================================================
% A very simple Navier-Stokes solver for a drop falling in a rectangular
% box. A forward in time, centered in space discretization is used.
% The density is advected by a front tracking scheme and a stretched grid
% is used, allowing us to concentrate the grid points in specific areas
%========================================================================
%domain size and physical variables
Lx=1.0;Ly=1.0; % Domain size
gx=0.0;gy=100.0; % Gravity
rho1=1.0; rho2=2.0; % Density of surrounding liquid and drop
m0=0.01; % Viscosity
sigma=10; % Surface tension
rro=rho1;
un=0;us=0;ve=0;vw=0;time=0.0; % Velocity boundary conditions
rad=0.15;xc=0.5;yc=0.7; % Initial drop size and location

% Numerical variables
nx=32;ny=32; % Numbers of cells
dt=0.0004;nstep=2000; % Time step size & number of time steps
maxit=200;maxError=0.001;beta=1.2; % Iteration parameters

% Zero various arrys
u=zeros(nx+1,ny+2);  v=zeros(nx+2,ny+1);  p=zeros(nx+2,ny+2); % Velocity / pressure arrays
ut=zeros(nx+1,ny+2); vt=zeros(nx+2,ny+1); % Temporary velocity arrays
tmp1=zeros(nx+2,ny+2); tmp2=zeros(nx+2,ny+2); % Arrays for velocity terms in pressure equation
uu=zeros(nx+1,ny+1); vv=zeros(nx+1,ny+1); % Arrays for plotting
fx=zeros(nx+2,ny+2); fy=zeros(nx+2,ny+2); % Surface tension forces

% Set the grid
dx=Lx/nx;dy=Ly/ny; % Cell size
for i=1:nx+2
    x(i)=dx*(i-1.5); % x coordinates
end
for j=1:ny+2
    y(j)=dy*(j-1.5); % y coordinates
end

% Set density in the domain and the drop
r=zeros(nx+2,ny+2)+rho1; % Set density of surrounding liquid
for i=2:nx+1
    for j=2:ny+1
        if ( (x(i)-xc)^2+(y(j)-yc)^2 < rad^2)
            r(i,j)=rho2; % Set density of drop
        end
    end
end

%================== SETUP THE FRONT ===================
Nf=100; % Number of front points
xf=zeros(1,Nf+2);yf=zeros(1,Nf+2); % Vectors of front point coordinates
uf=zeros(1,Nf+2);vf=zeros(1,Nf+2); % Vectors of front point velocities
tx=zeros(1,Nf+2); ty=zeros(1,Nf+2); % Surface tangent vectors

for l=1:Nf+2
    xf(l)=xc-rad*sin(2.0*pi*(l-1)/(Nf)); % Set front point x-coordinates
    yf(l)=yc+rad*cos(2.0*pi*(l-1)/(Nf)); % Set front point y-coordinates
end

%================== START TIME LOOP======================================
for is=1:nstep,is
    
    %------------------ FIND SURFACE TENSION --------------
    fx=zeros(nx+2,ny+2);fy=zeros(nx+2,ny+2);  % Set fx & fy to zero
    for l=1:Nf+1
        ds=sqrt((xf(l+1)-xf(l))^2+(yf(l+1)-yf(l))^2);
        tx(l)=(xf(l+1)-xf(l))/ds; % Tangent x-vectors
        ty(l)=(yf(l+1)-yf(l))/ds; % Tangent y-vectors
    end
    tx(Nf+2)=tx(2);ty(Nf+2)=ty(2);
    
	for l=2:Nf+1     % Distribute to the fixed grid
	    nfx=sigma*(tx(l)-tx(l-1));nfy=sigma*(ty(l)-ty(l-1));
        
        ip=floor(xf(l)/dx)+1; jp=floor((yf(l)+0.5*dy)/dy)+1;
        ax=xf(l)/dx-ip+1; ay=(yf(l)+0.5*dy)/dy-jp+1;
        fx(ip,jp)    =fx(ip,jp)+(1.0-ax)*(1.0-ay)*nfx/dx/dy;
        fx(ip+1,jp)  =fx(ip+1,jp)+ax*(1.0-ay)*nfx/dx/dy;
        fx(ip,jp+1)  =fx(ip,jp+1)+(1.0-ax)*ay*nfx/dx/dy;
        fx(ip+1,jp+1)=fx(ip+1,jp+1)+ax*ay*nfx/dx/dy;

        ip=floor((xf(l)+0.5*dx)/dx)+1; jp=floor(yf(l)/dy)+1;
        ax=(xf(l)+0.5*dx)/dx-ip+1; ay=yf(l)/dy-jp+1;	  
        fy(ip,jp)    =fy(ip,jp)+(1.0-ax)*(1.0-ay)*nfy/dx/dy;
        fy(ip+1,jp)  =fy(ip+1,jp)+ax*(1.0-ay)*nfy/dx/dy;
        fy(ip,jp+1)  =fy(ip,jp+1)+(1.0-ax)*ay*nfy/dx/dy;
        fy(ip+1,jp+1)=fy(ip+1,jp+1)+ax*ay*nfy/dx/dy;  
    end
    
    %---------------------------------------------------------
    % tangential velocity at boundaries
    u(1:nx+1,1)=2*us-u(1:nx+1,2);u(1:nx+1,ny+2)=2*un-u(1:nx+1,ny+1);
    v(1,1:ny+1)=2*vw-v(2,1:ny+1);v(nx+2,1:ny+1)=2*ve-v(nx+1,1:ny+1);
    
    for i=2:nx,for j=2:ny+1      % TEMPORARY u-velocity (accounting for advection, diffusion, and gravity terms)
            ut(i,j)=u(i,j)+dt*(-0.25*(((u(i+1,j)+u(i,j))^2-(u(i,j)+   ...
                u(i-1,j))^2)/dx+((u(i,j+1)+u(i,j))*(v(i+1,j)+         ...
                v(i,j))-(u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))/dy)+ ...
                (m0/(0.5*(r(i+1,j)+r(i,j))) )*(                          ...
                (u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+            ...
                (u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2 )+        ...
                fx(i,j)/(0.5*(r(i+1,j)+r(i,j)))-          ...
                (1.0 -rro/(0.5*(r(i+1,j)+r(i,j))) )*gx      );
        end,end
    
    for i=2:nx+1,for j=2:ny       % TEMPORARY v-velocity (accounting for advection, diffusion, and gravity terms)
            vt(i,j)=v(i,j)+dt*(-0.25*(((u(i,j+1)+u(i,j))*(v(i+1,j)+   ...
                v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j)))/dx+  ...
                ((v(i,j+1)+v(i,j))^2-(v(i,j)+v(i,j-1))^2)/dy)+        ...
                (m0/(0.5*(r(i,j+1)+r(i,j))) )*(                          ...
                (v(i+1,j)-2*v(i,j)+v(i-1,j))/dx^2+            ...
                (v(i,j+1)-2*v(i,j)+v(i,j-1))/dy^2 )+          ...
                fy(i,j)/(0.5*(r(i,j+1)+r(i,j)))-           ...
                (1.0 -rro/(0.5*(r(i,j+1)+r(i,j))) )*gy      );
        end,end
    %========================================================================
    %  Compute source term and the coefficient for p(i,j)
    rt=r; lrg=1000;
    rt(1:nx+2,1)=lrg;rt(1:nx+2,ny+2)=lrg;
    rt(1,1:ny+2)=lrg;rt(nx+2,1:ny+2)=lrg;
    
    for i=2:nx+1
        for j=2:ny+1
            tmp1(i,j)= (0.5/dt)*( (ut(i,j)-ut(i-1,j))/dx+(vt(i,j)-vt(i,j-1))/dy ); % Velocity term in pressure equation
            tmp2(i,j)=1.0/( (1./dx)*( 1./(dx*(rt(i+1,j)+rt(i,j)))+...  % Density term in pressure equation
                1./(dx*(rt(i-1,j)+rt(i,j)))  )+...
                (1./dy)*(1./(dy*(rt(i,j+1)+rt(i,j)))+...
                1./(dy*(rt(i,j-1)+rt(i,j)))   )   );
        end
    end
    
    for it=1:maxit	               % SOLVE FOR PRESSURE (iteratively; to ensure that pressure satisfies mass conservation equation in cell)
        for i=2:nx+1
            for j=2:ny+1
                oldArray=p;
                p(i,j)=(1.0-beta)*p(i,j)+beta* tmp2(i,j)*(...
                    (1./dx)*( p(i+1,j)/(dx*(rt(i+1,j)+rt(i,j)))+...
                    p(i-1,j)/(dx*(rt(i-1,j)+rt(i,j))) )+...
                    (1./dy)*( p(i,j+1)/(dy*(rt(i,j+1)+rt(i,j)))+...
                    p(i,j-1)/(dy*(rt(i,j-1)+rt(i,j))) ) - tmp1(i,j));
            end
        end
        if max(max(abs(oldArray-p))) <maxError
            break
        end
    end
    
    for i=2:nx
        for j=2:ny+1   % CORRECT THE u-velocity (accounting for the new pressure)
            u(i,j)=ut(i,j)-dt*(2.0/dx)*(p(i+1,j)-p(i,j))/(r(i+1,j)+r(i,j));
        end
    end
    
    for i=2:nx+1
        for j=2:ny   % CORRECT THE v-velocity (accounting for the new pressure)
            v(i,j)=vt(i,j)-dt*(2.0/dy)*(p(i,j+1)-p(i,j))/(r(i,j+1)+r(i,j));
        end
    end
    %================== ADVECT FRONT =====================
    for l=2:Nf+1
        ip=floor(xf(l)/dx)+1; jp=floor((yf(l)+0.5*dy)/dy)+1; % Nearest grid point indices
        ax=xf(l)/dx-ip+1;ay=(yf(l)+0.5*dy)/dy-jp+1;	 % Normalized distance to nearest grid points
        uf(l)=(1.0-ax)*(1.0-ay)*u(ip,jp)+ax*(1.0-ay)*u(ip+1,jp)+... % Front x-velocity (Eq. (3.5))
            (1.0-ax)*ay*u(ip,jp+1)+ax*ay*u(ip+1,jp+1);
        
        ip=floor((xf(l)+0.5*dx)/dx)+1; jp=floor(yf(l)/dy)+1;
        ax=(xf(l)+0.5*dx)/dx-ip+1;ay=yf(l)/dy-jp+1;
        vf(l)=(1.0-ax)*(1.0-ay)*v(ip,jp)+ax*(1.0-ay)*v(ip+1,jp)+... % Front y-velocity
            (1.0-ax)*ay*v(ip,jp+1)+ax*ay*v(ip+1,jp+1);
    end
    
    for i=2:Nf+1
        xf(i)=xf(i)+dt*uf(i); %MOVE THE FRONT
        yf(i)=yf(i)+dt*vf(i);
    end
    xf(1)=xf(Nf+1);yf(1)=yf(Nf+1);xf(Nf+2)=xf(2);yf(Nf+2)=yf(2);
    
    %------------ Add points to the front ------------
    xfold=xf;yfold=yf; j=1; % old array of front coordinates
    for l=2:Nf+1
        ds=sqrt( ((xfold(l)-xf(j))/dx)^2 + ((yfold(l)-yf(j))/dy)^2); % Distance between front points
        if (ds > 0.5) % Distance too large
            j=j+1;xf(j)=0.5*(xfold(l)+xf(j-1));yf(j)=0.5*(yfold(l)+yf(j-1)); % Coordinates of new front point
            j=j+1;xf(j)=xfold(l);yf(j)=yfold(l); % Coordinates of next (old) front point
        elseif (ds < 0.25) % Distance too small
            % DO NOTHING!
        else % Distance in desired limits
            j=j+1;xf(j)=xfold(l);yf(j)=yfold(l); % Coordinates of next (old) front point
        end
    end
    Nf=j-1;
    xf(1)=xf(Nf+1);yf(1)=yf(Nf+1);xf(Nf+2)=xf(2);yf(Nf+2)=yf(2);	% Ghost front points
    
    %------------ distribute gradient --------------
    fx=zeros(nx+2,ny+2);fy=zeros(nx+2,ny+2);  % Set fx & fy to zero
    for l=2:Nf+1
        nfx=-0.5*(yf(l+1)-yf(l-1))*(rho2-rho1); % x normal vector (front density gradient)
        nfy=0.5*(xf(l+1)-xf(l-1))*(rho2-rho1);  % y Normal vector (front density gradient)
        
        ip=floor(xf(l)/dx)+1; jp=floor((yf(l)+0.5*dy)/dy)+1; % Nearest grid point indices
        ax=xf(l)/dx-ip+1; ay=(yf(l)+0.5*dy)/dy-jp+1;  % Normalized distance to nearest grid points
        fx(ip,jp)    =fx(ip,jp)+(1.0-ax)*(1.0-ay)*nfx/dx/dy; % Weighted density gradient at front
        fx(ip+1,jp)  =fx(ip+1,jp)+ax*(1.0-ay)*nfx/dx/dy; % Weighted density gradient at front
        fx(ip,jp+1)  =fx(ip,jp+1)+(1.0-ax)*ay*nfx/dx/dy; % Weighted density gradient at front
        fx(ip+1,jp+1)=fx(ip+1,jp+1)+ax*ay*nfx/dx/dy; % Weighted density gradient at front
        
        ip=floor((xf(l)+0.5*dx)/dx)+1; jp=floor(yf(l)/dy)+1;
        ax=(xf(l)+0.5*dx)/dx-ip+1; ay=yf(l)/dy-jp+1;
        fy(ip,jp)    =fy(ip,jp)+(1.0-ax)*(1.0-ay)*nfy/dx/dy;
        fy(ip+1,jp)  =fy(ip+1,jp)+ax*(1.0-ay)*nfy/dx/dy;
        fy(ip,jp+1)  =fy(ip,jp+1)+(1.0-ax)*ay*nfy/dx/dy;
        fy(ip+1,jp+1)=fy(ip+1,jp+1)+ax*ay*nfy/dx/dy;
    end
    
    %------------ construct the density --------------
    r=zeros(nx+2,ny+2)+rho1;
    for iter=1:maxit
        oldArray=r;
        for i=2:nx+1
            for j=2:ny+1
                r(i,j)=0.25*(r(i+1,j)+r(i-1,j)+r(i,j+1)+r(i,j-1)+...  % Eq. (3.20)
                    dx*fx(i-1,j)-dx*fx(i,j)+...
                    dy*fy(i,j-1)-dy*fy(i,j));
            end
        end
        if max(max(abs(oldArray-r))) <maxError
            break
        end
    end
    %========================================================================
    time=time+dt                   % plot the results
    uu(1:nx+1,1:ny+1)=0.5*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
    vv(1:nx+1,1:ny+1)=0.5*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
    for i=1:nx+1,xh(i)=dx*(i-1);end;     for j=1:ny+1,yh(j)=dy*(j-1);end
    hold off,contour(x,y,flipud(rot90(r))),axis equal,axis([0 Lx 0 Ly]);
    hold on;quiver(xh,yh,flipud(rot90(uu)),flipud(rot90(vv)),'r');
    plot(xf(1:Nf),yf(1:Nf),'k','linewidth',5);pause(0.01)
end
