%***************************************************************************************************
%*   Solve 2D laplace equation by presented code.
%*   I take no responsibilities for any errors in the code or damage thereby.
%*   Please notify me at zolfaghari1992iut@gmail.com if the code is used in any type of application.
%***************************************************************************************************
%*   Developer   : Ali Zolfaghari Sichani (14-06-2018)
%***************************************************************************************************
%*   References  : 
%*   Computational Fluid Mechanics and Heat Transfer.
%*   by John C. Tannehill (Author), Dale Anderson (Author), Richard H. Pletcher (Author).
%***************************************************************************************************
%*   Laplace Equation in two-dimensional square domain. (solving by centeral difference scheme)   :   
%*   Uxx + Uyy = 0
%*   Inputs      :
%*   MType = 1 : GAUSS SEIDEL
%*   MType = 2 : SUCCESSIVE OVER RELAXATION
%*   MType = 3 : STRONGLY IMPLICIT PROCEDURE
%*   MType = 4 : MULTIGRID + GAUSS SEIDEL
%*   MType = 5 : MULTIGRID + STRONGLY IMPLICIT PROCEDURE
%*   M          (number of division of domain in x-direction   )
%*   N          (number of division of domain in y-direction   )
%*   Ws         (relaxation factor                             )
%*   MGITER     (max. allowable itrerations of gauss           )
%*   MAXERROR   (max. allowable error                          )
%*   MAXITER    (max. allowable itrerations                    )
%*   Outputs     :
%*   plot numerical and exact solution
%***************************************************************************************************


clear,clc,close all
format compact
format long


% Inputs
MType = 3;
M = 32;
N = 32;
MGITER = 3;
MAXITER = 200;
MAXERROR = 0.00001;
Ws = 0.8;


% Definition
M2 = M/2;
N2 = N/2;
M4 = M/4;
N4 = N/4;
M8 = M/8;
N8 = N/8;

dx = 1.0/M;
dy = 1.0/N;
dx2 = 1.0/M2;
dy2 = 1.0/N2;
dx4 = 1.0/M4;
dy4 = 1.0/N4;
dx8 = 1.0/M8;
dy8 = 1.0/N8;

x = zeros(1,(M+1)*(N+1));
y = zeros(1,(M+1)*(N+1));
B = zeros(1,(M+1)*(N+1));
UE = zeros(1,(M+1)*(N+1));
UI = zeros(1,(M+1)*(N+1));
UII = zeros(1,(M+1)*(N+1));

% Analytical
for i = 1:M+1
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        x(i) = (i-1)*dx;
        y(j) = (j-1)*dy;
        UE(ip) = (sin(pi*x(i))*sinh(pi*y(j)))/sinh(pi);
    end
end

% B.C. & I.C.
j = N+1;
for i = 2:M
    ip = (j-1)*(M+1)+i;
    B(ip) = sin(pi*x(i));
    UI(ip) = sin(pi*x(i));
    UII(ip) = sin(pi*x(i));
end

% Solver
Iteration = 1;
switch(MType)
    
    case 1
        % Numerical Method #1
        [ Ap,An,As,Ae,Aw ] = MatrixCalc( dx,dy,M,N );
        [ UII,Residual,Iteration ] = GSSolver( UII,UI,B,Ap,An,As,Ae,Aw,M,N,1.0,MAXITER,MAXERROR,1 );
        Error = sqrt(sum((UE-UII).^2))/((M+1)*(N+1));
        
    case 2
        % Numerical Method #2
        [ Ap,An,As,Ae,Aw ] = MatrixCalc( dx,dy,M,N );
        [ UII,Residual,Iteration ] = GSSolver( UII,UI,B,Ap,An,As,Ae,Aw,M,N,Ws,MAXITER,MAXERROR,1 );
        Error = sqrt(sum((UE-UII).^2))/((M+1)*(N+1));
        
    case 3
        % Numerical Method #3
        [ Ap,An,As,Ae,Aw ] = MatrixCalc( dx,dy,M,N );
        [ UII,Residual,Iteration ] = SIPSolver( UII,UI,B,Ap,An,As,Ae,Aw,M,N,MAXITER,MAXERROR,1 );
        Error = sqrt(sum((UE-UII).^2))/((M+1)*(N+1));
        
    case 4
        % Numerical Method #4
        
        [ Ap,An,As,Ae,Aw ] = MatrixCalc( dx,dy,M,N );
        [ Ap2,An2,As2,Ae2,Aw2 ] = MatrixCalc( dx2,dy2,M2,N2 );
        [ Ap4,An4,As4,Ae4,Aw4 ] = MatrixCalc( dx4,dy4,M4,N4 );
        [ Ap8,An8,As8,Ae8,Aw8 ] = MatrixCalc( dx8,dy8,M8,N8 );
        while( Iteration <= MAXITER )
            
            U2h = zeros(1,(M2+1)*(N2+1));
            U4h = zeros(1,(M4+1)*(N4+1));
            U8h = zeros(1,(M8+1)*(N8+1));
            U2g = zeros(1,(M2+1)*(N2+1));
            U4g = zeros(1,(M4+1)*(N4+1));
            U8g = zeros(1,(M8+1)*(N8+1));
            
            [ Uh,Res,Iter ] = GSSolver( UII,UI,B,Ap,An,As,Ae,Aw,M,N,1.0,MGITER,MAXERROR,1 );
            
            [ Rh ] = ResidualCalc( Uh,B,Ap,An,As,Ae,Aw,M,N );
            
            [ R2h ] = RestrictCalc( Rh,M,N );
            
            [ U2h,Res,Iter ] = GSSolver( U2h,U2g,R2h,Ap2,An2,As2,Ae2,Aw2,M2,N2,1.0,MGITER,MAXERROR,1 );
            
            [ R2h ] = ResidualCalc( U2h,R2h,Ap2,An2,As2,Ae2,Aw2,M2,N2 );
            
            [ R4h ] = RestrictCalc( R2h,M2,N2 );
            
            [ U4h,Res,Iter ] = GSSolver( U4h,U4g,R4h,Ap4,An4,As4,Ae4,Aw4,M4,N4,1.0,MGITER,MAXERROR,1 );
            
            [ R4h ] = ResidualCalc( U4h,R4h,Ap4,An4,As4,Ae4,Aw4,M4,N4 );
            
            [ R8h ] = RestrictCalc( R4h,M4,N4 );
            
            [ U8h,Res,Iter ] = GSSolver( U8h,U8g,R8h,Ap8,An8,As8,Ae8,Aw8,M8,N8,1.0,MAXITER,MAXERROR,2 );
            
            [ I4h ] = ProlongateCalc( U8h,M8,N8 );
            
            U4h = U4h+I4h;
            
            [ I2h ] = ProlongateCalc( U4h,M4,N4 );
            
            U2h = U2h+I2h;
            
            [ Ih ] = ProlongateCalc( U2h,M2,N2 );
            
            UII = Uh+Ih;
            
            Residual = sqrt(sum((UII-UI).^2))/((M+1)*(N+1));
            UI = UII;
            Iteration = Iteration+1;
        end
        Error = sqrt(sum((UE-UII).^2))/((M+1)*(N+1));
        
    case 5
        % Numerical Method #5
        
        [ Ap,An,As,Ae,Aw ] = MatrixCalc( dx,dy,M,N );
        [ Ap2,An2,As2,Ae2,Aw2 ] = MatrixCalc( dx2,dy2,M2,N2 );
        [ Ap4,An4,As4,Ae4,Aw4 ] = MatrixCalc( dx4,dy4,M4,N4 );
        [ Ap8,An8,As8,Ae8,Aw8 ] = MatrixCalc( dx8,dy8,M8,N8 );
        while( Iteration <= MAXITER )
            
            U2h = zeros(1,(M2+1)*(N2+1));
            U4h = zeros(1,(M4+1)*(N4+1));
            U8h = zeros(1,(M8+1)*(N8+1));
            U2g = zeros(1,(M2+1)*(N2+1));
            U4g = zeros(1,(M4+1)*(N4+1));
            U8g = zeros(1,(M8+1)*(N8+1));
            
            [ Uh,Res,Iter ] = SIPSolver( UII,UI,B,Ap,An,As,Ae,Aw,M,N,MGITER,MAXERROR,1 );
            
            [ Rh ] = ResidualCalc( Uh,B,Ap,An,As,Ae,Aw,M,N );
            
            [ R2h ] = RestrictCalc( Rh,M,N );
            
            [ U2h,Res,Iter ] = SIPSolver( U2h,U2g,R2h,Ap2,An2,As2,Ae2,Aw2,M2,N2,MGITER,MAXERROR,1 );
            
            [ R2h ] = ResidualCalc( U2h,R2h,Ap2,An2,As2,Ae2,Aw2,M2,N2 );
            
            [ R4h ] = RestrictCalc( R2h,M2,N2 );
            
            [ U4h,Res,Iter ] = SIPSolver( U4h,U4g,R4h,Ap4,An4,As4,Ae4,Aw4,M4,N4,MGITER,MAXERROR,1 );
            
            [ R4h ] = ResidualCalc( U4h,R4h,Ap4,An4,As4,Ae4,Aw4,M4,N4 );
            
            [ R8h ] = RestrictCalc( R4h,M4,N4 );
            
            [ U8h,Res,Iter ] = SIPSolver( U8h,U8g,R8h,Ap8,An8,As8,Ae8,Aw8,M8,N8,MAXITER,MAXERROR,2 );
            
            [ I4h ] = ProlongateCalc( U8h,M8,N8 );
            
            U4h = U4h+I4h;
            
            [ I2h ] = ProlongateCalc( U4h,M4,N4 );
            
            U2h = U2h+I2h;
            
            [ Ih ] = ProlongateCalc( U2h,M2,N2 );
            
            UII = Uh+Ih;
            
            Residual = sqrt(sum((UII-UI).^2))/((M+1)*(N+1));
            UI = UII;
            Iteration = Iteration+1;
        end
        Error = sqrt(sum((UE-UII).^2))/((M+1)*(N+1));
end

U = zeros(M+1,N+1);
UC = zeros(M+1,N+1);
for i = 1:M+1
    for j = 1:N+1
        ip = (j-1)*(M+1)+i;
        U(i,j) = UII(ip);
        UC(i,j) = UE(ip);
    end
end
% plot
figure(1);
contourf(U');title('Numerical');
colorbar('location','EastOutside')
figure(2);
contourf(UC');title('Analytical');
colorbar('location','EastOutside')
