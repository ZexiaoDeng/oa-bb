% ============================================
% 全局优化引论
% R. Horst, P.M. Pardalos, N.V. Thoai 著
% 黄红选 译
% 梁治安 校
% P178
% =============================================
%   min f( x )
%       x in D
%
%   D = { x: A*x <= b, x >= 0 } in R+^3
%
%   f(x) = -( abs( x1 + x2/2 + 2*x3/3) )^( 3/2 ) - x1^2
%   A = [  1,  1,  1   ; ...
%          1,  1, -1/4 ; ...
%         -2, -2,  1   ; ...
%          0,  0,  1   ; ] ;
%   b = [ 2 ; 1 ; 1 ; 3 ] ;
% 

clc ;
clear ;
close all ;

format long

path = './bt-1.3' ;
addpath( path ) ;

Aineq = [  1,  1,  1   ; ...
           1,  1, -1/4 ; ...
          -2, -2,  1   ; ...
           0,  0,  1   ; ] ;
bineq = [ 2 ; ...
          1 ; ...
          1 ; ...
          3 ; ] ;

rep.M = eye( 3 ) ;
rep.B = Aineq ;
rep.b = bineq ;
rep.l = zeros( 3, 1 ) ;
rep.u = [ inf ; inf ; inf ; ] ;

P  = eval( polyh( rep ) ) ;      % 求出多胞体 P 的 H-rep, V-rep, P-rep
CH = vrep( P ) ;                 % 获取多胞体 P 的 V-rep
A  = adj( P ) ;                  % 获取顶点对应的链表( 邻接表表示形式 )

plot( P ) ;
grid on ;
hold on ;

for i = 1: size( CH.V', 1 )
    fprintf( '%8.4f\t%8.4f\t%8.4f', ...
        CH.V( 1, i ), ...
        CH.V( 2, i ), ...
        CH.V( 3, i ) ) ; % 顶点
    fprintf( '%8d\t%8d', A{i}(1), A{i}(2) ) ;               % 顶点对应的链表
    fprintf( '\n' ) ;
end

% 包含可行集合的第一个单纯形
V1 = [ 0 ; 0 ; 0 ; ] ;
V2 = [ 4 ; 0 ; 0 ; ] ;
V3 = [ 0 ; 4 ; 0 ; ] ;
V4 = [ 0 ; 0 ; 4 ; ] ;

S1 = [ V1, V2, V3, V4 ; ] ;

rep.V = S1 ;
rep.D = zeros( 3 ) ;
P1    = polyh( rep, 'v' ) ;

opt.color = [0.5 0.6 0.5] ;
plot( P1, opt ) ;

% =================
% 式 ( 3.104 )
% =================
c1 = [ oracle( V1 ) ; ...
       oracle( V2 ) ; ...
       oracle( V3 ) ; ...
       oracle( V4 ) ; ] ;

Aineq1 = Aineq*S1   ;
bineq1 = bineq ;
Aeq1   = ones( 1, 4 ) ;
beq1   = 1 ;
lb1    = zeros( 4, 1 ) ;
ub1    = inf*ones( 4, 1 ) ;

ops = optimoptions( 'linprog'  , ...
                    'Algorithm', 'dual-simplex', ...
                    'display'  , 'none' ) ;
                
[ alpha   , ...
  fval    , ...
  exitflag ] = linprog( c1    ,        ...
                        Aineq1, bineq1, ...
                        Aeq1  , beq1  , ...
                        lb1   , ub1   , ops )
% =================
% 式 ( 3.105 )
% =================
x_opt = S1*alpha
f_opt = oracle( x_opt )

% 更新当前最优解
gamma1 = f_opt ;

% ======================================================
% 单纯形的辐射状分割
% 求多面体 D 与多面体 S1 顶点距离最大的点进行辐射状分割
% 其中多面体 D 中顶点的选择是随机的, 但不能 S1 的顶点重合
% ======================================================
w      = x_opt ;          % w = x_opt
lambda = alpha ;
S1p = simplicial_subdivision( w, alpha, S1 ) ;

% =================
% S11
% =================
c11 = [ oracle( S1p{ 1 }( : , 1 ) ) ; ...
        oracle( S1p{ 1 }( : , 2 ) ) ; ...
        oracle( S1p{ 1 }( : , 3 ) ) ; ...
        oracle( S1p{ 1 }( : , 4 ) ) ; ] ;

Aineq11 = Aineq*S1p{ 1 }   ;
bineq11 = bineq ;
Aeq11   = ones( 1, 4 ) ;
beq11   = 1 ;
lb11    = zeros( 4, 1 ) ;
ub11    = inf*ones( 4, 1 ) ;

ops = optimoptions( 'linprog'  , ...
                    'Algorithm', 'dual-simplex', ...
                    'display'  , 'none' ) ;
                
[ alpha   , ...
  fval    , ...
  exitflag ] = linprog( c11    ,        ...
                        Aineq11, bineq11, ...
                        Aeq11  , beq11  , ...
                        lb11   , ub11   , ops )

% =================
% S12
% =================
c12 = [ oracle( S1p{ 2 }( : , 1 ) ) ; ...
        oracle( S1p{ 2 }( : , 2 ) ) ; ...
        oracle( S1p{ 2 }( : , 3 ) ) ; ...
        oracle( S1p{ 2 }( : , 4 ) ) ; ] ;

Aineq12 = Aineq*S1p{ 2 }   ;
bineq12 = bineq ;
Aeq12   = ones( 1, 4 ) ;
beq12   = 1 ;
lb12    = zeros( 4, 1 ) ;
ub12    = inf*ones( 4, 1 ) ;

ops = optimoptions( 'linprog'  , ...
                    'Algorithm', 'dual-simplex', ...
                    'display'  , 'none' ) ;
                
[ alpha   , ...
  fval    , ...
  exitflag ] = linprog( c12    ,        ...
                        Aineq12, bineq12, ...
                        Aeq12  , beq12  , ...
                        lb12   , ub12   , ops )

% =================
% S13
% =================
c13 = [ oracle( S1p{ 3 }( : , 1 ) ) ; ...
        oracle( S1p{ 3 }( : , 2 ) ) ; ...
        oracle( S1p{ 3 }( : , 3 ) ) ; ...
        oracle( S1p{ 3 }( : , 4 ) ) ; ] ;

Aineq13 = Aineq*S1p{ 3 }   ;
bineq13 = bineq ;
Aeq13   = ones( 1, 4 ) ;
beq13   = 1 ;
lb13    = zeros( 4, 1 ) ;
ub13    = inf*ones( 4, 1 ) ;

ops = optimoptions( 'linprog'  , ...
                    'Algorithm', 'dual-simplex', ...
                    'display'  , 'none' ) ;
                
[ alpha   , ...
  fval    , ...
  exitflag ] = linprog( c13    ,        ...
                        Aineq13, bineq13, ...
                        Aeq13  , beq13  , ...
                        lb13   , ub13   , ops )
% =================
% 式 ( 3.105 )
% =================
x_opt = S1p{ 3 }*alpha
f_opt = oracle( x_opt )

% 更新当前最优解
gamma1 = f_opt ;


close all
plot( P ) ;
hold on
S2 = S1p{ 3 } ;
rep.V = S2 ;
P2    = polyh( rep, 'v' ) ;
opt.color = [0.5 0.2 0.1 ] ;
plot( P2, opt ) ;

% ======================================================
% 单纯形的辐射状分割
% ======================================================
w      = x_opt ;          % w = x_opt
lambda = alpha ;
S2p = simplicial_subdivision( w, alpha, S2 ) ;

% =================
% S21
% =================
c21 = [ oracle( S2p{ 1 }( : , 1 ) ) ; ...
        oracle( S2p{ 1 }( : , 2 ) ) ; ...
        oracle( S2p{ 1 }( : , 3 ) ) ; ...
        oracle( S2p{ 1 }( : , 4 ) ) ; ] ;

Aineq21 = Aineq*S2p{ 1 }   ;
bineq21 = bineq ;
Aeq21   = ones( 1, 4 ) ;
beq21   = 1 ;
lb21    = zeros( 4, 1 ) ;
ub21    = inf*ones( 4, 1 ) ;

ops = optimoptions( 'linprog'  , ...
                    'Algorithm', 'dual-simplex', ...
                    'display'  , 'none' ) ;

[ alpha   , ...
  fval    , ...
  exitflag ] = linprog( c21    ,        ...
                        Aineq21, bineq21, ...
                        Aeq21  , beq21  , ...
                        lb21   , ub21   , ops )

% =================
% S22
% =================
c22 = [ oracle( S2p{ 2 }( : , 1 ) ) ; ...
        oracle( S2p{ 2 }( : , 2 ) ) ; ...
        oracle( S2p{ 2 }( : , 3 ) ) ; ...
        oracle( S2p{ 2 }( : , 4 ) ) ; ] ;

Aineq22 = Aineq*S2p{ 2 }   ;
bineq22 = bineq ;
Aeq22   = ones( 1, 4 ) ;
beq22   = 1 ;
lb22    = zeros( 4, 1 ) ;
ub22    = inf*ones( 4, 1 ) ;

ops = optimoptions( 'linprog'  , ...
                    'Algorithm', 'dual-simplex', ...
                    'display'  , 'none' ) ;

[ alpha   , ...
  fval    , ...
  exitflag ] = linprog( c22    ,        ...
                        Aineq22, bineq22, ...
                        Aeq22  , beq22  , ...
                        lb22   , ub22   , ops )

% =================
% 式 ( 3.105 )
% =================
x_opt = S2p{ 1 }*alpha
f_opt = oracle( x_opt )

% 更新当前最优解
gamma2 = f_opt

% 目标函数信息 oracle
function f = oracle( x )
    f = -( abs( x(1) + x(2)/2 + 2*x(3)/3) )^( 3/2 ) - x(1)^2 ;
end









