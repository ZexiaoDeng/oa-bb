function Pi = simplicial_subdivision( w, lambda, P )
%   SIMPLICIAL_SUBDIVISION  计算单纯形的辐射状剖分
%
%       P = Union Pi(w) i in I
%       w in int( P )
%       w = sum( lambda_i*vi ), lambda_i >= 0, i = 1, ..., n
%       sum( lambda_i ) = 1, i = 1, ..., n
%
%   输入:
%       w      : 单纯形 P 的内点
%       lambda : 一个有限指标集, int( Pi )非空
%       P      : 原单纯形
%
%   输出:
%       Pi     : 子单纯形 cell
%                Pi = co( [ v0, v1, ..., vi-1, w, vi+1, ..., vn ] )
%
%    see also 
%       
%       全局优化引论, R. Horst, P.M. Pardalos, N.V. Thoai 著, 清华大学出版社, P148
%

Pi    = {} ;
icout = 1 ;

for idx = 1: size( P, 2 )
    if lambda( idx ) > 0
        Pi{ icout } = P ;
        Pi{ icout }( : , idx ) = w ;
        icout = icout + 1 ;
    else
        continue ;
    end
end

return ;

end




