function [CartesianR, CartesianV]=coe2rv(coe, mu, epsilon)
%coe2rv 根据经典轨道根数求地心惯性直角坐标系下的位置和速度分量
%输入参数 coe(6),mu:
% 	coe(1)半长轴a（米）:对于椭圆(含圆)和双曲线轨道,其值大于零;对于抛物线轨道,
%                      其值取为近星距.
% 	coe(2)偏心率e:对于圆轨道e=0;椭圆轨道0<e<1,抛物线轨道e=1,双曲线轨道e>1
% 	coe(3)轨道倾角i（弧度）:范围0<=i<=180度.
%	coe(4)升交点赤经Omega（弧度）:当轨道倾角为0时没有意义,可将其值取为0
%	coe(5)近拱点幅角omega（弧度）:当偏心率为0时没有意义,可将其值置为0
%	coe(6)真近点角f（弧度）:当偏心率为0时没有意义,可将其值取为omega+f,即纬度幅角
%   mu:中心引力场引力系数,万有引力常数与中心天体质量的乘积,单位m^3/s^2,
%      默认为地心引力系数3.98600441800e+14
%   epsilon:最小精度,默认为1.0e-12
%输出参数[CartesianR, CartesianV]:
%   位置列矢量CartesianR分量X、Y、Z和速度列矢量CartesianV的分量VX、VY、VZ,
%   单位：米、米每秒.
%注:输出参数为一个时,将位置和速度整合成一个6维列向量

% *             清华大学航天航空学院动力学与控制研究室博士生                 *
% *                蒋方华(jiangfh04@mails.thu.edu.cn)                      *
% *                 最近修改: 2009.2.24                                    *
%公式来源《航天器轨道理论》P.40，刘林著.

% Check inputs
if nargin < 3
  epsilon=1.0e-12;
  if nargin < 2
    mu=3.98600441800e+14;
    if nargin < 1
        error('ORBIT:coe2rv:NotEnoughInputs',...
              'Not enough input arguments.  See coe2rv.');
    end
  end
end
if(length(coe)~=6)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '经典轨道根数应为6个.  See coe2rv.');
end
if(mu<=0.0||epsilon<=eps)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '引力常数应为正,最小精度应当大于系统最小双精度值.  See coe2rv.');   
end
if(coe(1)<0.0||coe(2)<0.0)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '半长轴和偏心率都不能为负.  See coe2rv.');    
end
if(coe(3)<0.0||coe(3)>pi)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '轨道倾角必须在0至180度之间.  See coe2rv.');   
end
if((coe(2)*cos(coe(6)))<-1.0)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '不可能达到的双曲轨道.  See coe2rv.');    
end

p=coe(1)*abs(1-coe(2)^2);%半通径
if(abs(coe(2)-1.0)<epsilon)%如果是抛物线轨道,区别对待.
    p=2.0*coe(1);
end

sini=sin(coe(3));cosi=cos(coe(3));
sinO=sin(coe(4));cosO=cos(coe(4));
sino=sin(coe(5));coso=cos(coe(5));

%轨道平面法向单位矢量,即角动量单位矢量
HVector=[ sini*sinO;-sini*cosO; cosi];
      
%偏心率单位矢量,或叫Laplace矢量
PVector=[ cosO*coso-sinO*sino*cosi; sinO*coso+cosO*sino*cosi; sino*sini];

%半通径方向单位矢量,PVector,QVector,HVector构成右手坐标系
% QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
QVector=cross(HVector,PVector);

if((coe(2)*cos(coe(6)))+1.0<epsilon)
    warning('ORBIT:coe2rv:NotSuitableInputs',...
              '抛物或双曲轨道到达无穷远处.  See coe2rv.');
    r=p/epsilon;
else
    r=p/(1+coe(2)*cos(coe(6)));
end

CartesianR=r*(cos(coe(6))*PVector+sin(coe(6))*QVector);
CartesianV=sqrt(mu/p)*(-sin(coe(6))*PVector+(cos(coe(6))+coe(2))*QVector);
if nargout==1
    CartesianR(4:6)=CartesianV;
end