function [CartesianR, CartesianV]=coe2rv(coe, mu, epsilon)
%coe2rv ���ݾ�������������Ĺ���ֱ������ϵ�µ�λ�ú��ٶȷ���
%������� coe(6),mu:
% 	coe(1)�볤��a���ף�:������Բ(��Բ)��˫���߹��,��ֵ������;���������߹��,
%                      ��ֵȡΪ���Ǿ�.
% 	coe(2)ƫ����e:����Բ���e=0;��Բ���0<e<1,�����߹��e=1,˫���߹��e>1
% 	coe(3)������i�����ȣ�:��Χ0<=i<=180��.
%	coe(4)������ྭOmega�����ȣ�:��������Ϊ0ʱû������,�ɽ���ֵȡΪ0
%	coe(5)���������omega�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵ��Ϊ0
%	coe(6)������f�����ȣ�:��ƫ����Ϊ0ʱû������,�ɽ���ֵȡΪomega+f,��γ�ȷ���
%   mu:��������������ϵ��,���������������������������ĳ˻�,��λm^3/s^2,
%      Ĭ��Ϊ��������ϵ��3.98600441800e+14
%   epsilon:��С����,Ĭ��Ϊ1.0e-12
%�������[CartesianR, CartesianV]:
%   λ����ʸ��CartesianR����X��Y��Z���ٶ���ʸ��CartesianV�ķ���VX��VY��VZ,
%   ��λ���ס���ÿ��.
%ע:�������Ϊһ��ʱ,��λ�ú��ٶ����ϳ�һ��6ά������

% *             �廪��ѧ���캽��ѧԺ����ѧ������о��Ҳ�ʿ��                 *
% *                ������(jiangfh04@mails.thu.edu.cn)                      *
% *                 ����޸�: 2009.2.24                                    *
%��ʽ��Դ��������������ۡ�P.40��������.

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
              '����������ӦΪ6��.  See coe2rv.');
end
if(mu<=0.0||epsilon<=eps)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '��������ӦΪ��,��С����Ӧ������ϵͳ��С˫����ֵ.  See coe2rv.');   
end
if(coe(1)<0.0||coe(2)<0.0)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '�볤���ƫ���ʶ�����Ϊ��.  See coe2rv.');    
end
if(coe(3)<0.0||coe(3)>pi)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '�����Ǳ�����0��180��֮��.  See coe2rv.');   
end
if((coe(2)*cos(coe(6)))<-1.0)
    error('ORBIT:coe2rv:NotSuitableInputs',...
              '�����ܴﵽ��˫�����.  See coe2rv.');    
end

p=coe(1)*abs(1-coe(2)^2);%��ͨ��
if(abs(coe(2)-1.0)<epsilon)%����������߹��,����Դ�.
    p=2.0*coe(1);
end

sini=sin(coe(3));cosi=cos(coe(3));
sinO=sin(coe(4));cosO=cos(coe(4));
sino=sin(coe(5));coso=cos(coe(5));

%���ƽ�淨��λʸ��,���Ƕ�����λʸ��
HVector=[ sini*sinO;-sini*cosO; cosi];
      
%ƫ���ʵ�λʸ��,���Laplaceʸ��
PVector=[ cosO*coso-sinO*sino*cosi; sinO*coso+cosO*sino*cosi; sino*sini];

%��ͨ������λʸ��,PVector,QVector,HVector������������ϵ
% QVector=[-cosO*sino-sinO*coso*cosi;-sinO*sino+cosO*coso*cosi;coso*sini];
QVector=cross(HVector,PVector);

if((coe(2)*cos(coe(6)))+1.0<epsilon)
    warning('ORBIT:coe2rv:NotSuitableInputs',...
              '�����˫�������������Զ��.  See coe2rv.');
    r=p/epsilon;
else
    r=p/(1+coe(2)*cos(coe(6)));
end

CartesianR=r*(cos(coe(6))*PVector+sin(coe(6))*QVector);
CartesianV=sqrt(mu/p)*(-sin(coe(6))*PVector+(cos(coe(6))+coe(2))*QVector);
if nargout==1
    CartesianR(4:6)=CartesianV;
end