function [g,geq]=nonlcon(x)

[sigma, Q] = sol_TenBarTruss(x(1),x(2));
% max(abs(sigma))��X����ȫ�̤j��;
% sqrt(Q(3)^2+Q(4)^2)�p���I2�첾;

g(1) = max(abs(sigma))-250*10^6; %(�D�u��)�������Ag(1)<=0
g(2) = sqrt(Q(3)^2+Q(4)^2)-0.02; %g(2)<=0

geq=[];%�����������A�b�����D���S�������������A�G���Ŷ��X


