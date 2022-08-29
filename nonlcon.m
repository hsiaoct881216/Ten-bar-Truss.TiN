function [g,geq]=nonlcon(x)

[sigma, Q] = sol_TenBarTruss(x(1),x(2));
% max(abs(sigma))找出絕對值後最大值;
% sqrt(Q(3)^2+Q(4)^2)計算點2位移;

g(1) = max(abs(sigma))-250*10^6; %(非線性)拘束條件，g(1)<=0
g(2) = sqrt(Q(3)^2+Q(4)^2)-0.02; %g(2)<=0

geq=[];%等式拘束條件，在本問題中沒有等式拘束條件，故為空集合


