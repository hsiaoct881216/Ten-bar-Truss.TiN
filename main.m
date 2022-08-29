clear all
clc

% doc fmincon %�Q�� doc �d�� fmincon ���ϥΤ覡

x0=[0.01,0.01]; %�_�l�I

A=[]; %�u�ʤ�����������󪺫Y�Ưx�}
b=[]; %�u�ʤ�����������󪺫Y�ƦV�q AX <= b

Aeq=[]; %�u�ʵ���������󪺫Y�Ưx�}
beq=[]; %�u�ʵ���������󪺫Y�ƦV�q AeqX = beq

ub=[0.5;0.5]; %�]�p�Ŷ��� upper bounds
lb=[0.001;0.001]; %�]�p�Ŷ��� lower bounds

options = optimset ('display','off','Algorithm','sqp');%�t��k���ѼƳ]�w

[x,fval,exitflag]=fmincon(@(x)obj(x),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x),options);
%����fmincon�A��X�̨θ�,x; �̨ΥؼШ�ƭ�,fval; ���ı���,exitflag
%obj���ؼШ�ơAnonlcon �� (�D�u��) �������
