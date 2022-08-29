clear all
clc

% doc fmincon %ノ doc 琩高 fmincon ㄏノよΑ

x0=[0.01,0.01]; %癬﹍翴

A=[]; %絬┦ぃ单Α╇兵ン玒计痻皚
b=[]; %絬┦ぃ单Α╇兵ン玒计秖 AX <= b

Aeq=[]; %絬┦单Α╇兵ン玒计痻皚
beq=[]; %絬┦单Α╇兵ン玒计秖 AeqX = beq

ub=[0.5;0.5]; %砞璸丁 upper bounds
lb=[0.001;0.001]; %砞璸丁 lower bounds

options = optimset ('display','off','Algorithm','sqp');%簍衡猭把计砞﹚

[x,fval,exitflag]=fmincon(@(x)obj(x),x0,A,b,Aeq,beq,lb,ub,@(x)nonlcon(x),options);
%磅︽fmincon块程ㄎ秆,x; 程ㄎヘ夹ㄧ计,fval; Μ滥薄,exitflag
%objヘ夹ㄧ计nonlcon  (獶絬┦) ╇兵ン
