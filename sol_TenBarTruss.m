function [sigma, Q] = sol_TenBarTruss(r1, r2)
  
    % 定義各參數數值
    
    n1=[18.28 9.14];
    n2=[18.28 0];
    n3=[9.14 9.14];
    n4=[9.14 0];
    n5=[0 9.14];
    n6=[0 0];
     
    
    function A = area(x)   %桿件截面積
    A=pi*x^2;
    end    
    A1=area(r1);
    A2=area(r2);
    AM = [A1; A1; A1; A1; A1; A1; A2; A2; A2; A2];
    
    function L = long(xi,xj,yi,yj)  %桿件長度
    L = sqrt((xj-xi)^2+(yj-yi)^2);
    end
    L1 = long(n3(1),n5(1),n3(2),n5(2));
    L2 = long(n1(1),n3(1),n1(2),n3(2));
    L3 = long(n4(1),n6(1),n4(2),n6(2));
    L4 = long(n2(1),n4(1),n2(2),n4(2));
    L5 = long(n3(1),n4(1),n3(2),n4(2));
    L6 = long(n1(1),n2(1),n1(2),n2(2));
    L7 = long(n4(1),n5(1),n4(2),n5(2));
    L8 = long(n3(1),n6(1),n3(2),n6(2));
    L9 = long(n2(1),n3(1),n2(2),n3(2));
    L10 = long(n1(1),n4(1),n1(2),n4(2));
    LM = [L1; L2; L3; L4; L5; L6; L7; L8; L9; L10];
    

    function theta = theta(xi,xj,yi,yj)  %y座標變量
    theta = atand((yj-yi)/(xj-xi));
    end
    theta1 = theta(n3(1),n5(1),n3(2),n5(2));
    theta2 = theta(n1(1),n3(1),n1(2),n3(2)); 
    theta3 = theta(n4(1),n6(1),n4(2),n6(2));
    theta4 = theta(n2(1),n4(1),n2(2),n4(2));  
    theta5 = theta(n3(1),n4(1),n3(2),n4(2));  
    theta6 = theta(n1(1),n2(1),n1(2),n2(2));  
    theta7 = theta(n4(1),n5(1),n4(2),n5(2));  
    theta8 = theta(n3(1),n6(1),n3(2),n6(2));  
    theta9 = theta(n2(1),n3(1),n2(2),n3(2));  
    theta10 = theta(n1(1),n4(1),n1(2),n4(2));
    thetaM = [theta1; theta2; theta3; theta4; theta5; theta6; theta7; theta8; theta9; theta10];

    % 開一個空白的剛性矩陣 (stiffness matrix)
    E = 200*10^9;
    K = zeros(12,12);
  
    % 計算 stiffness matrix (可使用 add_element 函數)
    K = add_element(K,AM(1,1),E,LM(1,1),thetaM(1,1),3,5);
    K = add_element(K,AM(2,1),E,LM(2,1),thetaM(2,1),1,3);
    K = add_element(K,AM(3,1),E,LM(3,1),thetaM(3,1),4,6);
    K = add_element(K,AM(4,1),E,LM(4,1),thetaM(4,1),2,4);
    K = add_element(K,AM(5,1),E,LM(5,1),thetaM(5,1),3,4);
    K = add_element(K,AM(6,1),E,LM(6,1),thetaM(6,1),1,2);
    K = add_element(K,AM(7,1),E,LM(7,1),thetaM(7,1),4,5);
    K = add_element(K,AM(8,1),E,LM(8,1),thetaM(8,1),3,6);
    K = add_element(K,AM(9,1),E,LM(9,1),thetaM(9,1),2,3);
    K = add_element(K,AM(10,1),E,LM(10,1),thetaM(10,1),1,4);
    K_reduced = K(1:8,1:8);
    K_reaction = K(9:12,1:12);
    K_inv = inv(K_reduced);

    % 建立力矩陣
     F_matrix = [0; 0; 0; 1*10^7; 0; 0; 0; 1*10^7; 0; 0; 0; 0];
     F_reduced = F_matrix(1:8);

    % 建立空白位移矩陣
     Q = [];
     Q_reduced = [];
    % 計算位移量 (F = KQ)
     Q_reduced = K_inv*F_reduced;
     Q = [Q_reduced; 0; 0; 0; 0];

    % 建立空白應力矩陣
     sigma = [];
  
    % 計算應力 (stress) (可使用 compute_stress 函數)
     sigma1 = compute_stress (Q, E, LM(1,1), thetaM(1,1), 3, 5);
     sigma2 = compute_stress (Q, E, LM(2,1), thetaM(2,1), 1, 3);
     sigma3 = compute_stress (Q, E, LM(3,1), thetaM(3,1), 4, 6);
     sigma4 = compute_stress (Q, E, LM(4,1), thetaM(4,1), 2, 4);
     sigma5 = compute_stress (Q, E, LM(5,1), thetaM(5,1), 3, 4);
     sigma6 = compute_stress (Q, E, LM(6,1), thetaM(6,1), 1, 2);
     sigma7 = compute_stress (Q, E, LM(7,1), thetaM(7,1), 4, 5);
     sigma8 = compute_stress (Q, E, LM(8,1), thetaM(8,1), 3, 6);
     sigma9 = compute_stress (Q, E, LM(9,1), thetaM(9,1), 2, 3);
     sigma10 = compute_stress (Q, E, LM(10,1), thetaM(10,1), 1, 4);
     sigma = [sigma1; sigma2; sigma3; sigma4; sigma5; sigma6; sigma7; sigma8; sigma9; sigma10];

    % (optional) compute reactions
     R = K_reaction*Q;
  
end
