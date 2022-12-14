function K = add_element(K, A, E, L, theta, node1, node2)
    c = cosd(theta); s = sind(theta);
    temp = A*E/L*[c^2 c*s; c*s s^2];
    K((2*node1-1):(2*node1), (2*node1-1):(2*node1)) = K((2*node1-1):(2*node1), (2*node1-1):(2*node1)) + temp;
    K((2*node2-1):(2*node2), (2*node2-1):(2*node2)) = K((2*node2-1):(2*node2), (2*node2-1):(2*node2)) + temp;
    K((2*node1-1):(2*node1), (2*node2-1):(2*node2)) = K((2*node1-1):(2*node1), (2*node2-1):(2*node2)) - temp;
    K((2*node2-1):(2*node2), (2*node1-1):(2*node1)) = K((2*node2-1):(2*node2), (2*node1-1):(2*node1)) - temp;
end

