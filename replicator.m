function xdot = replicator(t,x,A)    
    f1 = A.W1(x);
    f2 = A.W2(x);
    f3 = A.W3(x);
    phi = f1.*x(1) + f2.*x(2) + f3.*x(3);
    
    xdot(1,1) = x(1)* ( f1 - phi ) ; 
    xdot(2,1) = x(2)* ( f2 - phi ) ;
    xdot(3,1) = x(3)* ( f3 - phi ) ;
    
end