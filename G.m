function G=gravitational(Y)
    Zg = 0.0196;   
    Xg = 0;
    Yg=  0;
    Zb = 0;   
    Xb = -0.611;
    Yb=  0;
    W=299;
    B=308;
    
G=[(W-B)*sin(Y(5)); 
    -(W-B)*cos(Y(5))*sin(Y(4));
    -(W-B)*cos(Y(5))*cos(Y(4)); 
    -(Yg*W-Yb*B)*cos(Y(5))*cos(Y(5))+(Zg*W-Zb*B)*cos(Y(5))*sin(Y(4)); 
    -(Zg*W-Zb*B)*sin(Y(5))+(Xg*W-Xb*B)*cos(Y(5))*cos(Y(4));
    -(Xg*W-Xb*B)*cos(Y(5))*sin(Y(4))-(Yg*W-Yb*B)*sin(Y(5))  ];