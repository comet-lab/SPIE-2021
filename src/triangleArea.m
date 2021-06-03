function A = triangleArea(p)
    p1 = p(1,:);
    p2 = p(2,:);
    p3 = p(3,:);
    
    a = norm(p1 - p2);
    b = norm(p2 - p3);
    c = norm(p3 - p1);
    
    s = (a + b + c) / 2;
    A = sqrt(s * (s - a) * (s - b) * (s - c));
end
