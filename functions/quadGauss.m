function gpArray = quadGauss(nGP, dim)
    
    gpArray1 = zeros(nGP, dim+1);
    if nGP == 1
        gpArray1(1,1) = 0.d0;
        gpArray1(1,2) = 2.d0;
    elseif nGP == 2
        gpArray1(:,1) = [-0.5773502691896257,    0.5773502691896257];
        gpArray1(:,2) = [1.d0    1.d0];
    elseif nGP == 3
        gpArray1(:,1) = [-0.7745966692414834     0.d0    0.7745966692414834];
        gpArray1(:,2) = [0.5555555555555556  0.8888888888888888  0.5555555555555556];
    elseif nGP == 4
        gpArray1(:,1) = [-0.861136311594053  -0.339981043584856  0.339981043584856   0.861136311594053];
        gpArray1(:,2) = [ 0.347854845137454  0.652145154862546   0.652145154862546   0.347854845137454];
    elseif nGP == 5
        gpArray1(:,1) = [-0.906179845938664  -0.538469310105683  0.d0    0.538469310105683   0.906179845938664];
        gpArray1(:,2) = [0.236926885056189   0.478628670499367   0.568888888888889   0.478628670499367   0.236926885056189];
    elseif nGP == 6
        gpArray1(:,1) = [-0.932469514203152  -0.661209386466264  -0.238619186083197...
        0.238619186083197   0.661209386466264   0.932469514203152];
        gpArray1(:,2) = [0.171324492379170   0.360761573048139   0.467913934572691...
        0.467913934572691   0.360761573048139   0.171324492379170];
    elseif nGP == 7
        gpArray1(:,1) = [-0.949107912342758  -0.741531185599394  -0.405845151377397   0.d0  0.405845151377397    0.741531185599394   0.949107912342758];
        gpArray1(:,2) = [0.129484966168870   0.279705391489277   0.381830050505119   0.417959183673469   0.381830050505119   0.279705391489277   0.129484966168870];
    elseif nGP == 8
        gpArray1(:,1) = [-0.960289856497536  -0.796666477413627  -0.525532409916329  -0.183434642495650   0.183434642495650  0.525532409916329   0.796666477413627   0.960289856497536];
        gpArray1(:,2) = [0.101228536290376   0.222381034453375   0.313706645877887   0.362683783378362   0.362683783378362   0.313706645877887   0.222381034453375   0.101228536290376];
    elseif nGP == 9
        gpArray1(:,1) = [-0.968160239507626  -0.836031107326636  -0.613371432700590  -0.324253423403809  0.d0   0.324253423403809    0.613371432700590   0.836031107326636     0.968160239507626];
        gpArray1(:,2) = [0.081274388361574   0.180648160694857   0.260610696402935   0.312347077040003   0.330239355001260   0.312347077040003   0.260610696402935   0.180648160694857   0.081274388361574];
    elseif nGP == 10
        gpArray1(:,1) = [-0.973906528517172  -0.865063366688985  -0.679409568299024  -0.433395394129247  -0.148874338981631  0.148874338981631   0.433395394129247   0.679409568299024   0.865063366688985   0.973906528517172];
        gpArray1(:,2) = [0.066671344308688   0.149451349150581   0.219086362515982   0.269266719309996   0.295524224714753   0.295524224714753   0.269266719309996   0.219086362515982   0.149451349150581   0.066671344308688];
    else
        for i = 1:nGP
        gpArray1(i,1) = i*2/nGP - 1/nGP - 1;
        gpArray1(i,2) = 2/nGP;
        end
    end
    
    if dim == 1 % One-dimensional
        gpArray = gpArray1;
    elseif dim == 2 % Two-dimensional
        gpArray = zeros(nGP^2, dim+1);
        gp = 1;
        for j=1:nGP
            for i=1:nGP
                gpArray(gp,1) = gpArray1(i,1);
                gpArray(gp,2) = gpArray1(j,1);
                gpArray(gp,3) = gpArray1(i,2)*gpArray1(j,2);
                gp            = gp+1;
            end
        end
    elseif dim == 3 % Three-dimensional
        gpArray = zeros(nGP^3, dim+1);
        gp = 1 ;
        for k = 1:nGP
            for j = 1:nGP
                for i = 1:nGP
                    gpArray(gp,1) = gpArray1(i,1) ;
                    gpArray(gp,2) = gpArray1(j,1) ;
                    gpArray(gp,3) = gpArray1(k,1) ;
                    gpArray(gp,4) = gpArray1(i,2) * gpArray1(j,2) * gpArray1(k,2) ;
                    gp            = gp + 1 ;
                end
            end
        end
    end
end