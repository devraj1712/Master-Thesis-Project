% FELP - Gaussian Quadature Points and Weights
% Roger Sauer - 12.12.05


% 1D Quadrature
if jgaus == 1
   gfa = 0.0;
   xig1u(1,1) = gfa ;
   xig1u(1,2) = 2.d0 ;
end

% 2 GP Quadrature
if jgaus == 2
  gfa = 1.d0/sqrt(3.d0) ;
  xig1u(1,1) = -gfa ;
  xig1u(1,2) = 1.d0 ;
  xig1u(2,1) = gfa  ;
  xig1u(2,2) = 1.d0 ;
    
% 3 GP Quadrature
elseif jgaus == 3
  gfa = sqrt(3.d0/5.d0) ;
  xig1u(1,1) = -gfa      ;
  xig1u(1,2) = 5.d0/9.d0 ;
  xig1u(2,1) = 0.d0      ;
  xig1u(2,2) = 8.d0/9.d0 ;
  xig1u(3,1) = gfa       ;
  xig1u(3,2) = 5.d0/9.d0 ;

% 5 GP Quadrature
elseif jgaus == 5
  xig1u(1,1) = -0.906179845938664 ;
  xig1u(1,2) =  0.236926885056189 ;
  xig1u(2,1) = -0.538469310105683 ;
  xig1u(2,2) =  0.478628670499366 ;
  xig1u(3,1) =  0                 ;   
  xig1u(3,2) =  0.568888888888889 ;  
  xig1u(4,1) =  0.538469310105683 ;
  xig1u(4,2) =  0.478628670499366 ;       
  xig1u(5,1) =  0.906179845938664 ;
  xig1u(5,2) =  0.236926885056189 ;
  
% 10 GP Quadrature
elseif jgaus == 10
  xig1u(1,1)  = -0.973906528517172 ;
  xig1u(1,2)  =  0.066671344308688 ;
  xig1u(2,1)  = -0.865063366688985 ;
  xig1u(2,2)  =  0.149451349150581 ;
  xig1u(3,1)  = -0.679409568299024 ;
  xig1u(3,2)  =  0.219086362515982 ;   
  xig1u(4,1)  = -0.433395394129247 ;
  xig1u(4,2)  =  0.269266719309996 ;       
  xig1u(5,1)  = -0.148874338981631 ;
  xig1u(5,2)  =  0.295524224714753 ;
  xig1u(6,1)  =  0.148874338981631 ;
  xig1u(6,2)  =  0.295524224714753 ;
  xig1u(7,1)  =  0.433395394129247 ;
  xig1u(7,2)  =  0.269266719309996 ; 
  xig1u(8,1)  =  0.679409568299024 ;
  xig1u(8,2)  =  0.219086362515982 ;
  xig1u(9,1)  =  0.865063366688985 ;
  xig1u(9,2)  =  0.149451349150581 ;
  xig1u(10,1) =  0.973906528517172 ;
  xig1u(10,2) =  0.066671344308688 ;  

else
  %break  
  for i = 1:jgaus
    xig1u(i,1) = i*2/jgaus-1/jgaus-1 ;
    xig1u(i,2) = 2/jgaus ;
  end
end

% Guass points along v
jgausv = jgaus; %(q_a+1);
for i = 1:jgausv
    xig1v(i,1) = i*2/jgausv-1/jgausv-1 ;
    xig1v(i,2) = 2/jgaus ;
end       

% 2D Quadrature
gp = 0 ;
for j = 1:jgausv
  for i = 1:jgaus
    gp = gp + 1 ;
    xig2(gp,1) = xig1u(i,1) ;
    xig2(gp,2) = xig1v(j,1) ;     
    xig2(gp,3) = xig1u(i,2)*xig1v(j,2) ;
  end
end