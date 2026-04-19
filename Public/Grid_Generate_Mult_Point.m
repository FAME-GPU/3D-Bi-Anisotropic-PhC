function [ Bodycenter_coef, Standard_coef,Edgecenter_coef_x,Edgecenter_coef_y, Edgecenter_coef_z,Facecenter_coef_x,Facecenter_coef_y,Facecenter_coef_z] = Grid_Generate_Mult_Point( Dim )
    n1 = Dim.n1; m1 = Dim.m1;
    n2 = Dim.n2; m2 = Dim.m2;
    n3 = Dim.n3; m3 = Dim.m3;  
    %% Create the Yee's scheme point
    % Creating the grids on the vertices
    r1_coef = (Dim.rmc_1D.x : n1-1)'/n1;
    r2_coef = (Dim.rmc_1D.y : n2-1)'/n2;
    r3_coef = (Dim.rmc_1D.z : n3-1)'/n3;
    
    r1_coef_shift = ((0 : n1-1)' + .5)/n1;%边心
    r2_coef_shift = ((0 : n2-1)' + .5)/n2;
    r3_coef_shift = ((0 : n3-1)' + .5)/n3;
    %% Constructing standard  points
    Standard_coef(:,1) = kron( ones( m3*m2 , 1) , r1_coef );
    Standard_coef(:,2) = kron( ones( m3, 1) , kron( r2_coef , ones( m1, 1)  )  );
    Standard_coef(:,3) = kron( r3_coef, ones( m2*m1, 1) );
    %% Constructing points at cell centers
    Bodycenter_coef(:,1) = kron( ones( n3*n2 , 1) , r1_coef_shift );
    Bodycenter_coef(:,2) = kron( ones( n3 , 1) , kron( r2_coef_shift , ones( n1 , 1)  )  );
    Bodycenter_coef(:,3) = kron( r3_coef_shift , ones( n2*n1 , 1) );
    %% Constructing points on edge for the electric field and on face for the magnetic field
    % Constructing points on edge for the electric field
    Edgecenter_coef_x(:,1) = kron( ones( m3*m2 , 1) , r1_coef_shift );
    Edgecenter_coef_x(:,2) = kron( ones( m3 , 1) , kron( r2_coef       , ones( n1 , 1)  )  );
    Edgecenter_coef_x(:,3) = kron( r3_coef       , ones( m2*n1 , 1) );

    Edgecenter_coef_y(:,1) = kron( ones( m3*n2 , 1) , r1_coef       );
    Edgecenter_coef_y(:,2) = kron( ones( m3 , 1) , kron( r2_coef_shift , ones( m1 , 1)  )  );
    Edgecenter_coef_y(:,3) = kron( r3_coef       , ones( n2*m1 , 1) );

    Edgecenter_coef_z(:,1) = kron( ones( n3*m2 , 1) , r1_coef       );
    Edgecenter_coef_z(:,2) = kron( ones( n3 , 1) , kron( r2_coef       , ones( m1 , 1)  )  ); 
    Edgecenter_coef_z(:,3) = kron( r3_coef_shift , ones( m2*m1 , 1) );

    % Constructing points on face for the magnetic field
    Facecenter_coef_x(:,1) = kron( ones(  n3*n2 , 1) , r1_coef       );
    Facecenter_coef_x(:,2) = kron( ones( n3 , 1) , kron( r2_coef_shift , ones( m1 , 1)  )  );
    Facecenter_coef_x(:,3) = kron( r3_coef_shift , ones(  n2*m1 , 1) );

    Facecenter_coef_y(:,1) = kron( ones(  n3*m2 , 1) , r1_coef_shift );
    Facecenter_coef_y(:,2) = kron( ones( n3 , 1) , kron( r2_coef       , ones( n1 , 1)  )  );
    Facecenter_coef_y(:,3) = kron( r3_coef_shift , ones(  m2*n1 , 1) );

    Facecenter_coef_z(:,1) = kron( ones(  m3*n2 , 1) , r1_coef_shift );
    Facecenter_coef_z(:,2) = kron( ones( m3 , 1) , kron( r2_coef_shift , ones( n1 , 1)  )  ); 
    Facecenter_coef_z(:,3) = kron( r3_coef       , ones(  n2*n1 , 1) );
end