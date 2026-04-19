function mtx = C(row,col,delta,BC)
    % Curl operator for approximate electric field e
    mtx = [  DistOp.O_x(row,col)          , -DistOp.C_13(row,col,delta,BC),  DistOp.C_12(row,col,delta,BC)  ;
             DistOp.C_23(row,col,delta,BC),  DistOp.O_y(row,col)          , -DistOp.C_21(row,col,delta,BC)  ;
            -DistOp.C_32(row,col,delta,BC),  DistOp.C_31(row,col,delta,BC),  DistOp.O_z(row,col)            ];
end