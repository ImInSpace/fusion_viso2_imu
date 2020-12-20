syms qx qy qz qw 

eul = [ atan2( 2*(qx.*qy+qw.*qz), qw.^2 + qx.^2 - qy.^2 - qz.^2 ), ...
            asin( -2*(qx.*qz-qw.*qy) ), ...
            atan2( 2*(qy.*qz+qw.*qx), qw.^2 - qx.^2 - qy.^2 + qz.^2 )];
        
print_vec(eul,'eul',[],[],1)