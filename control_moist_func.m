function litter_moist = control_moist_func (moist_in)

moist_coeff = [1.1,  2.4,  0.29 ] ;
moistcont_min = 0.25 ;

moistfunc_result = -moist_coeff(1) * moist_in * moist_in +...
    moist_coeff(2)* moist_in - moist_coeff(3);
litter_moist = max( moistcont_min, min( 1, moistfunc_result ) );

end