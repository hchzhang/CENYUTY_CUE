function [som_temp,litt_temp] = control_temp_func (temp_in, frozen_respiration_func,tem,x)

soil_Q10=0.69; % SOC mixed with fresh litter 
% soil_Q10=1.11;  % SOC mixed with pre-incubated litter
tsoil_ref = 30; % Fresh litter

litt_Q10 = 0.69; % fresh litter
% litt_Q10=0.52; % Inhibited litter
tlitt_ref = 30;

Q10 = 10. ;
ZeroCelsius = 273.15;
%frozen_respiration_func = 1; 
%temp_in in celsius 
%soil_Q10=2;
switch(frozen_respiration_func)
    
    case(0) % this is the standard ORCHIDEE state
        som_result(:) = exp( soil_Q10 * ( temp_in - (ZeroCelsius+tsoil_ref)) / Q10 );
        som_result(:) = min( 1.0, som_result(:) );
        
    case(1)  % cutoff respiration when T < -1C
        % SOM
        if (temp_in > ZeroCelsius ) % normal as above
            som_result = exp( soil_Q10 * ( temp_in - (ZeroCelsius+tsoil_ref) ) / 10. );
        elseif (temp_in > ZeroCelsius - 1. )
            som_result = (temp_in - (ZeroCelsius - 1.)) * ...
                exp( soil_Q10 * ( ZeroCelsius - (ZeroCelsius+tsoil_ref) ) / 10. );
        else
            som_result = 0.0;
        end
        
%         tempfunc_result = max(min( 1.0, tempfunc_result ), 0);
        som_result = max(som_result,0);
        % litter
        if (temp_in > ZeroCelsius ) % normal as above
            litt_result = exp( litt_Q10 * ( temp_in - (ZeroCelsius+tlitt_ref) ) / 10. );
        elseif (temp_in > ZeroCelsius - 1. )
            litt_result = (temp_in - (ZeroCelsius - 1.)) * ...
                exp( litt_Q10 * ( ZeroCelsius - (ZeroCelsius+tlitt_ref) ) / 10. );
        else
            litt_result = 0.0;
        end
        
%         tempfunc_result = max(min( 1.0, tempfunc_result ), 0);
        litt_result = max(litt_result,0);
        
    case(2)  % cutoff respiration when T < -3C
        if (temp_in > ZeroCelsius )
            som_result = exp( soil_Q10 * ( temp_in - (ZeroCelsius+tsoil_ref) ) / 10. );
        elseif (temp_in > ZeroCelsius - 3. )
            som_result = ((temp_in - (ZeroCelsius - 3.))/3.)...
                * exp( soil_Q10 * ( ZeroCelsius - (ZeroCelsius+tsoil_ref) ) / 10. );
        else
            som_result = 0.0;
        end
        
    case(3)  % q10 = 100 when below zero
        if  (temp_in(:) > ZeroCelsius )
            som_result = exp( soil_Q10 * ( temp_in - (ZeroCelsius+tsoil_ref) ) / 10. );
        else
            som_result = exp( 4.605 * ( temp_in - (ZeroCelsius) ) / 10.) ...
                * exp( soil_Q10 * ( -tsoil_ref ) / 10. );
        end
        
    case(4)  % q10 = 1000 when below zero
        if (temp_in(:) > ZeroCelsius )
            som_result = exp( soil_Q10 * ( temp_in - (ZeroCelsius+tsoil_ref) ) / 10. );
        else
            som_result = exp( 6.908 * ( temp_in - (ZeroCelsius) ) / 10.) ...
                * exp( soil_Q10 * ( -tsoil_ref ) / 10. );
        end
        
end
% litter_temp = max(min( 1.0, tempfunc_result ), 0);
som_temp = max(som_result,0);
litt_temp = max(litt_result,0);

end