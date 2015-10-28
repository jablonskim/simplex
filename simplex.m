## Copyright (C) 2015 jablonskim
## 
## This program is free software; you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 3 of the License, or
## (at your option) any later version.
## 
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## -*- texinfo -*- 
## @deftypefn {Function File} {@var{retval} =} simplex (@var{input1}, @var{input2})
##
## @seealso{}
## @end deftypefn

## Author: jablonskim <jablonskim@LENOVO-PC>
## Created: 2015-10-27

function [x, fval, comment] = simplex (f, A, b)
    n = length(f);
    m = length(b);
    fval = 0; # TODO
    
    # Wsp. funkcji celu
    c = f';
    
    # Indeksy zmiennych bazowych
    x_b = find_brd(A, n, m);
        
    while true
        # Wsp. bazowe funkcji celu
        c_b = c(x_b);
        
        # Wektor wskaznikow optymalnosci
        z = c_b * A;
        z_c = z - c;
        
        # Znalezione rozwiazanie
        [val, k_in] = min(z_c);
        if val >= 0
            x = zeros(n, 1);
            for i = 1:m
                x(x_b(i)) = b(i);
            end
            fval = sum(x' .* c);
            comment = 'OK';
            return;
        end
        
        # Zadanie nie posiada skonczonego RO
        for j = 1:n
            if z_c(j) < 0 && max(A(:, j)) < 0
                x = [];
                comment = 'Brak skonczonego RO.';
                return;
            end    
        end
      
        # Znalezienie indeksu zmiennej usuwanej
        out_values = [];
        
        for i = 1:m
            if A(i, k_in) <= 0
                out_values = [out_values; Inf];
            else
                out_values = [out_values; b(i) / A(i, k_in)];
            end
        end
        
        [tmp, k_out] = min(out_values);
        
        div_val = A(k_out, k_in);
        
        new_A = zeros(m, n);
        new_b = zeros(m, 1);
        
        new_A(k_out, :) = A(k_out, :) / div_val;
        new_b(k_out) = b(k_out) / div_val;
        
        for i = 1:m
            if i != k_out
                new_A(i, :) = A(i, :) - (A(k_out, :) * A(i, k_in) / div_val);
                new_b(i) = b(i) - (A(i, k_in) * b(k_out) / div_val);
            end
        end
        
        A = new_A;
        b = new_b;
        
        x_b(k_out) = k_in;
        
    end
end

function [brd] = find_brd (A, n, m)
    brd = [];
    
    for j = 1:m
        for i = 1:n
            is_ok = true;
            
            for k = 1:m
                if (k == j) && (A(k, i) != 1)
                    is_ok = false;
                    break;
                elseif (k != j) && (A(k, i) != 0)
                    is_ok = false;
                    break;
                end
            end
            
            if is_ok == true
                brd = [brd; i];
            end
        end
    end
end