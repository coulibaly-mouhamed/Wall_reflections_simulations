 function [rhs2] =  DtN(phi_hat,p)
% Approximation to \phi_z (i.e. Dirichelt-to-Neumann operator)
    switch p.DtN_method

        case 1 % Wave equation approximation
            w = p.d.*ifft2(p.KxiKy.*phi_hat);
            A = fft2(w); 
            As = conj(A(p.shift1,p.shift2));
            rhs2      =   - (p.KxmiKy.*A/2 + p.KxiKy.*As/2);

        case 3 % Exact DtN for flat bottom
            rhs2  = p.abs_K.*tanh(p.abs_K.*p.h).*phi_hat;
    end

end