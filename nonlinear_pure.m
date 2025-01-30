classdef nonlinear_pure
    properties
        n % Neighborhood Size (including self)
        c_v % Cost of vasculature production
        c_i % Cost of immunology production
        s % Archetype 1 growth benefit
        ben_v % Function handle for vasculature benefit
        ben_i % Function handle for immunology benefit
    end

    methods
        function obj = nonlinear_pure(n, c_v, c_i, s, ben_v, ben_i)
            obj.n = n;
            obj.c_v = c_v;
            obj.c_i = c_i;
            obj.s = s;
            obj.ben_v = @(j) ben_v(j, obj.n);
            obj.ben_i = @(j) ben_i(j, obj.n);
        end

        function pr = f_ij(obj, i, j, x)
            x1 = x(1);
            x2 = x(2);
            x3 = 1 - x1 - x2;
            c1 = nchoosek(obj.n-1, i);
            c2 = nchoosek(obj.n-1-i, j);
            pr = c1 * c2 * x1^i * x2^j * x3^(obj.n - 1 - i - j);
        end

        function W1 = W1(obj, x)
            W1 = 0;
            for i = 0:obj.n - 1
                for j = 0:obj.n - 1 - i
                    term = obj.f_ij(i, j, x) * (obj.ben_v(i + j + 1) + obj.ben_i(j));
                    W1 = W1 + term;
                end
            end
            W1 = W1 - obj.c_v + obj.s;
        end

        function W2 = W2(obj, x)
            W2 = 0;
            for i = 0:obj.n - 1
                for j = 0:obj.n - 1 - i
                    term = obj.f_ij(i, j, x) * (obj.ben_v(i + j + 1) + obj.ben_i(j+1));
                    W2 = W2 + term;
                end
            end
            W2 = W2 - obj.c_v - obj.c_i;
        end

        function W3 = W3(obj, x)
            W3 = 0;
            for i = 0:obj.n - 1
                for j = 0:obj.n - 1 - i
                    term = obj.f_ij(i, j, x) * (obj.ben_v(i + j) + obj.ben_i(j));
                    W3 = W3 + term;
                end
            end
        end

    end
end
