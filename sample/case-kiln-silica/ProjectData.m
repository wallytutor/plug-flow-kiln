classdef ProjectData
    methods (Access = public, Static = true)
        function [val] = thickness_coating(z)
            % Coating thickness over length [m].
            cl = 6.0;
            t0 = 0.0254 * 3.00;
            t1 = 0.0254 * 0.01;
            t2 = 0.0254 * 0.01;

            m0 = (t1 - t0) / cl;
            m1 = (t2 - t1) / (z(end)  - cl);

            rng0 = z <  cl;
            rng1 = z >= cl;

            z0 = z(rng0);
            z1 = z(rng1) - cl;

            val = zeros(size(z));
            val(rng0) = t0 + m0 .* z0;
            val(rng1) = t1 + m1 .* z1;
        end

        function [val] = thickness_refractory(z)
            % Refractory thickness over length [m].
            val = NA;
        end

        function [val] = thickness_shell(z)
            % Shell thickness over length [m].
            val = NA;
        end

        function [val] = conductivity_coating(z, T)
            % Coating thermal conductivity over length [W/(m².K)].
            val = NA;
        end

        function [val] = conductivity_refractory(z, T)
            % Refractory thermal conductivity over length [W/(m².K)].
            val = NA;
        end

        function [val] = conductivity_shell(z, T)
            % Shell thermal conductivity over length [W/(m².K)].
            val = NA;
        end

        function [val] = environment_htc(z)
            % Environment heat transfer coefficient [W/(m.K)].
            val = 10.0 * ones(size(z));
        end

        function [val] = environment_temperature(z)
            % Environment temperature [K].
            val = 313.15 * ones(size(z));
        end
    end
end
