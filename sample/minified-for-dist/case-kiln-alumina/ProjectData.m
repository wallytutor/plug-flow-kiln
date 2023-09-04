classdef ProjectData
    methods (Access = public, Static = true)
        function [val] = thickness_coating(z)
            % Coating thickness over length [m].
            val = 0.001 * ones(size(z));
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
            val = 1.8;
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
            dist = @(x, c, sig) exp((-(x - c).^2) / (2 * sig.^2));

            % Informed centers are at 2, 4, 7, 10, 11, 20, and 27 m
            % but these values do not match the observations so some
            % of the coolers where slighly shifted below.
            val = 5.0;
            val += 20.0 * dist(z, 2.0,  1.0);
            val += 20.0 * dist(z, 5.0,  1.0);
            val += 20.0 * dist(z, 8.0,  1.0);
            val += 20.0 * dist(z, 11.5, 1.0);
            val += 20.0 * dist(z, 12.5, 1.0);
            val += 20.0 * dist(z, 21.0, 1.0);
            val += 20.0 * dist(z, 27.0, 1.0);
        end

        function [val] = environment_temperature(z)
            % Environment temperature [K].
            val = 313.15 * ones(size(z));
        end
    end
end
