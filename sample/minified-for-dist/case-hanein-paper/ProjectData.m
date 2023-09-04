classdef ProjectData
    methods (Access = public, Static = true)
        function [val] = thickness_coating(z)
            % Coating thickness over length [m].
            val = NA;
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
            val = 0.2475 .* (1.0 + 5.85e-04 .* T);
        end

        function [val] = conductivity_refractory(z, T)
            % Refractory thermal conductivity over length [W/(m².K)].
            val = 0.2475 .* (1.0 + 5.85e-04 .* T);
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
