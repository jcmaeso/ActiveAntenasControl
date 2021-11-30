classdef ActiveAntenna
    %ACTIVEANTENNA Summary of this class goes here
    %   Detailed explanation goes here
    properties
        element_separation
        n_elements
        n_bits
        amplitude_resolution
        attenuation %Sets amplitude discretization (Min attenuation = 0 or MAX VAL)
        positions
    end

    methods
        function obj = ActiveAntenna(n_elements,element_separation,n_bits,amplitude_resolution,atenuation)
            %ACTIVEANTENNA Construct an instance of this class
            %   Detailed explanation goes here
            obj.element_separation = element_separation;
            obj.n_elements = n_elements;
            obj.n_bits = n_bits;
            obj.amplitude_resolution = amplitude_resolution;
            obj.positions = linspace(-(n_elements-1)/2,(n_elements-1)/2,n_elements);
            obj.attenuation = boolean(atenuation); 
            % If attenuation = 0, max bit val  is minimum atenuation, otherwase min bit val is minimum atenuation
        end

        function [phases_dis,phases_real,phases_chip] = aim2dir(obj,direction,initial_phase)
            %aim2dir Generates phase coeficients for a given direction
            %   phases_dis: Discretized phases to given chip
            %   phases_real: Real phases synthethised
            %   initial_phase_state
            if nargin < 3
                initial_phase = 0;
            end
            prog_phase = rad2deg(2*pi*obj.element_separation*sind(direction));
            phases_real = wrapTo360(prog_phase.*(0:(obj.n_elements-1))+initial_phase);
            phases_dis = obj.discretize_phase(phases_real);
            phases_chip = phases_dis./(2^obj.n_bits);
        end

        function [amp_dis,amp_real,amp_chip] = amps_cosine_over_pedestal(obj,H,N)
            amps = 1+H*cos(pi*obj.positions/(N-1)).^2;
            [amp_dis,amp_real,amp_chip] = obj.normalize_and_discretize_amplitude(amps);
        end

        function [amp_dis,amp_real,amp_chip] = amps_parabolic(obj,delta)
            amps = 1-(1-delta)*(obj.positions/3.5).^2;
            [amp_dis,amp_real,amp_chip] = obj.normalize_and_discretize_amplitude(amps);
        end

        function [amp_dis,amp_real,amp_chip] = amps_triangular(obj,step)
            amps = zeros(1,obj.n_elements);
            if mod(obj.n_elements,2) == 0
                amps(1:(obj.n_elements/2)) =  1+step*obj.positions(1:obj.n_elements/2)/obj.n_elements;
                amps((obj.n_elements/2+1):end) = flip(amps(1:(obj.n_elements/2))); 
            else
                amps(1:round(obj.n_elements/2)) =  1+step*obj.positions(1:round(obj.n_elements/2))/obj.n_elements;
                amps((round(obj.n_elements/2)+1):end) =  flip(amps(1:round(obj.n_elements/2)));
            end
            [amp_dis,amp_real,amp_chip] = obj.normalize_and_discretize_amplitude(amps);
        end

        function [amp_dis,amp_real,amp_chip] = amps_cheb(obj,sll)
            amps = chebwin(obj.n_elements,sll);
            [amp_dis,amp_real,amp_chip] = obj.normalize_and_discretize_amplitude(amps);
        end

        function [amp_dis,amp_real,amp_chip] = amps_taylor(obj,sll,n_sll)
            amps = taylorwin(obj.n_elements,n_sll,-sll);
            [amp_dis,amp_real,amp_chip] = obj.normalize_and_discretize_amplitude(amps);
        end
    end
    methods (Access = private)
        function [pha_val,pha_chip] = discretize_phase(obj,phases)
            %discretize_phase Discretization of the given phase list to
            %best suited state
            % INPUT:
            %   phases: phases list to be discretized
            % OUTPUT:
            %   pha_val: Discretized value of the phase
            %   pha_chip: Chip bit state for phase
            scale = 0:(360/(2^obj.n_bits)):(360-360/(2^obj.n_bits));
            pha_chip = zeros(size(phases));
            pha_val = zeros(size(phases));
            for i = 1:length(phases)
                [~,pha_chip(i)] = min(abs(scale-phases(i)));
                pha_val(i) = scale(pha_chip(i));
            end
        end

        function [amp_dis,amp_real,amp_chip] = normalize_and_discretize_amplitude(obj,amps)
            amp_real = 20*log10(amps)-max(20*log10(amps));
            amp_dis = obj.discretize_amplitude(amp_real);
            if obj.attenuation
                amp_chip = (2^obj.n_bits) + amp_dis/obj.amplitude_resolution; 
            else
                amp_chip = -amp_dis/obj.amplitude_resolution; 
            end
        end

        function [amp_val,amp_chip] = discretize_amplitude(obj,amplitudes)
            bits = 0:(2^obj.n_bits-1);
            scale = bits*(-obj.amplitude_resolution);
            amp_chip = zeros(size(amplitudes));
            amp_val = zeros(size(amplitudes));
            for i = 1:length(amplitudes)
                [~,amp_chip(i)] = min(abs(scale-amplitudes(i)));
                amp_val(i) = scale(amp_chip(i));
            end
        end
    end
end

