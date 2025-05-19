using TNRKit
import TNRKit: step!
using LinearAlgebra
using LoggingExtras


abstract type TNRScheme end


function run_sing!(scheme, trscheme, criterion;
              finalize_beginning=true, verbosity=1)
    data = []
    singular_values = []

    LoggingExtras.withlevel(; verbosity) do
        @infov 1 "Starting simulation\n $(scheme)\n"
        if finalize_beginning
            push!(data, scheme.finalize!(scheme))
        end

        steps = 0
        crit = true

        t = @elapsed while crit
            @infov 2 "Step $(steps + 1), data[end]: $(!isempty(data) ? data[end] : "empty")"
            step!(scheme, trscheme)
            push!(data, scheme.finalize!(scheme))
            
            # Calculate singular values
            T = scheme.T
            U, S, V = tsvd(T; trunc=trscheme)
            push!(singular_values, S.data)
            steps += 1
            crit = criterion(steps, data)
        end

        #@infov 1 "Simulation finished\n $(stopping_info(criterion, steps, data))\n Elapsed time: $(t)s\n Iterations: $steps"
        # @infov 1 "Elapsed time: $(t)s"
    end
    return data, singular_values
end
