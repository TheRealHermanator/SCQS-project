using Plots
BTRG = Dict("symmetric" => 487.134, "Not-symmetric" => 958.923)
TRG = Dict("symmetric" => 527.397, "Not-symmetric" => 933.966)
HOTRG = Dict("symmetric" => 1669, "Not-symmetric" => 3306)
for (key, value) in BTRG
    if key == "symmetric"
        scatter!((value, 1), label=key, yticks=(0:0.2:1,["BTRG", "TRG", "HOTRG"]), legend=:bottomright, mc=:red)
    else
        scatter!((value, 1), label=key, yticks=(0:0.2:1,["BTRG", "TRG", "HOTRG"]), legend=:bottomright, mc=:blue)
    end
end
for (key, value) in TRG
    if key == "symmetric"
        scatter!((value, 1.5), mc=:red, label="")
    else
        scatter!((value, 1.5), mc=:blue, label="")
    end
end
for (key, value) in HOTRG
    if key == "symmetric"
        scatter!((value, 2), mc=:red, label="")
    else
        scatter!((value, 2), mc=:blue, label="")
    end
end
annotate!(200, 1, text("BTRG", :red,11))
annotate!(200, 1.5, text("TRG", :red,11))
annotate!(225, 2, text("HOTRG", :red,11))
title!("Time benchmarks")
xlabel!("Time in ms")
savefig("timebench.png")  