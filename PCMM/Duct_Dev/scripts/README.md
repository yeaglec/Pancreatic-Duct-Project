Output from first run: 

julia> priors = Dict(
           "seg_lin"                   => (0.01,   0.1),
           "home_lin"                  => (0.0005, 0.005),
           "lumenal_pressure_strength" => (0.005,  0.05))
Dict{String, Tuple{Float64, Float64}} with 3 entries:
  "home_lin"                  => (0.0005, 0.005)
  "seg_lin"                   => (0.01, 0.1)
  "lumenal_pressure_strength" => (0.005, 0.05)

julia> book  = Set(["weight","distance","monad_id","sim_id"])
Set{String} with 4 elements:
  "weight"
  "distance"
  "monad_id"
  "sim_id"

julia> pcols = [c for c in names(df) if !(c in book) && any(k -> endswith(c, k), keys(priors))]
3-element Vector{String}:
 "user_parameters/seg_lin"
 "user_parameters/home_lin"
 "user_parameters/lumenal_pressure_strength"

julia> w = Float64.(weights); w ./= sum(w)
50-element Vector{Float64}:
 0.021249014670100125
 0.01687660359763608
 0.018513569105615548
 0.015943292992039452
 0.028596324566662498
 0.020196481843166487
 0.016446072501160376
 0.019167343777497337
 0.016570461002091023
 0.017954392227611488
 0.015437418719324818
 0.02876483529331828
 ⋮
 0.015089065189920382
 0.015065441848678904
 0.018915794231706914
 0.020370008267375863
 0.016168613894322965
 0.018189253053203045
 0.02361618657253772
 0.017428838076106134
 0.02109554428441015
 0.0163203289041065
 0.022108956986802986

julia> @printf("\n%-28s %11s %11s %11s %9s\n","param","post_mean","post_std","prior_std","std_frac")

param                          post_mean    post_std   prior_std  std_frac

julia> for c in pcols
           x = Float64.(df[!, c]); m = sum(w .* x); sd = sqrt(max(sum(w .* (x .- m).^2), 0.0))
           k = first(kk for kk in keys(priors) if endswith(c, kk)); a,b = priors[k]; ps = (b-a)/sqrt(12)
           @printf("%-28s %11.5g %11.5g %11.5g %9.2f\n", c, m, sd, ps, ps>0 ? sd/ps : NaN)
       end
user_parameters/seg_lin         0.049294    0.020503    0.025981      0.79
user_parameters/home_lin      0.00091632  0.00031915    0.001299      0.25
user_parameters/lumenal_pressure_strength    0.034101    0.011351     0.01299      0.87

julia> println("  std_frac -> 1: unconstrained;  << 1: constrained")
  std_frac -> 1: unconstrained;  << 1: constrained

julia> println("\nCorrelations:")

Correlations:

julia> X = hcat([Float64.(df[!, c]) for c in pcols]...)
50×3 Matrix{Float64}:
 0.0472036  0.00053092   0.0381764
 0.0525137  0.000795216  0.039032
 0.0251867  0.000752689  0.0262082
 0.0649503  0.00127098   0.0265785
 0.0718898  0.000684418  0.0127525
 0.026088   0.000642704  0.0335495
 0.0729325  0.00146073   0.0400489
 0.0328212  0.000734545  0.0165155
 0.030411   0.000911866  0.0270712
 0.0402256  0.000597586  0.0311424
 0.0611326  0.000918215  0.0355118
 0.0254251  0.000588176  0.0444736
 ⋮                       
 0.0643475  0.000998273  0.0328493
 0.0399497  0.000904269  0.0315356
 0.076006   0.00133248   0.0468942
 0.0743719  0.00105569   0.0467903
 0.0531713  0.00104242   0.042989
 0.0402574  0.00143093   0.0259945
 0.0528201  0.000588143  0.0439627
 0.0655273  0.000773503  0.0292679
 0.0485601  0.00180607   0.0377479
 0.0688048  0.00153739   0.0384534
 0.038576   0.000832989  0.0471269

julia> show(stdout, "text/plain", DataFrame(cor(X), pcols)); println()
3×3 DataFrame
 Row │ user_parameters/seg_lin  user_parameters/home_lin  user_parameters/lumenal_pressure_strength 
     │ Float64                  Float64                   Float64                                   
─────┼──────────────────────────────────────────────────────────────────────────────────────────────
   1 │                1.0                       0.345294                                   0.192678
   2 │                0.345294                  1.0                                        0.355235
   3 │                0.192678                  0.355235                                   1.0

julia> d = Float64.(df.distance); bi = argmin(d)
10

julia> @printf("\ndistances  min=%.4g  median=%.4g  max=%.4g\n", minimum(d), median(d), maximum(d))

distances  min=0.0009638  median=0.04119  max=0.09988

julia> @printf("best particle: monad %d, distance %.4g\n", df.monad_id[bi], d[bi])
best particle: monad 357, distance 0.0009638

julia> for c in pcols; @printf("   %-28s = %.5g\n", c, df[bi, c]); end
   user_parameters/seg_lin      = 0.040226
   user_parameters/home_lin     = 0.00059759
   user_parameters/lumenal_pressure_strength = 0.031142

julia> # best particle's metrics vs target (reads output files, runs no sims):
       s = bm_summary_statistic(df.monad_id[bi]; target_cell_count = 155)
Computing SimulationPopulationTimeSeries for Simulation 712...done.
Computing SimulationPopulationTimeSeries for Simulation 713...done.
Dict{String, Float64} with 14 entries:
  "max_indent_depth" => 0.255035
  "perimeter_t0"     => 1788.04
  "area_frac_change" => -0.0906191
  "max_displacement" => 106.276
  "IC"               => 1.41116
  "eval_timestep"    => 53.0
  "IC_tf"            => 1.41116
  "IC_t0"            => 1.07978
  "area_tf"          => 2.14266e5
  "perimeter_tf"     => 1949.26
  "lobe_amp"         => 0.0843376
  "indent_extent"    => 0.339844
  "area_t0"          => 2.35618e5
  "roughness"        => 0.0381489

julia> for (k, tgt) in observed_target
           haskey(s, k) && @printf("%-20s sim=%.4f  target=%.4f  resid=%+.4f\n", k, s[k], tgt, s[k]-tgt)
       end
IC                   sim=1.4112  target=1.4100  resid=+0.0012
area_frac_change     sim=-0.0906  target=-0.0950  resid=+0.0044