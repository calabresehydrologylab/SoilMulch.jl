module SoilMulch

# Packages [deps]
using Distributions
using DataFrames
using Statistics

export rainfall_poisson, soil_leakage, plant_transpiration, soil_evaporation
export cover_fraction, mulch_leakage, mulch_evaporation
export soil_water_balance, soil_mulch_water_balance, df_dt_day
export sol_swb, sol_mswb, sol_swb_crop, sol_mswb_crop


# ------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------

# Rainfall Poisson -------------------------------------------------------------
"""
`rainfall_poisson(n, α, λ)`

Generates rainfall using Poisson process.

# Arguments
- `n::Integer`: number of events.
- `α::Float64`: mean rain depth.
- `λ::Float64`: rain event interval (1/d).

If rainfall needs to be similated in timescale smaller than daily, multiply `λ` by dt.

# Example
```
n = 30
α = 0.75
λ = 0.45
rain = rainfall_poisson(n, α, λ)
```
"""
function rainfall_poisson(n, α, λ)
    n = Int(n)
    rain = zeros(n)
    for i in 1:n
        r1 = rand(Uniform(0, 1))
        if r1 < λ
            r2 = rand(Uniform(0, 1))
            rain[i] = -α * log((1 - r2))
        else
            rain[i] = 0
        end
    end
    return rain
end


# Leakage ----------------------------------------------------------------------
"""
`soil_leakage(s, ks, β)`

Compute the soil leakage.

# Arguments
- `s::Float64`: soil moisture.
- `ks::Float64`: saturated hydraulic conductivity.
- `β::Float64`: exponent of hydraulic conductivity curve.
"""
function soil_leakage(s, ks, β)
    lk = ks * s^β
    return lk
end


# Transpiration ----------------------------------------------------------------
"""
`plant_transpiration(s, sw, sstar, ct, eto)`

Compute plant transpiration.

# Arguments
- `s::Float64`: soil moisture.
- `sw::Float64`: soil moisture at wilting point.
- `sstar::Float64`: soil moisture below the field capacity.
- `ct::Float64`: baseline transpiration factor.
- `eto::Float64`: reference evapotranspiration.
"""
function plant_transpiration(s, sw, sstar, ct, eto)
    # Compute soil moisture factor
    if s <= sw
        fs = 0
    elseif s > sw && s <= sstar
        fs = (s - sw)/(sstar - sw)
    else
        fs = 1
    end
    # Compute transpiration
    tp = fs * ct * eto
    return tp
end


# Soil evaporation -------------------------------------------------------------
"""
`soil_evaporation(s, sw, sstar, ct, eto)`

Compute soil evaporation.

# Arguments
- `s::Float64`: soil moisture.
- `sh::Float64`: soil moisture at wilting point.
- `ce::Float64`: baseline evaporation coefficient.
- `eto::Float64`: reference evapotranspiration.
"""
function soil_evaporation(s, sh, ce, eto)
    # Compute soil moisture factor
    if s <= sh
        fs = 0
    else
        fs = (s - sh)/(1 - sh)
    end
    ev = fs * ce * eto
    return ev
end


# Canopy interception ----------------------------------------------------------
"""
`canopy_interception(rain, cws, dt)`

Computes the amount of rain intercepted by the canopy.

# Arguments
- `rain::Float64`: rainfall.
- `cws::Float64`: canopy water storage capacity.
- `dt::Float64`: time interval.

Considering that `cws` is expressed in cm/day, but the rainfall might
be in a smaller timescale, `dt` will adjust `cws` to compute the interception
properly.
"""
function canopy_interception(rain, cws, dt)
    if rain <= (cws * dt)
        ci = rain
    else
        ci = (cws * dt)
    end
    return ci
end


# Cover fraction ---------------------------------------------------------------
"""
`cover_fraction(cf, t, tr, cmax, r, γ, tsen)`

Compute crop cover fraction.

# Arguments
- `cf::Float64`: crop cover fraction.
- `t::Float64`: time [days].
- `tr::Float64`: plant transpiration.
- `cmax::Float64`: maximum crop cover fraction.
- `r::Float64`: parameter of cover fraction model.
- `γ::Float64`: slope of increase senescence after tsen.
- `tsen::Float64`: senescence time.

If plant transpiration is not available, use ETo.
"""
function cover_fraction(cf, t, tr, cmax, r, γ, tsen)
    if (t - tsen) < 0
        hvs = 0
    else
        hvs = 1
    end
    dtc = (r * tr * cf * (1 - (cf/cmax))) - ((γ * (t - tsen) * hvs) * cf^2)
    return dtc
end


# Mulch leakage ----------------------------------------------------------------
"""
`mulch_leakage(ϑ, ϑh, km, g)`

Compute leakage from the mulch layer.

# Arguments
- `ϑ::Float64`: crop cover fraction.
- `ϑh::Float64`: threshold that stops mulch leakage and evaporation.
- `km::Float64`: percolation rate.
- `g::Float64`: exponent of mulching percolation.
"""
function mulch_leakage(ϑ, ϑh, km, g)
    if ϑ <= ϑh
        lm = 0
    else
        lm = km * ϑ^g
    end
    return lm
end


# Mulch evaporation ------------------------------------------------------------
"""
`mulch_evaporation(ϑ, ϑh, ce, eto)`

Compute crop cover fraction.

# Arguments
- `ϑ::Float64`: crop cover fraction.
- `ϑh::Float64`: threshold that stops mulch leakage and evaporation.
- `ce::Float64`: baseline evaporation coefficient.
- `eto::Float64`: reference evapotranspiration.
"""
function mulch_evaporation(ϑ, ϑh, ce, eto)
    if ϑ <= ϑh
        fs = 0
    else
        fs = (ϑ - ϑh) / (1 - ϑh)
    end
    em = fs * ce * eto
    return em
end


# Mulch evaporation reduction factor -------------------------------------------
"""
`mulch_res_zm(μ, zm)`

Computes soil evaporation resistance factor due to the mulching layer thickness.

# Arguments
- `μ::Float64`: an empirical constant to compute resistance factor.
- `zm::Float64`: mulching layer thickness [cm].

Reference: Jalota, S., Prihar, S., 1990. Effect of straw mulch on evaporation
reduction in relation to rates of mulching and evaporativity. J. Indian Soc.
Soil Sci. 38, 728–730.

Suggested value for `μ` = 0.11.
"""
function mulch_res_zm(μ, zm)
    etom = exp.(-μ * zm)
    return etom
end


# Soil water balance -----------------------------------------------------------
"""
`soil_water_balance(rain, eto, cf, s, sh, sw, sstar, ks, β, n, zr, rth, ce, ct, dt)`

Computes the soil water balance under no mulching conditions.

# Arguments
- `rain::Float64`: rainfall.
- `eto::Float64`: reference evapotranspiration.
- `cf::Float64`: crop cover fraction.
- `s::Float64`: soil moisture.
- `sh::Float64`: soil moisture at wilting point.
- `sw::Float64`: soil moisture at wilting point.
- `sstar::Float64`: soil moisture below the field capacity.
- `ks::Float64`: saturated hydraulic conductivity.
- `β::Float64`: exponent of hydraulic conductivity curve.
- `n::Float64`: soil porosity.
- `zr::Float64`: root depth.
- `cws::Float64`: canopy water storage capacity.
- `ce::Float64`: baseline evaporation coefficient.
- `ct::Float64`: baseline transpiration factor.
- `dt::Float64`: time interval.
"""
function soil_water_balance(rain, eto, cf, s, sh, sw, sstar, ks, β, n, zr, cws, ce, ct, dt)
    # Soil water storage capacity
    nzr = n * zr

    # Throughfall
    ci = canopy_interception(rain, cws, dt) * cf
    tf = rain - ci

    s = s + tf / nzr

    # Check if there is runoff
    if s > 1
        q = (s - 1) * nzr
        s = 1
    else
        q = 0
    end

    # Compute plant transpiration
    tr = plant_transpiration(s, sw, sstar, ct, eto) * cf * dt
    if s >= sstar
        trns = tr
    else
        trns = 0
    end
    s = s - tr / nzr

    # Compute soil evaporation
    ev = soil_evaporation(s, sh, ce, eto) * (1 - cf) * dt
    s = s - ev / nzr

    # Compute leakage
    lk = soil_leakage(s, ks, β) * dt
    s = s - lk / nzr

    # Write outputs
    y = (s = s, rain = rain, ci = ci, q = q, tr = tr, trns = trns, ev = ev, lk = lk)
    return y
end


# Soil water balance with mulching ---------------------------------------------

"""
`soil_mulch_water_balance(rain, eto, cf, ϑ, s, sh, sw, sstar, ks, β, n, zr, cws,
                          ce, ct, μ, ϕ, zm, ϑh, km, g, dt)`

Computes the mulching and soil water balance.

# Arguments
- `rain::Float64`: rainfall.
- `eto::Float64`: reference evapotranspiration.
- `cf::Float64`: crop cover fraction.
- `ϑ::Float64`: mulching moisture.
- `s::Float64`: soil moisture.
- `sh::Float64`: soil moisture at wilting point.
- `sw::Float64`: soil moisture at wilting point.
- `sstar::Float64`: soil moisture below the field capacity.
- `ks::Float64`: saturated hydraulic conductivity.
- `β::Float64`: exponent of hydraulic conductivity curve.
- `n::Float64`: soil porosity.
- `zr::Float64`: root depth.
- `cws::Float64`: canopy water storage capacity.
- `ce::Float64`: baseline evaporation coefficient.
- `ct::Float64`: baseline transpiration factor.
- `ϕ::Float64`: mulching porosity.
- `ϑh::Float64`: threshold that stops mulch leakage and evaporation.
- `km::Float64`: mulching percolation rate.
- `g::Float64`: exponent of mulching percolation.
- `dt::Float64`: time interval.
"""
function soil_mulch_water_balance(rain, eto, cf, ϑ, s, sh, sw, sstar, ks, β, n,
                                  zr, cws, ce, ct, μ, ϕ, zm, ϑh, km, g, dt)
    # Soil and mulch water storage capacity
    pzm = ϕ * zm
    nzr = n * zr

    # Throughfall
    ci = canopy_interception(rain, cws, dt) * cf
    tf = rain - ci

    ϑ = ϑ + tf / pzm

    if ϑ > 1
        q = (ϑ - 1) * pzm
        ϑ = 1
    else
        q = 0
    end

    # Mulching leakage
    lm = mulch_leakage(ϑ, ϑh, km, g) * dt
    ϑ = ϑ - lm / pzm

    # Mulching evaporation
    em = mulch_evaporation(ϑ, ϑh, ce, eto) * dt
    ϑ = ϑ - em / pzm

    # Soil layer
    s = s + lm / nzr

    if s > 1
        q = q + (s - 1) * nzr
        s = 1
    else
        q = q + 0
    end

    # Compute plant transpiration
    tr = plant_transpiration(s, sw, sstar, ct, eto) * cf * dt
    if s >= sstar
        trns = tr
    else
        trns = 0
    end
    s = s - tr / nzr

    # Compute soil evaporation
    etom = eto * mulch_res_zm(μ, zm)
    evr = soil_evaporation(s, sh, ce, etom) * (1 - cf) * dt - em
    if evr <= 0
        ev = 0
    else
        ev = evr
    end

    s = s - ev / nzr

    # Compute leakage
    lk = soil_leakage(s, ks, β) * dt
    s = s - lk / nzr

    # Write outputs
    y = (s = s, rain = rain, ci = ci, q = q, tr = tr, trns = trns, ev = ev,
         lk = lk, ϑ = ϑ, em = em, lm = lm)
    return y
end


# Convert results to daily timescale -------------------------------------------
"""
`df_dt_day(df)`

Convert model output from `dt` to daily timescale.

# Arguments
-`df::DataFrame`: a data frame.
"""
function df_dt_day(df)
    df1 = combine(groupby(df, :Days), :Rain=>sum, :CI=>sum, :Q=>sum, :s=>mean,
                  :Lk=>sum, :Es=>sum, :Tr=>sum, :Trns=>sum)
    # Check if there is mulching layer
    if "ϑ" in names(df)
        df2 = combine(groupby(df, :Days), :ϑ=>mean, :Lm=>sum, :Em=>sum)
        df1 = rightjoin(df1, df2, on=:Days)
    end

    # Check if there is crop cover fraction
    if "CF" in names(df)
        df3 = combine(groupby(df, :Days), :CF=>last)
        df1 = rightjoin(df1, df3, on=:Days)
    end

    # Get columns name and remove the operation
    cn = names(df1)
    cn2 = []
    for i in 1:length(cn)
        push!(cn2, split(cn[i], "_")[1])
    end

    rename!(df1, Symbol.(cn2))
    return df1
end


# ------------------------------------------------------------------------------
# Solutions
# ------------------------------------------------------------------------------

# Solution crop cover fraction -------------------------------------------------
"""
`cover_fraction_sol(cf, t, tr, cmax, r, γ, tsen)`

Compute numerical solution for the crop cover fraction.

# Arguments
- `c::Float64`: crop cover fraction.
- `t::Float64`: time [days].
- `tr::Float64`: plant transpiration.
- `cmax::Float64`: maximum crop cover fraction.
- `r::Float64`: parameter of cover fraction model.
- `γ::Float64`: slope of increase senescence after tsen.
- `tsen::Float64`: senescence time.

If plant transpiration is not available, use ETo.
"""
function cover_fraction_sol(cf, t, tr, cmax, r, γ, tsen)
    # Solver cover model for the growing season
    cday = zeros(length(eto))
    cday[1] = cf
    for i in 2:length(cday)
        cday[i] = cday[i-1] + cover_fraction(cday[i-1], t[i], tr[i], cmax, r, γ, tsen) * 1
    end
    return cday
end


# Solve soil water balance -----------------------------------------------------
"""
`sol_swb(rain, eto, cf, s, p)`

Solves the soil water balance function.

# Arguments
- `p::NamedTuple`: p is a named tuple with the parameters mostly considered as
constant.

For the parameters descriptions, please check the help of
`soil_water_balance` function.

"""
function sol_swb(rain, eto, cf, s, p)
    # Unzip parameter
    sh = p.sh
    sw = p.sw
    sstar = p.sstar
    ks = p.ks
    β = p.β
    n = p.n
    zr = p.zr
    cws = p.cws
    ce = p.ce
    ct = p.ct
    dt = p.dt

    nr = length(rain) + 1
    sres = zeros(nr)
    cires = zeros(nr)
    qres = zeros(nr)
    eres = zeros(nr)
    tres = zeros(nr)
    trnsres = zeros(nr)
    lres = zeros(nr)

    # Initial condition for soil moisture
    sres[1] = s

    # Create vector for time
    ndays = Int((nr - 1) * dt)
    days = [1:ndays;]
    days = repeat(days, inner = Int(1 / dt))
    days = vcat(0, days)
    daysc = [(1+dt):dt:(ndays+1);] .- dt

    # Adjust ETo
    if length(eto) == 1
        eto = repeat([eto], length(rain))
    else
        eto = eto
    end

    # Adjust cover fraction
    if length(cf) == 1
        cf = repeat([cf], length(rain))
    else
        cf = cf
    end

    for i in 1:length(rain)
        swb = soil_water_balance(rain[i], eto[i], cf[i], sres[i],
                                 sh, sw, sstar, ks, β,
                                 n, zr, cws, ce, ct, dt)
        sres[i+1] = swb.s
        cires[i+1] = swb.ci
        qres[i+1] = swb.q
        eres[i+1] = swb.ev
        tres[i+1] = swb.tr
        trnsres[i+1] = swb.trns
        lres[i+1] = swb.lk
    end
    df = DataFrame(Days = days, DaysC = vcat(0, daysc), Rain = vcat(0, rain),
                   CI = cires, Q = qres, s = sres, Lk = lres, Es = eres,
                   Tr = tres, Trns = trnsres)
    return df

end


# Solve soil water balance coupled with crop -----------------------------------
"""
`sol_swb_crop(rain, eto, cf, s, p)`

Solves the soil water balance coupled with the crop cover fraction.

# Arguments
- `p::NamedTuple`: p is a named tuple with the parameters mostly considered as
constant.

For the parameters descriptions, please check the help of
`soil_water_balance` function.

"""
function sol_swb_crop(rain, eto, cf, s, p)
    # Unzip parameter
    sh = p.sh
    sw = p.sw
    sstar = p.sstar
    ks = p.ks
    β = p.β
    n = p.n
    zr = p.zr
    cws = p.cws
    ce = p.ce
    ct = p.ct
    cmax = p.cmax
    r = p.r
    γ = p.γ
    tsen = p.tsen
    dt = p.dt

    nr = length(rain) + 1
    sres = zeros(nr)
    cires = zeros(nr)
    qres = zeros(nr)
    eres = zeros(nr)
    tres = zeros(nr)
    trnsres = zeros(nr)
    lres = zeros(nr)
    cfres = zeros(nr)

    # Initial condition for soil moisture
    sres[1] = s
    cfres[1] = cf

    # Create vector for time
    ndays = Int((nr - 1) * dt)
    days = [1:ndays;]
    days = repeat(days, inner = Int(1 / dt))
    days = vcat(0, days)
    daysc = [(1+dt):dt:(ndays+1);] .- dt

    # Adjust ETo
    if length(eto) == 1
        eto = repeat([eto], length(rain))
    else
        eto = eto
    end

    for i in 1:length(rain)
        swb = soil_water_balance(rain[i], eto[i], cfres[i], sres[i],
                                 sh, sw, sstar, ks, β,
                                 n, zr, cws, ce, ct, dt)
        sres[i+1] = swb.s
        cires[i+1] = swb.ci
        qres[i+1] = swb.q
        eres[i+1] = swb.ev
        tres[i+1] = swb.tr
        trnsres[i+1] = swb.trns
        lres[i+1] = swb.lk
        cfres[i+1] = cfres[i] + cover_fraction(cfres[i], daysc[i], swb.tr,
                                               cmax, r, γ, tsen)
    end
    df = DataFrame(Days = days, DaysC = vcat(0, daysc), Rain = vcat(0, rain),
                   CI = cires, Q = qres, s = sres, Lk = lres, Es = eres,
                   Tr = tres, Trns = trnsres, CF = cfres)
    return df

end


# Solve mulching soil water balance --------------------------------------------
"""
`sol_mswb(rain, eto, cf, ϑ, s, p)`

Solves the soil water balance function.

# Arguments
- `p::NamedTuple`: p is a named tuple with the parameters mostly considered as
constant.

For the parameters descriptions, please check the help of
`soil_mulch_water_balance` function.

"""
function sol_mswb(rain, eto, cf, ϑ, s, p)

    # Unzip parameter
    sh = p.sh
    sw = p.sw
    sstar = p.sstar
    ks = p.ks
    β = p.β
    n = p.n
    zr = p.zr
    cws = p.cws
    ce = p.ce
    ct = p.ct
    μ = p.μ
    ϕ = p.ϕ
    zm = p.zm
    ϑh = p.ϑh
    km = p.km
    g = p.g
    dt = p.dt

    # Vectors to receive outputs
    nr = length(rain) + 1
    sres = zeros(nr)
    ϑres = zeros(nr)
    cires = zeros(nr)
    qres = zeros(nr)
    eres = zeros(nr)
    tres = zeros(nr)
    trnsres = zeros(nr)
    lres = zeros(nr)
    lmres = zeros(nr)
    emres = zeros(nr)

    # Initial condition for soil moisture
    sres[1] = s
    ϑres[1] = ϑ

    # Create vector for time
    ndays = Int((nr - 1) * dt)
    days = [1:ndays;]
    days = repeat(days, inner = Int(1 / dt))
    days = vcat(0, days)
    daysc = [(1+dt):dt:(ndays+1);] .- dt

    # Adjust ETo
    if length(eto) == 1
        eto = repeat([eto], length(rain))
    else
        eto = eto
    end

    # Adjust cover fraction
    if length(cf) == 1
        cf = repeat([cf], length(rain))
    else
        cf = cf
    end

    for i in 1:length(rain)
        swb = soil_mulch_water_balance(rain[i], eto[i], cf[i], ϑres[i], sres[i],
                                       sh, sw, sstar, ks, β, n, zr, cws, ce, ct,
                                       μ, ϕ, zm, ϑh, km, g, dt)
        sres[i+1] = swb.s
        ϑres[i+1] = swb.ϑ
        cires[i+1] = swb.ci
        qres[i+1] = swb.q
        eres[i+1] = swb.ev
        tres[i+1] = swb.tr
        trnsres[i+1] = swb.trns
        lres[i+1] = swb.lk
        emres[i+1] = swb.em
        lmres[i+1] = swb.lm
    end
    df = DataFrame(Days = days, DaysC = vcat(0, daysc), Rain = vcat(0, rain),
                   CI = cires, Q = qres, s = sres, Lk = lres, Es = eres,
                   Tr = tres, Trns = trnsres, ϑ = ϑres, Lm = lmres, Em = emres)
    return df

end


# Solve mulching soil water balance coupled with crop --------------------------
"""
`sol_mswb_crop(rain, eto, cf, ϑ, s, p)`

Solves the mulchg and soil water balance coupled with the crop cover fraction.

# Arguments
- `p::NamedTuple`: p is a named tuple with the parameters mostly considered as
constant.

For the parameters descriptions, please check the help of
`soil_mulch_water_balance` function.

"""
function sol_mswb_crop(rain, eto, cf, ϑ, s, p)

    # Unzip parameter
    sh = p.sh
    sw = p.sw
    sstar = p.sstar
    ks = p.ks
    β = p.β
    n = p.n
    zr = p.zr
    cws = p.cws
    ce = p.ce
    ct = p.ct
    μ = p.μ
    ϕ = p.ϕ
    zm = p.zm
    ϑh = p.ϑh
    km = p.km
    g = p.g
    cmax = p.cmax
    r = p.r
    γ = p.γ
    tsen = p.tsen
    dt = p.dt

    # Vectors to receive outputs
    nr = length(rain) + 1
    sres = zeros(nr)
    ϑres = zeros(nr)
    cires = zeros(nr)
    qres = zeros(nr)
    eres = zeros(nr)
    tres = zeros(nr)
    trnsres = zeros(nr)
    lres = zeros(nr)
    lmres = zeros(nr)
    emres = zeros(nr)
    cfres = zeros(nr)

    # Initial condition for soil moisture
    sres[1] = s
    ϑres[1] = ϑ
    cfres[1] = cf

    # Create vector for time
    ndays = Int((nr - 1) * dt)
    days = [1:ndays;]
    days = repeat(days, inner = Int(1 / dt))
    days = vcat(0, days)
    daysc = [(1+dt):dt:(ndays+1);] .- dt

    # Adjust ETo
    if length(eto) == 1
        eto = repeat([eto], length(rain))
    else
        eto = eto
    end

    for i in 1:length(rain)
        swb = soil_mulch_water_balance(rain[i], eto[i], cfres[i], ϑres[i], sres[i],
                                       sh, sw, sstar, ks, β, n, zr, cws, ce, ct,
                                       μ, ϕ, zm, ϑh, km, g, dt)
        sres[i+1] = swb.s
        ϑres[i+1] = swb.ϑ
        cires[i+1] = swb.ci
        qres[i+1] = swb.q
        eres[i+1] = swb.ev
        tres[i+1] = swb.tr
        trnsres[i+1] = swb.trns
        lres[i+1] = swb.lk
        emres[i+1] = swb.em
        lmres[i+1] = swb.lm
        cfres[i+1] = cfres[i] + cover_fraction(cfres[i], daysc[i], swb.tr,
                                               cmax, r, γ, tsen)
    end
    df = DataFrame(Days = days, DaysC = vcat(0, daysc), Rain = vcat(0, rain),
                   CI = cires, Q = qres, s = sres, Lk = lres, Es = eres,
                   Tr = tres, Trns = trnsres, ϑ = ϑres, Lm = lmres, Em = emres,
                   CF = cfres)
    return df

end


end
