const EP_COLOR = :red
const FBA_COLOR = :blue
const HR_COLOR = :black

## ----------------------------------------------------------------------------------
ϕ(x, μ, σ) = inv(σ*sqrt(2π))*exp(-((x - μ)^2)/(2σ^2))
_unbig(x) = round(Float64(x); digits = 15)
function _trunc_samples(μ, σ, lb, ub; xbins = 1000)
    
    local_max = clamp(μ, lb, ub)
    max_range = max(abs(local_max - lb), abs(local_max - ub))
    log_range = range(log10(max_range), log10(σ/xbins), length = div(xbins, 2))
    xs = [
        local_max .- 10.0.^(log_range);
        local_max;
        local_max .+ 10.0.^(reverse(log_range))
    ] 
    Txs = [big(x) for x in xs if lb <= x <= ub]
    Tpdf = ϕ.(Txs, μ, σ)
    Z = sum(@view(Tpdf[1:end - 1]) .* diff(Txs))
    Tpdf .= Tpdf ./ Z
    _unbig.(Txs), _unbig.(Tpdf)
end

## ----------------------------------------------------------------------------------
function plot_marginal!(p, μ::Real, σ::Real, lb::Real, ub::Real;
        av = nothing, h = 1, lw = 5, color = :black, xbins = 1000,
        label = "", normalize = false, kwargs...
    )
    σ = σ < 1e-25 ? 1e-25 : σ

    xs, ys = _trunc_samples(μ, σ, lb, ub; xbins)
    nans = findfirst(isnan, ys)
    if !isnothing(nans)
        
        global_max_ = isnothing(av) ? clamp(μ, lb, ub) : av
        xs, ys = _trunc_samples(global_max_, σ, lb, ub; xbins)
    end
    
    margin_ = abs(ub - lb) * 0.1
    ys .= !normalize ? ys ./ maximum(ys) : ys
    plot!(p, xs, ys;
        xlim = [lb - margin_, ub + margin_],
        lw, label, color, kwargs...
    )
    return plot!(p; kwargs...)
end

## ----------------------------------------------------------------------------------
# Single out
function plot_marginal!(p, metnet::ChU.MetNet, out::Union{ChU.FBAout, ChU.EPout}, ider; 
        color = out isa ChU.FBAout ? FBA_COLOR : EP_COLOR,
        label = ChU.rxns(metnet, ider), 
        kwargs...
    )
        
    lb_ = ChU.lb(metnet, ider)
    ub_ = ChU.ub(metnet, ider)

    μ_ = ChU.μ(metnet, out, ider);
    av_ = ChU.av(metnet, out, ider);
    σ_ = sqrt(ChU.σ(metnet, out, ider));
    plot_marginal!(p, μ_, σ_, lb_, ub_; av = av_, color, label, kwargs...)
end


# TODO: add normalization option
function plot_marginal!(p, metnet::ChU.MetNet, out::ChU.HRout, ider; 
        h = :ignored, color = HR_COLOR, 
        label = ChU.rxns(metnet, ider), kwargs...)

    lb_ = ChU.lb(metnet, ider)
    ub_ = ChU.ub(metnet, ider)
    margin_ = abs(ub_ - lb_) * 0.1
    hist = hists(metnet, out, ider)
    plot!(p, normalize(hist, mode = :pdf); 
        xaxis = [lb_ - margin_, ub_ + margin_],
        label, color, kwargs...
    )
end

function plot_marginal!(p, metnet::ChU.MetNet, outs::Vector, ider::ChU.IDER_TYPE; kwargs...)
    pdf_maxval_ = ChU.pdf_maxval(metnet, outs, ider)
    pdf_maxval_ = pdf_maxval_ != -1 ? pdf_maxval_ : 1
    for out in outs
        plot_marginal!(p, metnet, out, ider; h = pdf_maxval_ * 1.1, kwargs...)
    end
    return p
end

## ----------------------------------------------------------------------------------
plot_marginal(args...; kwargs...) = plot_marginal!(plot(), args...; kwargs...)