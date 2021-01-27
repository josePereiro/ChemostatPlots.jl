const EP_COLOR = :red
const FBA_COLOR = :blue
const HR_COLOR = :black

function plot_marginal!(p, μ::Real, σ::Real, lb::Real, ub::Real;
    av = nothing, h = 1, lw = 5, color = :black,
    label = "", normalize = false, kwargs...
)
    
    sσ = sqrt(σ)
    tN = Truncated(Normal(μ, sσ), lb, ub) 
    global_max_ = isnothing(av) ? clamp(μ, lb, ub) : av
    margin_ = abs(ub - lb) * 0.1
    if sσ == 0.0 || isinf(pdf(tN, global_max_))
        # delta
        plot!(p, [lb - margin_, ub + margin_], [0.0, 0.0]; label = "")
        plot!(p, [global_max_, global_max_], [0.0, h]; 
            lw, label, color, kwargs...
        )
    else
        # normal
        xs = range(lb, ub; length = 1000)
        ys = [pdf(tN, x) for x in xs]
        ys .= !normalize ? ys ./ maximum(ys) : ys
        plot!(p, xs, ys;
            xlim = [lb - margin_, ub + margin_],
            lw, label, color, kwargs...
        )
    end
    return plot!(p; kwargs...)
end

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
    σ_ = ChU.σ(metnet, out, ider);
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
