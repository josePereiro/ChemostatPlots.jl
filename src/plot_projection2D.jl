function plot_projection2D!(p, proj; l = 10, spkwargs...)

    def_spkwargs = (;label = "", alpha = 0.5)

    flx1s = sort(collect(keys(proj)))
    rflx1s = reverse(flx1s)
    flx2lbs = [first(proj[flx1]) for flx1 in flx1s]
    flx2ubs = [last(proj[flx1]) for flx1 in rflx1s]

    s = Shape([flx1s; rflx1s], [flx2lbs; flx2ubs])
    plot!(p, [s]; label = "", def_spkwargs..., spkwargs...)

    return p
end

## ----------------------------------------------------------------------------
function plot_projection2D!(p::AbstractPlot, model::ChU.MetNet, 
        ider1::ChU.IDER_TYPE, ider2::ChU.IDER_TYPE; 
        l = 10, spkwargs...
    )

    proj = ChLP.projection2D(model, ider1, ider2; l)
    plot_projection2D!(p, proj, l, spkwargs...)

end