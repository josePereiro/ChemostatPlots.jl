function plot_projection2D!(p, proj; l = 10, spkwargs...)

    def_spkwargs = (;label = "", lw = 8, alpha = 0.2, color = :black)
    for (flx1, (flx2_lb, flx2_ub)) in proj
        plot!(p, [flx1, flx1], [flx2_lb, flx2_ub]; def_spkwargs..., spkwargs...)
    end
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