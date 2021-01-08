
from mtuq.graphics.beachball import\
    plot_beachball, misfit_vs_depth

from mtuq.graphics.uq.lune import\
    plot_misfit_lune, plot_likelihood_lune, plot_marginal_lune, plot_misfit_mt_lune, plot_magnitude_lune

from mtuq.graphics.uq.double_couple import\
    plot_misfit_dc, plot_likelihood_dc, plot_marginal_dc

from mtuq.graphics.uq.vw import\
    plot_misfit_vw, plot_likelihood_vw, plot_marginal_vw

from mtuq.graphics.uq.force import\
    plot_misfit_force, plot_likelihood_force, plot_marginal_force, plot_force_amplitude

from mtuq.graphics.uq.origin_depth import\
    plot_misfit_depth, plot_likelihood_depth, plot_marginal_depth

from mtuq.graphics.uq.origin_xy import\
    plot_misfit_xy, plot_mt_xy

from mtuq.graphics.waveform import\
    plot_waveforms1, plot_waveforms2, plot_data_greens


# use Helvetica if available
for fontname in ['Helvetica', 'Arial']:
    try:
        from matplotlib.font_manager import find_font
        find_font(fontname, fallback_to_default=False)
        matplotlib.rcParams['font.sans-serif'] = fontname
        matplotlib.rcParams['font.family'] = "sans-serif"
        break
    except:
        continue
