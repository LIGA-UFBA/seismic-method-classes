import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import animation
from IPython.display import display, HTML


__all__ = [
    'plot_animated_wavefield',
    'plot_wavefield_snaps',
    'plot_shotrecord'
]

def plot_animated_wavefield(u, rec, geometry, src_pos, model, title=None, rec_pos=None, scale_factor=1/10, wave_scale_factor=1, 
                            intf_positions=None):
    """
    Reproduz uma animação do campo de ondas se propagando.

    Parameters
    ----------
    u : np.ndarray
        Dados do campo de ondas (propriedade `.data`).
    rec : np.ndarray
        Dados dos receptores (propriedade `.data`).
    geometry : AquisitionGeometry
        Geometria usada na modelagem.
    src_pos : np.ndarray
        Posição da fonte.
    model : ElasticModel or ISOSeismicModel
    title : str
        Título da plotagem.
    rec_pos : np.ndarray
        Posições dos receptores.
    scale_factor : float
        Exagero da amplitude do sismograma.
    wave_scale_factor : float
        Exagero da amplitude do campo de ondas.
    intf_positions: List[List[(float, float)]]
        Posições das interfaces
    """
    shape_pad = np.array(model.shape)
    origin_pad = model.origin
    extent_pad = tuple([s*(n-1) for s, n in zip(model.spacing, shape_pad)])

    t0 = geometry.t0
    tn = geometry.tn
    time_factor = geometry.nt/u.shape[0]
    
    # Note: flip direction of second dimension to make the plot positive downwards
    plt_extent = [origin_pad[0], origin_pad[0] + extent_pad[0],
                origin_pad[1] + extent_pad[1], origin_pad[1]]
    plt_extent2 = [model.origin[0], model.origin[0] + model.domain_size[0],
              tn, t0]
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 10), sharex=False)
    ax = axes[1]
    sr_ax = axes[0]

    snapshots = np.linspace(0, u.shape[0]-1, num=u.shape[0]-1, dtype=np.int32)

    udata = u
    amax = np.max(np.abs(udata[-2, :, :])) * wave_scale_factor
    
    snapshot = snapshots[0]
    wave_plot = ax.imshow(udata[snapshot, model.nbl:-model.nbl, model.nbl:-model.nbl].T,
                            cmap="seismic", vmin=-amax,
                vmax=+amax, extent=plt_extent)
    # wave_plot.autoscale()
    sr_plot = sr_ax.imshow(rec*0,
                            cmap="gray", vmin=-np.max(rec) * scale_factor,
                vmax=+np.max(rec) * scale_factor, extent=plt_extent2)
    # sr_plot.autoscale()

    nbl_color = 'xkcd:dark blue'

    if intf_positions:
        for intf in intf_positions:

            ax.plot(
                intf[0], intf[1], color=nbl_color, ls='--'
            )
    
    ax.plot(src_pos[0], src_pos[1], 'red', linestyle='None', marker='*',
            markersize=8, label="Source")
    
    if type(rec_pos) != type(None):
        ax.plot(rec_pos[0], rec_pos[1], 'green', linestyle='None', marker='v',
                markersize=8, label="Receiver")
        
    time_text = ax.text(
        plt_extent[0] + (plt_extent[1] - plt_extent[0])*0.05, 
        plt_extent[2] + (plt_extent[3] - plt_extent[2])*0.05
    , 't=...', color=nbl_color)
    # ax.set_title(title)

    ax.grid()
    ax.tick_params('both', length=4, width=0.5, which='major', labelsize=11)

    ax.set_xlabel("x (m)", fontsize=12)
    ax.set_ylabel("z (m)", fontsize=12)
    
    sr_ax.set_xlabel("x (m)", fontsize=12)
    sr_ax.set_ylabel("Time (ms)", fontsize=12)
    
    if title:
        fig.suptitle(title)

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(wave_plot, cax=cax)

    def update(snapshot):
        t = (snapshot*geometry.dt*time_factor + geometry.t0)
        wave_plot.set_array(udata[snapshot, model.nbl:-model.nbl, model.nbl:-model.nbl].T)
        mask = np.vstack([(geometry.time_axis.time_values <= t)]*rec.shape[1]).T
        sr_plot.set_array(rec * mask)
        time_text.set_text("t=%.2fms" % t)

    anim = animation.FuncAnimation(fig, update, frames=snapshots, blit=False)

    plt.close(fig)
    display(HTML(anim.to_jshtml()))



def plot_wavefield_snaps(u, geometry, src_pos, model, title=None, rec_pos=None, t0=0, tn=-1, wave_scale_factor=1, 
                            intf_positions=None, annotation_callback=None, figsize=None):
    """
    Exibe o campo de onda em 4 instantes distribuídos em um período determinado por `t0` e `tn`

    Parameters
    ----------
    u : np.ndarray
        Dados do campo de ondas (propriedade `.data`).
    geometry : AquisitionGeometry
        Geometria usada na modelagem.
    src_pos : np.ndarray
        Posição da fonte.
    model : SeismicModel
    title : str
        Título da plotagem.
    rec_pos : np.ndarray
        Posições dos receptores.
    t0 : float
        Tempo inicial desejado.
    tn : float
        Tempo final desejado.
    wave_scale_factor : float
        Exagero da amplitude do campo de ondas.
    intf_positions: List[List[(float, float)]]
        Posições das interfaces
    """
    shape_pad = np.array(model.shape)
    origin_pad = model.origin
    extent_pad = tuple([s*(n-1) for s, n in zip(model.spacing, shape_pad)])

    # Note: flip direction of second dimension to make the plot positive downwards
    plt_extent = [origin_pad[0], origin_pad[0] + extent_pad[0],
                origin_pad[1] + extent_pad[1], origin_pad[1]]
    
    nrows = 2
    ncols = 2

    fig, axes = plt.subplots(nrows, ncols, figsize=(10, 12)) if figsize == None else plt.subplots(nrows, ncols, figsize=figsize)
    fig.suptitle(title, size=20, y=1.0)
    isfirst = True
    last_ax = None
    last_plot = None

    t0 = t0 - geometry.t0
    tn = tn - geometry.t0
    time_factor = geometry.nt/u.shape[0]
    dt = geometry.dt * time_factor

    snapshots = None
    if tn < 0:
        idx_min = max(int(np.round(t0/(dt))), 0)
        snapshots = np.linspace(idx_min, u.shape[0] - 1, num=4, dtype=np.int32)
    else:
        idx_min = max(int(np.round(t0/dt)), 0)
        idx_max = min(int(np.round(tn/dt)), u.shape[0] - 1)
        snapshots = np.linspace(idx_min, idx_max, num=4, dtype=np.int32)


    udata = u
    amax = np.max(np.abs(udata[-1, :, :])) * wave_scale_factor

    nbl_color = 'xkcd:dark blue'

    for r in range(nrows):
        for c in range(ncols):
            
        # for count, ax in enumerate(ax_row.ravel()):
            count = r * nrows + c
            ax = axes[r][c]

            snapshot = snapshots[count]
            last_plot = ax.imshow(udata[snapshot, model.nbl:-model.nbl, model.nbl:-model.nbl].T,
                                cmap="seismic", vmin=-amax,
                    vmax=+amax, extent=plt_extent)
            ax.plot(src_pos[0], src_pos[1], 'red', linestyle='None', marker='*',
                    markersize=8, label="Source")
            
            if type(rec_pos) != type(None):
                ax.plot(rec_pos[0], rec_pos[1], 'green', linestyle='None', marker='v',
                        markersize=8, label="Receiver")

            if intf_positions:
                for intf in intf_positions:
                    ax.plot(
                        intf[0], intf[1], color=nbl_color, ls='--'
                    )
                
            ax.grid()
            ax.tick_params('both', length=4, width=0.5, which='major', labelsize=11)

            ax.set_title("t=%.2fms" % (snapshot*dt + geometry.t0), fontsize=12)
            ax.set_xlabel("x (m)", fontsize=12)
            isfirst and ax.set_ylabel("z (m)", fontsize=12)
            isfirst = False
            last_ax = ax

        isfirst = True
        # divider = make_axes_locatable(last_ax)
    # divider = make_axes_locatable(axes)
    # cax = divider.append_axes("right", size="5%", pad=0.05)
    # plt.colorbar(last_plot, cax=cax)
    fig.colorbar(last_plot, ax=axes, fraction=.05)
    
    if annotation_callback:
        annotation_callback(fig, axes)



def plot_shotrecord2(rec, model, t0, tn, colorbar=True, title=None, scale_factor = 1/10, annotation_callback=None):
    scale = np.max(rec) * scale_factor
    extent = [model.origin[0]*1e-3, model.origin[0]*1e-3 + 1e-3*model.domain_size[0],
              1e-3*tn, t0*1e-3]

    plot = plt.imshow(rec, vmin=-scale, vmax=scale, cmap='gray', extent=extent)
    plt.xlabel('X position (km)')
    plt.ylabel('Time (s)')
    if title:
        plt.title(title)

    plt.grid()
    if annotation_callback:
        annotation_callback()
    
    # Create aligned colorbar on the right
    if colorbar:
        ax = plt.gca()
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)
        plt.colorbar(plot, cax=cax)
    plt.show()


def plot_shotrecord3(rec, model, t0, tn, colorbar=True, title=None, scale_factor = 1/10, annotation_points=None):
    def annotation_callback():
        i = 1
        def mark_point(x, y):
            nonlocal i
            plt.annotate(i, (x,y), color='red', size=10, 
                        bbox={ 'boxstyle': 'square', 'ec': 'red', 'fc': 'white' })
            i += 1
        
        for p in annotation_points:
            mark_point(p[0], p[1])

    if annotation_points:
        return plot_shotrecord2(rec, model, t0, tn, colorbar, title, scale_factor, annotation_callback=annotation_callback)
    return plot_shotrecord2(rec, model, t0, tn, colorbar, title, scale_factor)

plot_shotrecord = plot_shotrecord3