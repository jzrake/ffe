import os


def load_tseries_data(rundir, recache=False, tmin=0.01, tmax=None):
    """Load a pickle of a time series from the disk if it exists, otherwise load the
    data from the text file and create the pickle

    """
    import pickle
    import numpy as np
    import h5py

    pickle_name = '%s/tseries.pk' % rundir

    if os.path.isfile(pickle_name) and not recache:
        print "load from pickle...", rundir
        return pickle.load(open(pickle_name, 'r'))
    else:
        print "no pickle, caching data"

    h5f = h5py.File("%s/analysis.h5" % rundir, 'r')
    ffe_dat = np.loadtxt("%s/ffe.dat"% rundir)

    tk = [ ]
    k_X = { }

    k_X['electric'] = [ ]
    k_X['magnetic'] = [ ]
    k_X['helicity'] = [ ]

    k_E, k_B, k_H = k_X['electric'], k_X['magnetic'], k_X['helicity']

    h5f.close()

    tk = np.array(tk)
    pickle.dump([tk, k_X, ffe_dat], open(pickle_name, 'w'))
    return tk, k_X, ffe_dat



class TimeSeriesPlot(object):


    def __init__(self, figure):
        """Should not modify the figure instance, i.e. do not leave any new axes
        attached to it, as add_subplot would do

        """
        from matplotlib import lines
        from matplotlib import axes

        self.magnetic = lines.Line2D([], [])
        self.electric = lines.Line2D([], [])
        self.helicity = lines.Line2D([], [])
        self.monopole = lines.Line2D([], [])

        self.ax1 = axes.Axes(figure, [0.1, 0.1, 0.8, 0.8])
        self.ax1.add_line(self.magnetic)
        self.ax1.add_line(self.electric)
        self.ax1.add_line(self.helicity)
        self.ax1.add_line(self.monopole)
        self.configure()


    def get_axes(self):
        return [self.ax1]


    def set_data_source(self, data_source):
        """Heavier re-plot; may re-load data from the disk

        """
        rundir = data_source.get_rundir()
        if rundir is None: return
        tk, k_X, ffe_dat = load_tseries_data(rundir)

        t  = ffe_dat[:,0]
        Ue = ffe_dat[:,1]
        Ub = ffe_dat[:,2]
        Hm = ffe_dat[:,3]
        Bm = ffe_dat[:,4]

        #Ub /= Ub[0]

        self.magnetic.set_data(t, Ub)
        self.electric.set_data(t, Ue)
        self.helicity.set_data(t, abs(Hm))
        self.monopole.set_data(t, abs(Bm))
        self.ax1.relim(visible_only=True)
        self.ax1.autoscale_view()


    def configure(self, **configuration):
        """Light re-plot; will only issue fast reconfigurations

        """
        self.magnetic.set_color('blue')
        self.electric.set_color('green')
        self.helicity.set_color('red')
        self.monopole.set_color('purple')
        self.magnetic.set_visible(configuration.get('magnetic', True))
        self.electric.set_visible(configuration.get('electric', True))
        self.helicity.set_visible(configuration.get('helicity', False))
        self.monopole.set_visible(configuration.get('monopole', False))
        self.ax1.set_xscale('log' if configuration.get('logx', False) else 'linear')
        self.ax1.set_yscale('log' if configuration.get('logy', False) else 'linear')
        self.ax1.relim(visible_only=True)
        self.ax1.autoscale_view()


    def create_configuration(self, user_config):
        user_config.add_boolean("magnetic", True)
        user_config.add_boolean("electric", True)
        user_config.add_boolean("helicity", False)
        user_config.add_boolean("monopole", False)
        user_config.add_boolean("log x", False, alias='logx')
        user_config.add_boolean("log y", False, alias='logy')



class SingleImagePlot(object):
    def __init__(self, figure):
        from matplotlib import image
        from matplotlib import axes

        self.ax1 = axes.Axes(figure, [0.0, 0.0, 1.0, 1.0])
        self.image = image.AxesImage(self.ax1)
        self.image.set_data([[],[]])
        self.ax1.add_image(self.image)

        self.chkpt = None
        self.field = 'B3'

        self.configure()


    def get_axes(self):
        return [self.ax1]


    def reload_data(self, chkpt, field):
        import h5py
        import numpy as np
        if chkpt == self.chkpt and field == self.field: return
        if chkpt is None or field is None:
            data = np.array([[0,0], [0,0]])
            Li = 1
            Lj = 1
        else:
            which = dict(E='electric', B='magnetic', J='electric_current')[field[0]]
            h5f = h5py.File(chkpt, 'r')
            if len(h5f[which][field].shape) == 2:
                data = h5f[which][field][:,:].T
            else:
                data = h5f[which][field][:,:,0].T

            try:
                Li = h5f['sim']['domain_size[1]'][0]
                Lj = h5f['sim']['domain_size[2]'][0]
            except KeyError:
                Li = data.shape[1]
                Lj = data.shape[0]

            h5f.close()
        self.image.set_data(data)
        self.image.set_extent([0, Li, 0, Lj])
        if not self.fix_clim:
            self.image.set_clim(data.min(), data.max())
        self.ax1.set_aspect('equal')
        self.chkpt = chkpt
        self.field = field


    def set_data_source(self, data_source):
        """Heavier re-plot; may re-load data from the disk

        """
        rundir = data_source.get_rundir()
        chkpt = data_source.get_checkpoint()
        self.reload_data(chkpt, self.field)


    def configure(self, **configuration):
        self.fix_clim = configuration.get('fix_clim', False)
        self.image.set_cmap(configuration.get('cmap', 'jet'))
        self.ax1.xaxis.set_visible(configuration.get('draw_axes', False))
        self.ax1.yaxis.set_visible(configuration.get('draw_axes', False))
        field = configuration.get('field', 'B3')
        self.reload_data(self.chkpt, field)


    def create_configuration(self, user_config):
        import matplotlib
        cmaps = matplotlib.cm.datad.keys()
        user_config.add_boolean('Fix limits', default=False, alias='fix_clim')
        user_config.add_choices("Color map", cmaps, default='jet', alias='cmap')
        user_config.add_choices("Field",
                                ['E1', 'E2', 'E3',
                                 'B1', 'B2', 'B3',
                                 'J1', 'J2', 'J3'], default='B3', alias='field', radio=True)



class ProfilePlot(object):
    def __init__(self, figure):
        from matplotlib import lines
        from matplotlib import axes

        self.ax1 = axes.Axes(figure, [0.1, 0.1, 0.8, 0.8])
        self.B1 = lines.Line2D([], [])
        self.ax1.add_line(self.B1)

        self.chkpt = None
        self.field = 'B3'

        self.configure()


    def get_axes(self):
        return [self.ax1]


    def reload_data(self, chkpt, field):
        import h5py
        import numpy as np

        if chkpt is None or field is None:
            return

        which = dict(E='electric', B='magnetic', J='electric_current')[field[0]]

        h5f = h5py.File(chkpt, 'r')
        shape = h5f[which][field].shape
        if len(shape) == 2:
            x = np.linspace(0, 1, shape[0])
            y = h5f[which][field][:,shape[1]/2]
        h5f.close()

        self.B1.set_data(x, y)
        self.ax1.relim(visible_only=True)
        self.ax1.autoscale_view()
        self.chkpt = chkpt
        self.field = field


    def set_data_source(self, data_source):
        """Heavier re-plot; may re-load data from the disk

        """
        rundir = data_source.get_rundir()
        chkpt = data_source.get_checkpoint()
        self.reload_data(chkpt, self.field)


    def configure(self, **configuration):
        self.field = configuration.get('field', 'B3')
        self.reload_data(self.chkpt, self.field)


    def create_configuration(self, user_config):
        user_config.add_choices("Field",
                                ['E1', 'E2', 'E3',
                                 'B1', 'B2', 'B3',
                                 'J1', 'J2', 'J3'], default='B3', alias='field', radio=True)


class PowerSpectrumPlot(object):
    def __init__(self, figure):
        from matplotlib import lines
        from matplotlib import axes

        self.ax1 = axes.Axes(figure, [0.1, 0.1, 0.8, 0.8])
        self.electric = lines.Line2D([], [], lw=2, color='green')
        self.magnetic = lines.Line2D([], [], lw=2, color='blue')
        self.helicity = lines.Line2D([], [], lw=2, color='red')

        self.ax1.add_line(self.electric)
        self.ax1.add_line(self.magnetic)
        self.ax1.add_line(self.helicity)

        self.ax1.set_xscale('log')
        self.ax1.set_yscale('log')

        self.configure()


    def get_axes(self):
        return [self.ax1]


    def set_data_source(self, data_source):
        import h5py
        rundir = data_source.get_rundir()
        group = data_source.get_analysis_group()
        if rundir is None or group is None: return

        analysis_h5 = os.path.join(rundir, 'analysis.h5')
        h5f = h5py.File(analysis_h5, 'r')

        for which in ['electric', 'magnetic', 'helicity']:
            k = h5f[group][which]['binlocX'][:]
            P = h5f[group][which]['binval'][:]
            getattr(self, which).set_data(k, P)

        self.ax1.relim()
        self.ax1.autoscale_view()

        h5f.close()


    def configure(self, **configuration):
        pass


    def create_configuration(self, user_config):
        pass



class AlphaHistogramPlot(object):
    def __init__(self, figure):
        from matplotlib import lines
        from matplotlib import axes

        self.alp = lines.Line2D([], [], lw=1, color='k')
        self.ax1 = axes.Axes(figure, [0.1, 0.1, 0.8, 0.8])
        self.ax1.add_line(self.alp)
        self.ax1.set_xscale('linear')
        self.ax1.set_yscale('log')

        self.data_source = None

        self.downsample_factor = 0
        self.configure()


    def get_axes(self):
        return [self.ax1]


    def set_data_source(self, data_source=None):
        import h5py

        if data_source is None: data_source = self.data_source
        if data_source is None: return
        self.data_source = data_source

        if self.checkpoint:
            chkpt = data_source.get_checkpoint()

            if chkpt is None:
                return
            else:
                group = 'spectra'
                h5f = h5py.File(chkpt, 'r')

        else:
            group = data_source.get_analysis_group()

            if group is None:
                return
            else:
                rundir = data_source.get_rundir()
                analysis_h5 = os.path.join(rundir, 'analysis.h5')
                h5f = h5py.File(analysis_h5, 'r')

        try:
            d = lambda A: self.downsample(A, self.downsample_factor)
            alp_e = h5f[group]['alpha-hist']['binedgX'][:]
            alp_x = h5f[group]['alpha-hist']['binlocX'][:]
            alp_y = h5f[group]['alpha-hist']['binval'][:]

            dxa = alp_e[1:] - alp_e[:-1]
            tota = (dxa * alp_y).sum() # not correct, cow binmode set to 'counts'
            self.alp.set_data(d(alp_x), d(alp_y) / tota)

        except KeyError: # There was no data
            print "There was no alpha data"
            self.alp.set_data([], [])

        h5f.close()

        self.ax1.relim()
        self.ax1.autoscale_view()


    def configure(self, **configuration):
        self.checkpoint = configuration.get('checkpoint', False)
        self.downsample_factor = int(configuration.get('downsample', 0))
        self.set_data_source()


    def create_configuration(self, user_config):
        user_config.add_boolean('checkpoint', False)
        user_config.add_choices('downsample', ['0', '1', '2', '3', '4'])


    def downsample(self, A, factor):
        """
        Combine 2^factor adjacent freqency bins
        """
        if factor == 0:
            return A
        else:
            return self.downsample(0.5 * (A[0:-2:2] + A[1:-1:2]), factor-1)
