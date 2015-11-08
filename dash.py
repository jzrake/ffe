#!/usr/bin/env python


import os
import argparse

try: import wx
except ImportError:
    print "Problem: wxPython is a dependency"



class HasListeners(object):
    def __init__(self):
        self.listeners = [ ]

    def add_listener(self, listener):
        self.listeners.append(listener)

    def notify_listeners(self):
        for listener in self.listeners:
            listener.receive_change(self)



class UserConfiguration(wx.Panel, HasListeners):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        HasListeners.__init__(self)

        self.widgets = { }
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.SetSizer(self.sizer)


    def add_boolean(self, label, default, alias=None):
        if alias is None: alias = label
        cb = wx.CheckBox(self, label=label)
        cb.SetValue(default)
        cb.Bind(wx.EVT_CHECKBOX, self.value_changed)
        self.sizer.Add(cb)
        self.sizer.Layout()
        self.widgets[alias] = cb


    def add_choices(self, label, choices, default=None, alias=None, radio=False):
        if alias is None: alias = label

        if radio:
            cb = wx.RadioBox(self, label=label, choices=choices, style=wx.RA_SPECIFY_ROWS)
            cb.SetSelection(choices.index(default))
            cb.Bind(wx.EVT_RADIOBOX, self.value_changed)
        else:
            cb = wx.ComboBox(self, style=wx.CB_READONLY, name=label, choices=choices)
            cb.Bind(wx.EVT_COMBOBOX, self.value_changed)
            cb.SetValue(default if default is not None else choices[0])

        self.sizer.Add(cb)
        self.sizer.Layout()
        self.widgets[alias] = cb


    def value_changed(self, event):
        self.notify_listeners()


    def get_configuration(self):
        cfg = dict()
        for alias in self.widgets:
            try:
                cfg[alias] = self.widgets[alias].GetValue()
            except AttributeError:
                n = self.widgets[alias].GetSelection()
                cfg[alias] = self.widgets[alias].GetString(n)
        return cfg



class DataSourcePanel(wx.Panel, HasListeners):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        HasListeners.__init__(self)

        self.listbox_run = wx.ListBox(self)
        self.listbox_chk = wx.ListBox(self)
        self.listbox_ana = wx.ListBox(self)

        self.listbox_run.Bind(wx.EVT_LISTBOX, self.evt_listbox_run)
        self.listbox_chk.Bind(wx.EVT_LISTBOX, self.evt_listbox_chk)
        self.listbox_ana.Bind(wx.EVT_LISTBOX, self.evt_listbox_ana)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.listbox_run, 1, wx.GROW)
        sizer.Add(self.listbox_chk, 1, wx.GROW)
        sizer.Add(self.listbox_ana, 1, wx.GROW)
        self.SetSizer(sizer)


    def load_run(self, rundir, load_chk=True, load_ana=True):
        runid = os.path.basename(rundir)

        has_ffe_dat = os.path.isfile(os.path.join(rundir, 'ffe.dat'))
        has_analysis_h5 = os.path.isfile(os.path.join(rundir, 'analysis.h5'))

        if not (has_ffe_dat and has_analysis_h5):
            print "Not a valid run directory", rundir
        elif self.listbox_run.FindString(runid) != -1:
            self.make_run_active(rundir) # will reload checkpoints
        else:
            self.listbox_run.Append(runid, clientData=rundir)


    def make_run_active(self, rundir):
        import pickle
        import h5py


        analysis_h5 = os.path.join(rundir, "analysis.h5")
        manifest_pk = os.path.join(rundir, "manifest.pk")

        prev_chk = self.listbox_chk.GetSelection()
        prev_ana = self.listbox_ana.GetSelection()

        self.listbox_chk.Clear()
        self.listbox_ana.Clear()

        for thing in os.listdir(rundir):
            if thing.startswith('chkpt.') and thing.endswith('.h5'):
                self.listbox_chk.Append(thing)

        try:
            groups = pickle.load(open(manifest_pk, 'r'))
        except:
            print "caching analysis groups for", rundir
            groups = [ ]
            h5f = h5py.File(analysis_h5, "r")
            for group in h5f:
                print group
                groups.append(group)
            h5f.close()
            pickle.dump(groups, open(manifest_pk, 'w'))

        for group in groups:
            self.listbox_ana.Append(group)

        if prev_ana != -1 and prev_ana < self.listbox_ana.GetCount():
            self.listbox_ana.SetSelection(prev_ana)
        if prev_chk != -1 and prev_chk < self.listbox_chk.GetCount():
            self.listbox_chk.SetSelection(prev_chk)


    def get_rundir(self):
        n = self.listbox_run.GetSelection()
        if n < 0:
            return None
        else:
            return self.listbox_run.GetClientData(n)


    def get_checkpoint(self):
        rundir = self.get_rundir()
        if rundir is None:
            return None
        else:
            n = self.listbox_chk.GetSelection()
            if n < 0:
                return None
            else:
                return os.path.join(rundir, self.listbox_chk.GetString(n))
            

    def get_analysis_group(self):
        rundir = self.get_rundir()
        if rundir is None:
            return None
        else:
            n = self.listbox_ana.GetSelection()
            if n < 0:
                return None
            else:
                return self.listbox_ana.GetString(n)



    """Callbacks bound to things

    """
    def evt_listbox_run(self, event):
        rundir = self.get_rundir()
        if rundir is not None:
            self.make_run_active(rundir)
            self.notify_listeners()


    def evt_listbox_chk(self, event):
        self.notify_listeners()


    def evt_listbox_ana(self, event):
        self.notify_listeners()



class MainWindow(wx.Frame):

    def __init__(self, directories):
        wx.Frame.__init__(self, None, -1, 'FFE Dashboard')

        self.directories = directories
        self.import_matplotlib()
        self.create_panel()


    def import_matplotlib(self):
        import matplotlib

        matplotlib.use('WXAgg')

        from matplotlib.figure import Figure
        from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg
        from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg

        self.Figure = Figure
        self.FigureCanvas = FigureCanvasWxAgg
        self.NavigationToolbar = NavigationToolbar2WxAgg


    def create_panel(self):
        import matplotlib.pyplot as plt

        panel = wx.Panel(self)

        figure = self.Figure([6, 6])
        canvas = self.FigureCanvas(panel, -1, figure)
        navbar = self.NavigationToolbar(canvas)


        import plots
        self.plots = { }
        self.plots[      'TimeSeriesPlot'] = plots.      TimeSeriesPlot(figure)
        self.plots[     'SingleImagePlot'] = plots.     SingleImagePlot(figure)
        self.plots[         'ProfilePlot'] = plots.         ProfilePlot(figure)
        self.plots[   'PowerSpectrumPlot'] = plots.   PowerSpectrumPlot(figure)
        self.plots[  'AlphaHistogramPlot'] = plots.  AlphaHistogramPlot(figure)

        for plot in self.plots.values():
            user_configuration = UserConfiguration(panel)
            user_configuration.Hide()
            user_configuration.add_listener(self)
            plot.create_configuration(user_configuration)
            plot.user_configuration = user_configuration


        panel.data_source = DataSourcePanel(panel)
        panel.canvas = canvas
        panel.navbar = navbar
        panel.button_refresh = wx.Button(panel, label="Refresh directories")
        panel.button_recache = wx.Button(panel, label="Recache current run")
        panel.radio_plots = wx.RadioBox(panel, label="", choices=self.plots.keys(),
                                        style=wx.RA_SPECIFY_ROWS)



        panel.button_refresh.Bind(wx.EVT_BUTTON, self.refresh_directories)
        panel.button_recache.Bind(wx.EVT_BUTTON, self.recache_current_run)
        panel.radio_plots.Bind(wx.EVT_RADIOBOX, self.on_radio_plots)
        panel.data_source.add_listener(self)

        vsizer = wx.BoxSizer(wx.VERTICAL)
        hsizer = wx.BoxSizer(wx.HORIZONTAL)

        vsizer.Add(panel.navbar)
        vsizer.Add(panel.button_refresh)
        vsizer.Add(panel.button_recache)
        vsizer.Add(panel.radio_plots)

        for plot in self.plots.values():
            vsizer.Add(plot.user_configuration, 1, wx.EXPAND)

        hsizer.Add(panel.data_source, 1, wx.EXPAND)
        hsizer.Add(vsizer, 0, wx.EXPAND)
        hsizer.Add(panel.canvas, 4, wx.EXPAND)
        panel.SetSizer(hsizer)
        hsizer.Fit(self)


        self.figure = figure
        self.panel = panel
        self.plot = None
        self.on_radio_plots()
        self.refresh_directories()


    def switch_active_plot(self, plot_name):
        
        plot = self.plots[plot_name]
        plot.set_data_source(self.panel.data_source)

        if self.plot is not None: self.plot.user_configuration.Hide()
        plot.user_configuration.Show()

        self.panel.Layout()

        for ax in self.figure.axes: self.figure.delaxes(ax)
        for ax in plot.get_axes(): self.figure.add_axes(ax)
        self.panel.canvas.draw()
        self.plot = plot


    def receive_change(self, broadcaster):

        if self.plot is None: return

        if broadcaster is self.panel.data_source:
            self.plot.set_data_source(self.panel.data_source)
            self.panel.canvas.draw()
        elif broadcaster is self.plot.user_configuration:
            cfg = self.plot.user_configuration.get_configuration()
            self.plot.configure(**cfg)
            self.panel.canvas.draw()


    def on_radio_plots(self, event=None):
        n = self.panel.radio_plots.GetSelection()
        plot_name = self.panel.radio_plots.GetString(n)
        self.switch_active_plot(plot_name)


    def refresh_directories(self, event=None):
        for directory in self.directories:
            #for runid in os.listdir(directory):
            #self.panel.data_source.load_run(os.path.join(directory, runid))
            self.panel.data_source.load_run(directory)


    def recache_current_run(self, event=None):
        rundir = self.panel.data_source.get_rundir()
        try: os.remove(os.path.join(rundir, 'tseries.pk'))
        except OSError: pass
        try: os.remove(os.path.join(rundir, 'manifest.pk'))
        except OSError: pass



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="FFE run viewer --- needs: numpy, h5py, matplotlib, wxPython")
    parser.add_argument("rundir", nargs="*")

    args = parser.parse_args()

    app = wx.App(False)
    app.frame = MainWindow(args.rundir)
    app.frame.Show()
    app.MainLoop()
