#! /Users/lacquema/Oracle.env/bin/python3


### --- Packages --- ###

# Transverse packages
import sys

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QApplication
from PyQt6.QtCore import QEvent

# Matplotlib packages
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure

import time


### --- Canvas Generating --- ###

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self):
        self.fig = Figure()
        # self.fig.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.15)

        super(MplCanvas, self).__init__(self.fig)


### --- Plot Window Generating --- ###

class WidgetPlot(QWidget):

    def __init__(self, plot_function, xlim, ylim, zlim, azim, elev, xlabel, ylabel, zlabel, title, legend):
        """
        Initialize the WidgetPlot with plot parameters, layout, canvas, and toolbar.
        Stores plot state and connects toolbar actions to custom handlers.
        """
        super().__init__()

        # Layout
        self.Layout = QVBoxLayout()

        # Plot
        self.plot_function = plot_function
        self.xlim = xlim
        self.ylim = ylim
        self.zlim = zlim
        self.azim = azim
        self.elev = elev
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.zlabel = zlabel
        self.title = title
        self.legend = legend

        self.last_azim = -60
        self.last_elev = 30
        
        # Canvas initialization
        self.Canvas = MplCanvas()

        # Toolbar
        self.Toolbar = NavigationToolbar(self.Canvas, self)
        self.Layout.addWidget(self.Toolbar)
        
        # Plot on the Canvas
        self.Layout.addWidget(self.Canvas)

        # Widget container
        self.setLayout(self.Layout)

        # Store toolbar and plot state
        self.toolbar_state = {}
        self.plot_state = {}

        # Historique des états
        self.labels = {}
        self.history = []
        self.history_index = 0

        # Override the "home", "back" and "forward" buttons in the toolbar
        self.Toolbar._actions['home'].triggered.connect(self.restore_plot_to_initial_state)
        self.Toolbar._actions['back'].triggered.connect(self.undo_plot_state)
        self.Toolbar._actions['forward'].triggered.connect(self.redo_plot_state)

        self.Toolbar._actions['pan'].triggered.connect(self.auto_refresh_after_pan_or_zoom)
        self.Toolbar._actions['zoom'].triggered.connect(self.auto_refresh_after_pan_or_zoom)

        # Flags to track events
        self.button_release_event_triggered = False

        # Connect events
        self.draw_event = self.Canvas.mpl_connect('draw_event', self.on_draw_event)
        self.Canvas.mpl_connect('button_release_event', self.on_button_release_event)

    def auto_refresh_after_pan_or_zoom(self):
        """
        Automatically refresh the plot after a pan or zoom action is completed.
        If pan or zoom is still active, do nothing. Otherwise, refresh the plot.
        """
        if self.Toolbar._actions['pan'].isChecked() or self.Toolbar._actions['zoom'].isChecked():
            return
        else:
            self.refresh_plot()

    def plotting(self):
        """
        Clear the figure and call the user-provided plot function to draw the plot.
        """
        self.clear_figure()
        self.plot_function()
        self.Canvas.fig.subplots_adjust(left=0.15, right=0.93, top=0.9, bottom=0.12)

    def refresh_plot(self):
        """
        Refresh the plot.
        Closes any open figure options, redraws the plot, restores the current state, and updates toolbar buttons.
        """
        self.close_figure_options()
        self.plotting()
        self.restore_plot(index=self.history_index)
        self.Canvas.draw()
        if self.detect_plot_state_change():
            self.save_plot_state()
        if self.detect_plot_labels_change():
            self.save_plot_labels()
        self.update_toolbar_buttons()

    def reset_plot(self):
        """
        Reset the plot state history and clear the current state.
        """
        self.labels = {}
        self.history = []
        self.history_index = 0
        self.update_toolbar_buttons()

    def clear_figure(self):
        """
        Clear all axes from the current figure.
        """
        # for i in range(len(self.Canvas.fig.axes)):
        #     self.Canvas.fig.delaxes(self.Canvas.fig.axes[0])
        # for ax in self.Canvas.fig.axes:
        #     ax.cla()
        #     print(ax.get_xlim())
        self.Canvas.fig.clear()

    def on_draw_event(self, event):
        """
        Handle the Matplotlib draw_event.
        Save plot labels and check if the plot state should be saved.
        """
        if self.detect_plot_labels_change():
            self.save_plot_labels()
        self.update_toolbar_buttons()


    def on_button_release_event(self, event):
        """
        Handle the Matplotlib button_release_event.
        Detect if azim/elev have changed and call a callback if so.
        """        
        if self.detect_plot_state_change():
            self.save_plot_state()
            self.update_toolbar_buttons()

    def look_plot_labels(self):
        ax = self.Canvas.fig.axes[0]
        self.current_labels = {}
        if self.title:
            self.current_labels['title'] = ax.get_title()
        if self.xlabel:
            self.current_labels['xlabel'] = ax.get_xlabel()
        if self.ylabel:
            self.current_labels['ylabel'] = ax.get_ylabel()
        if self.zlabel and hasattr(ax, 'get_zlabel'):
            self.current_labels['zlabel'] = ax.get_zlabel()
        if self.legend and ax.get_legend():
            self.current_labels['legend'] = [text.get_text() for text in ax.get_legend().get_texts()]
        return self.current_labels

    def detect_plot_labels_change(self):

        self.current_labels = self.look_plot_labels()

        if len(self.history) == 0 or self.dicos_differs(self.current_labels, self.labels):
            return True
        else:
            return False
        
    def save_plot_labels(self):

        self.labels = self.current_labels
        # print('save labels')
        self.update_toolbar_buttons()

    def dicos_differs(self, dico1, dico2):
        for key in dico1:
            if key not in dico2 or dico1[key] != dico2[key]:
                return True
        return False
            
    def look_plot_state(self):
        ax = self.Canvas.fig.axes[0]
        self.current_state = {}
        if self.xlim:
            self.current_state['xlim'] = ax.get_xlim()
        if self.ylim:
            self.current_state['ylim'] = ax.get_ylim()
        # Save zlim, azimute and elevation if the plot is 3D and requested
        if self.zlim and hasattr(ax, 'get_zlim'):
            self.current_state['zlim'] = ax.get_zlim()
        if self.azim and hasattr(ax, 'azim'):
            self.current_state['azim'] = ax.azim
        if self.elev and hasattr(ax, 'elev'):
            self.current_state['elev'] = ax.elev
        return self.current_state

    def detect_plot_state_change(self):
        """
        Check if the necessary events have been triggered and save the plot state if appropriate.
        Resets event flags after saving.
        """
        self.current_state = self.look_plot_state()
    
        if len(self.history) == 0 or (self.dicos_differs(self.current_state, self.history[self.history_index]) and self.dicos_differs(self.current_state, self.history[0])):
            return True
        else:
            return False
        
    def save_plot_state(self):
        self.remove_undone_states()
        self.history.append(self.current_state)
        self.history_index = len(self.history) - 1
        if self.Toolbar._actions['pan'].isChecked():
            self.Toolbar._actions['pan'].trigger()
        if self.Toolbar._actions['zoom'].isChecked():
            self.Toolbar._actions['zoom'].trigger()
        # print('save state', self.history_index, len(self.history))
        self.update_toolbar_buttons()

    def remove_undone_states(self):
        """
        Remove any states from the history that were undone (after an undo action).
        """
        if self.history_index != len(self.history)-1:
            self.history = self.history[:self.history_index+1]

    def restore_plot(self, index):
        """
        Restore the plot state (axes limits, labels, legend) from the history at the given index.
        """
        if not self.history or index >= len(self.history):
            return
        self.current_state = self.history[index]
        if not self.Canvas.fig.axes:
            return
        ax = self.Canvas.fig.axes[0]
        if 'xlim' in self.current_state:
            ax.set_xlim(self.current_state['xlim'])
        if 'ylim' in self.current_state:
            ax.set_ylim(self.current_state['ylim'])
        # Restore zlim, azimute and elevation if the plot is 3D
        if 'zlim' in self.current_state and hasattr(ax, 'set_zlim'):
            ax.set_zlim(self.current_state['zlim'])
        if 'azim' in self.current_state and hasattr(ax, 'azim'):
            ax.azim = self.current_state['azim']
        if 'elev' in self.current_state and hasattr(ax, 'elev'):
            ax.elev = self.current_state['elev']
        # Restore labels and legend
        if 'xlabel' in self.labels:
            ax.set_xlabel(self.labels['xlabel'])
        if 'ylabel' in self.labels:
            ax.set_ylabel(self.labels['ylabel'])
        if 'zlabel' in self.labels and hasattr(ax, 'set_zlabel'):
            ax.set_zlabel(self.labels['zlabel'])
        if 'title' in self.labels:
            ax.set_title(self.labels['title'])
        if 'legend' in self.labels and self.labels['legend']:
            handles, _ = ax.get_legend_handles_labels()
            if handles:
                ax.legend(handles, self.labels['legend'])

    def undo_plot_state(self):
        """
        Revert to the previous plot state in the history (undo action).
        """
        if self.history_index > 0: 
            self.history_index -= 1
            # print('undo ', self.history_index, len(self.history))
            self.plotting()
            self.restore_plot(index=self.history_index)
            self.Canvas.draw()
            self.update_toolbar_buttons()

    def redo_plot_state(self):
        """
        Advance to the next plot state in the history (redo action).
        """
        if self.history_index < len(self.history) - 1:
            self.history_index += 1
            # print('redo ', self.history_index, len(self.history))
            self.plotting()
            self.restore_plot(index=self.history_index)
            self.Canvas.draw()
            self.update_toolbar_buttons()

    def restore_plot_to_initial_state(self):
        """
        Restore the plot and toolbar to the initial state (home action).
        """
        self.remove_undone_states()
        if self.history[0]!= self.history[-1]:
            self.history.append(self.history[0])
        self.history_index = len(self.history)-1
        # print('home ', self.history_index, len(self.history))
        self.plotting()
        self.restore_plot(index=0)
        self.Canvas.draw()
        self.update_toolbar_buttons()

    def update_toolbar_buttons(self):
        """
        Update the enabled/disabled state of the undo and redo toolbar buttons.
        """
        # print('back ',self.history_index > 0)
        self.Toolbar._actions['back'].setEnabled(self.history_index > 0)
        # print('forward ',self.history_index < len(self.history) - 1)
        self.Toolbar._actions['forward'].setEnabled(self.history_index < len(self.history) - 1)

    def adapt_tight(self):
        """
        Apply tight layout to the figure multiple times to optimize spacing.
        """
        for i in range(10):
            self.Canvas.fig.tight_layout()

    def define_closeEvent_figure_options(self, event):
        for widget in QApplication.instance().allWidgets():
            if widget.isWindow() and widget.windowTitle() == 'Figure options':
                widget.destroyed.connect(lambda: print('figure options closed'))

    def close_figure_options(self):
        """
        Close any open 'Figure options' windows in the application.
        """
        for widget in QApplication.instance().allWidgets():
            if widget.isWindow() and widget.isVisible():
                if widget.windowTitle() == 'Figure options':
                    widget.close()

            


if __name__=="__main__":
    app = QApplication(sys.argv)
    Plot = WidgetPlot()
    Window = QMainWindow()
    Window.setCentralWidget(Plot)
    Window.show()
    app.exec()