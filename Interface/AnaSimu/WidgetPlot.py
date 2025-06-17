#! /Users/lacquema/Oracle.env/bin/python3


### --- Packages --- ###

# Transverse packages
import sys

# PyQt packages
from PyQt6.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QApplication

# Matplotlib packages
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


### --- Canvas Generating --- ###

class MplCanvas(FigureCanvasQTAgg):

    def __init__(self):
        self.fig = Figure()
        # self.fig.subplots_adjust(left=0.2, right=0.8, top=0.9, bottom=0.15)

        super(MplCanvas, self).__init__(self.fig)


### --- Plot Window Generating --- ###

class WidgetPlot(QWidget):

    def __init__(self, plot_function):
        super().__init__()

        # Layout
        self.Layout = QVBoxLayout()

        # Plot
        self.plot_function = plot_function

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
        self.history = []
        self.history_index = 0

        # Override the "home", "back" and "forward" buttons in the toolbar
        self.Toolbar._actions['home'].triggered.connect(self.restore_plot_to_initial_state)
        self.Toolbar._actions['back'].triggered.connect(self.undo_plot_state)
        self.Toolbar._actions['forward'].triggered.connect(self.redo_plot_state)

        # Flags to track events
        self.draw_event_triggered = False
        self.button_release_event_triggered = False

        # Connect events
        self.Canvas.mpl_connect('draw_event', self.on_draw_event)
        self.Canvas.mpl_connect('button_release_event', self.on_button_release_event)

        

    def plotting(self):
        self.clear_figure()
        self.plot_function()
        self.Canvas.draw()
        # if len(self.history)==0: self.save_plot_state()

    def refresh_plot(self):
        if self.isVisible():
            self.close_figure_options()
            self.plotting()
            self.restore_plot_state(index=self.history_index)
            self.update_toolbar_buttons()

    def reset_history(self):
        """Reset the history and clear the current state."""
        if self.isVisible():
            print('reset_event', self.history_index, len(self.history))
            # self.history = [self.history[0]]
            self.history = []
            self.history_index = 0
            self.update_toolbar_buttons()
            print('reset_accomplished', self.history_index, len(self.history))
        # print("History reset")

    def clear_figure(self):
        # for i in range(len(self.Canvas.fig.axes)):
        #     self.Canvas.fig.delaxes(self.Canvas.fig.axes[0])
        # for ax in self.Canvas.fig.axes:
        #     ax.cla()
        #     print(ax.get_xlim())
        self.Canvas.fig.clear()

    def close_figure_options(self):
        for widget in QApplication.instance().allWidgets():
                if widget.isWindow() and widget.isVisible():
                    if widget.windowTitle() == 'Figure options':
                        widget.close()

    def save_plot_labels(self):
        self.labels = {
        'title': self.Canvas.fig.axes[0].get_title(),
        'legend': [legend.get_text() for legend in self.Canvas.fig.axes[0].get_legend().get_texts()] if self.Canvas.fig.axes[0].get_legend() else None,
        'xlabel': self.Canvas.fig.axes[0].get_xlabel(),
        'ylabel': self.Canvas.fig.axes[0].get_ylabel(),
        }

    def on_draw_event(self, event):
        """Handle draw_event."""
        print('draw_event', self.history_index, len(self.history))
        self.save_plot_labels()
        self.draw_event_triggered = True
        self.check_and_save_state()
        # print(self.history)
        # self.save_plot_state()

    def on_button_release_event(self, event):
        """Handle button_release_event."""
        self.button_release_event_triggered = True

    def check_and_save_state(self):
        """Check if both events have been triggered and save the state."""
        # if self.draw_event_triggered and self.button_release_event_triggered:
        if self.Toolbar._actions['pan'].isChecked():
            condition = self.draw_event_triggered and self.button_release_event_triggered
            # print('pan')
        else:
            condition = self.draw_event_triggered
            # print('other')
        if condition:
            self.save_plot_state()
            # Reset flags after saving
            self.draw_event_triggered = False
            self.button_release_event_triggered = False

    def save_plot_state(self, event=None):
        """Save the current state of the toolbar and plot to history."""
        state = {
            'xlim': self.Canvas.fig.axes[0].get_xlim(),
            'ylim': self.Canvas.fig.axes[0].get_ylim(),
        }
        # Save zlim if the plot is 3D
        if hasattr(self.Canvas.fig.axes[0], 'get_zlim'):
            state['zlim'] = self.Canvas.fig.axes[0].get_zlim()
        if len(self.history) == 0: 
            self.remove_undone_states()
            self.history.append(state)
            self.history_index = len(self.history) - 1
            # print('save ', self.history_index, len(self.history))
            # print(state)
            if self.Toolbar._actions['pan'].isChecked(): self.Toolbar._actions['pan'].trigger()
            if self.Toolbar._actions['zoom'].isChecked(): self.Toolbar._actions['zoom'].trigger()
            self.update_toolbar_buttons()
        else: 
            if (state['xlim'] != self.history[self.history_index]['xlim']) or \
               (state['ylim'] != self.history[self.history_index]['ylim']) or \
               ('zlim' in state and state['zlim'] != self.history[self.history_index].get('zlim', None)):
               if (state['xlim'] != self.history[0]['xlim']) or \
                  (state['ylim'] != self.history[0]['ylim']) or \
                  ('zlim' in state and state['zlim'] != self.history[0].get('zlim', None)):
                    self.remove_undone_states()
                    self.history.append(state)
                    self.history_index = len(self.history) - 1
                    # print('save ', self.history_index, len(self.history))
                    # print(state)
                    if self.Toolbar._actions['pan'].isChecked(): self.Toolbar._actions['pan'].trigger()
                    if self.Toolbar._actions['zoom'].isChecked(): self.Toolbar._actions['zoom'].trigger() 
                    self.update_toolbar_buttons()       

    def remove_undone_states(self):
        """Remove any states that were undone."""
        if self.history_index != len(self.history)-1:
            self.history = self.history[:self.history_index+1]

    def restore_plot_state(self, index, limites=True, labels=True):
        """Restore the current state of the toolbar and plot."""
        state = self.history[index]
        ax = self.Canvas.fig.axes[0]
        if limites==True:
            ax.set_xlim(state['xlim'])
            ax.set_ylim(state['ylim'])
            # Restore zlim if the plot is 3D
            if 'zlim' in state and hasattr(ax, 'set_zlim'):
                ax.set_zlim(state['zlim'])
        if labels==True:
            ax.set_xlabel(self.labels['xlabel'])
            ax.set_ylabel(self.labels['ylabel'])
            ax.set_title(self.labels['title'])
            if self.labels['legend']: ax.legend(labels=self.labels['legend'])
        self.Canvas.draw()
        # self.adapt_tight()

    def undo_plot_state(self):
        """Revert to the previous state in history."""
        if self.history_index > 0: 
            self.history_index -= 1
            # print('undo ', self.history_index, len(self.history))
            self.restore_plot_state(index=self.history_index)
            self.update_toolbar_buttons()
            # self.adapt_tight()

    def redo_plot_state(self):
        """Advance to the next state in history."""
        if self.history_index < len(self.history) - 1:
            self.history_index += 1
            # print('redo ', self.history_index, len(self.history))
            self.restore_plot_state(index=self.history_index)
            self.update_toolbar_buttons()
            # self.adapt_tight()

    def restore_plot_to_initial_state(self):
        """Restor the toolbar to its initial state."""
        self.remove_undone_states()
        if self.history[0]!= self.history[-1]:
            self.history.append(self.history[0])
        self.history_index = len(self.history)-1
        # print('home ', self.history_index, len(self.history))
        self.restore_plot_state(index=0)
        self.update_toolbar_buttons()
        # self.adapt_tight()

    def update_toolbar_buttons(self):
        """Update the state of the undo and redo buttons."""
        # print('back ',self.history_index > 0)
        self.Toolbar._actions['back'].setEnabled(self.history_index > 0)
        # print('forward ',self.history_index < len(self.history) - 1)
        self.Toolbar._actions['forward'].setEnabled(self.history_index < len(self.history) - 1)

    def adapt_tight(self):
        for i in range(15):
            self.Canvas.fig.tight_layout()
            


if __name__=="__main__":
    app = QApplication(sys.argv)
    Plot = WidgetPlot()
    Window = QMainWindow()
    Window.setCentralWidget(Plot)
    Window.show()
    app.exec()