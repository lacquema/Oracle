#!/Users/lacquema/Oracle.env/bin/python3

import sys
import os
from PyQt6.QtWidgets import QApplication, QMainWindow

# Ajout des chemins des modules internes
sys.path.append(os.path.dirname(__file__))
sys.path.append(os.path.join(os.path.dirname(__file__), 'AnaSimu'))
sys.path.append(os.path.join(os.path.dirname(__file__), 'NewSimu'))
sys.path.append(os.path.join(os.path.dirname(__file__), 'ContSimu'))

# Importation des fenêtres de l'application
from WindowMenu import WindowMenuClass
from NewSimu.WindowSetNewSimu import WindowSetNewSimu
from ContSimu.WindowSetContSimu import WindowSetContSimu
from AnaSimu.WindowSetAnaSimu import WindowSetAnaSimu


class MainWindow(QMainWindow):
    """Fenêtre principale gérant le menu et la navigation entre les fenêtres de simulation."""

    def __init__(self):
        super().__init__()

        # Initialisation du menu principal
        self.menu_window = WindowMenuClass()
        self.menu_window.show()

        # Initialisation des fenêtres secondaires
        self.new_simu_window = WindowSetNewSimu()
        self.cont_simu_window = WindowSetContSimu()
        self.ana_simu_window = WindowSetAnaSimu()

        # Connexion des signaux de fermeture des fenêtres secondaires
        self.new_simu_window.SignalCloseWindowSetNewSimu.connect(self.reopen_menu)
        self.cont_simu_window.SignalCloseWindowSetContSimu.connect(self.reopen_menu)
        self.ana_simu_window.SignalCloseWindowSetAnaSimu.connect(self.reopen_menu)
        self.ana_simu_window.ReSignalCloseWindowMain.connect(self.reopen_menu)

        # Connexion des boutons du menu aux méthodes d'ouverture des fenêtres
        self.menu_window.BtnNew.clicked.connect(lambda: self.open_window(self.new_simu_window))
        self.menu_window.BtnContinue.clicked.connect(lambda: self.open_window(self.cont_simu_window))
        self.menu_window.BtnAnalyse.clicked.connect(lambda: self.open_window(self.ana_simu_window))

    def open_window(self, window):
        """Ferme le menu et ouvre la fenêtre spécifiée."""
        self.menu_window.close()
        window.show()

    def reopen_menu(self):
        """Ferme toutes les fenêtres et réaffiche le menu principal."""
        QApplication.closeAllWindows()
        self.menu_window.show()


if __name__ == "__main__":
    app = QApplication(sys.argv)  # Création de l'application PyQt
    main_window = MainWindow()  # Instanciation et affichage de la fenêtre principale
    sys.exit(app.exec())  # Lancement de l'application
