#!/Users/lacquema/Oracle.env/bin/python3

import sys
import os
from PyQt6.QtWidgets import (
    QMainWindow, QVBoxLayout, QLabel, QStatusBar, QWidget, QApplication, QPushButton
)
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt


class WindowMenuClass(QMainWindow):
    """Fenêtre principale du menu, permettant de naviguer entre les simulations."""

    def __init__(self):
        super().__init__()

        # Définition du répertoire du fichier
        self.DirPath = os.path.dirname(__file__)

        # Configuration de la fenêtre
        self.setWindowTitle("Oracle 1.0")
        self.setFixedSize(850, 450)

        # Ajout d'une image de fond
        self.setStyleSheet(f"background-image: url({self.DirPath}/Items/LoadingBackground.png)")

        # Initialisation du layout principal
        self.init_ui()

    def init_ui(self):
        """Initialisation des widgets et du layout."""
        layout = QVBoxLayout()
        layout.addSpacing(20)

        # Création des boutons
        self.BtnNew = self.create_button("New Simulation")
        layout.addWidget(self.BtnNew, alignment=Qt.AlignmentFlag.AlignHCenter)

        layout.addSpacing(20)

        self.BtnContinue = self.create_button("Continue")
        layout.addWidget(self.BtnContinue, alignment=Qt.AlignmentFlag.AlignHCenter)

        layout.addSpacing(95)

        self.BtnAnalyse = self.create_button("Analyse")
        layout.addWidget(self.BtnAnalyse, alignment=Qt.AlignmentFlag.AlignHCenter)

        layout.addSpacing(40)  # Ajustement de l'espacement

        # Ajout de la barre de statut avec crédits
        self.init_status_bar()

        # Conteneur central
        container = QWidget()
        container.setLayout(layout)
        self.setCentralWidget(container)

    def create_button(self, text):
        """Crée un bouton stylisé avec une taille définie."""
        button = QPushButton(text)
        button.setFixedSize(200, 40)
        button.setStyleSheet(
            "QPushButton {background-color: grey; color: white; font: italic 15px;}"
        )
        return button

    def init_status_bar(self):
        """Initialise la barre de statut avec les crédits."""
        status_bar = QStatusBar(self)
        status_bar.addWidget(QLabel(" Version 1.0, 2024, IPAG, Hervé Beust, Antoine Lacquement"))
        self.setStatusBar(status_bar)


if __name__ == "__main__":
    app = QApplication(sys.argv)  # Création de l'application PyQt
    main_menu = WindowMenuClass()  # Affichage du menu principal
    sys.exit(app.exec())  # Exécution de l'application
