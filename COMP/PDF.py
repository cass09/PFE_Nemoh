import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from PIL import Image
import os

nom_pdf = "ConvergenceMesh_barge_Nb1.pdf"
dossier_graphes = "./Mesh"
fichiers_graphes = sorted([
    f for f in os.listdir(dossier_graphes)
    if f.startswith("N_") and f.endswith(".png")
])

with PdfPages(nom_pdf) as pdf:
    for fichier in fichiers_graphes:
        chemin = os.path.join(dossier_graphes, fichier)
        image = Image.open(chemin)
        fig, ax = plt.subplots(figsize=(8, 6))
        ax.imshow(image)
        ax.axis('off')
        ax.set_title(fichier, fontsize=10)
        pdf.savefig(fig, bbox_inches='tight')
        plt.close(fig)

# nom_pdf = "Graphs_errors.pdf"
# dossier_graphes = "."
# fichiers_graphes = sorted([
#     f for f in os.listdir(dossier_graphes)
#     if f.startswith("Errors_") and f.endswith(".png")
# ])

# with PdfPages(nom_pdf) as pdf:
#     for fichier in fichiers_graphes:
#         chemin = os.path.join(dossier_graphes, fichier)
#         image = Image.open(chemin)
#         fig, ax = plt.subplots(figsize=(8, 6))
#         ax.imshow(image)
#         ax.axis('off')
#         ax.set_title(fichier, fontsize=10)
#         pdf.savefig(fig, bbox_inches='tight')
#         plt.close(fig)
