import numpy as np

def get_panel_centers_cylindrical(meshfile_path):
    """
    Lit un fichier de maillage au format 'meshfile' et renvoie un tableau des
    coordonnées cylindriques (r, theta, z) des centres des panneaux.

    Paramètres
    ----------
    meshfile_path : str
        Chemin vers le fichier meshfile.

    Retour
    ------
    cylindrical_centers : np.ndarray
        Tableau de forme (N, 3) contenant (r, theta, z) pour chaque panneau.
    """
    with open(meshfile_path, "r") as f:
        lines = f.readlines()

    point_dict = {}
    panels = []

    reading_nodes = True
    reading_panels = False

    for line in lines:
        if not line.strip():
            continue
        parts = line.strip().split()

        if reading_nodes:
            if parts[0] == '0':
                # Fin de la liste des noeuds
                reading_nodes = False
                reading_panels = True
                continue
            # Lecture noeud : ID x y z
            ID = int(parts[0])
            coords = list(map(float, parts[1:4]))
            point_dict[ID] = coords
        elif reading_panels:
            if parts[0] == '0':
                # Fin des panneaux
                break
            # Lecture panneau : indices des noeuds (sans ID panneau)
            ids = list(map(int, parts))
            panels.append(ids)

    cylindrical_centers = []
    areas = []

    for panel in panels:
        try:
            verts = np.array([point_dict[pid] for pid in panel])
        except KeyError as e:
            print(f"⚠️ ID de nœud non trouvé : {e}. Panneau ignoré.")
            continue

        # Centre
        center = verts.mean(axis=0)
        x, y, z = center
        r = np.sqrt(x**2 + y**2)
        theta = np.arctan2(y, x)
        cylindrical_centers.append([r, theta, z])

        # Aire du panneau rectangle
        if len(verts) == 4:
            AB = verts[1] - verts[0]
            AD = verts[3] - verts[0]
            area = np.linalg.norm(np.cross(AB, AD))
        else:
            area = 0.0  # ou gérer différemment si triangles ou autres
        areas.append(area)

    return np.array(cylindrical_centers), np.array(areas)


# meshfile = "../wec_inputs/PIT3_test/Cylinder.dat"
# centers_cyl, areas = get_panel_centers_cylindrical(meshfile)
# print("Aire totale :", np.sum(areas))
# print("Nombre de centres :", len(centers_cyl))
# # print("Centres cylindriques :")
# # print(centers_cyl)
# theta = centers_cyl[:,1]
# print("Theta min (rad) :", theta.min())
# print("Theta max (rad) :", theta.max())
# print("Theta min (deg) :", np.degrees(theta.min()))
# print("Theta max (deg) :", np.degrees(theta.max()))