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
    First=True
    for line in lines:
        if not line.strip():
            continue
        parts = line.strip().split()
        if First:
            sym=int(parts[1])
            First=False
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

        def compute_one_panel(verts_panel):
            center = verts_panel.mean(axis=0)
            x, y, z = center
            r = np.sqrt(x**2 + y**2)
            theta = np.arctan2(y, x)
            if len(verts_panel) == 4 and not np.allclose(verts_panel[0], verts_panel[1]):
                AB = verts_panel[1] - verts_panel[0]
                AD = verts_panel[3] - verts_panel[0]
                area = np.linalg.norm(np.cross(AB, AD))
            else:
                AB = verts_panel[1] - verts_panel[0]
                AC = verts_panel[2] - verts_panel[0]
                area = 0.5 * np.linalg.norm(np.cross(AB, AC))
            return [r, theta, z], area

        # Panneau original
        r_theta_z, area = compute_one_panel(verts)
        cylindrical_centers.append(r_theta_z)
        areas.append(area)
        # Si symétrie, ajouter le panneau miroir
        if sym == 1:
            verts_mirror = verts.copy()
            verts_mirror[:, 1] *= -1  # y -> -y
            verts_mirror = verts_mirror[::-1]  # inverser l’ordre des points pour conserver orientation normale
            r_theta_z_mirror, area_mirror = compute_one_panel(verts_mirror)
            cylindrical_centers.append(r_theta_z_mirror)
            areas.append(area_mirror)
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

def OCmesh(R, Ntheta, Nz, depth):
    # DO i=1,cyldZ
    #         DO j=1,cyldTheta
    #             IF (cyldZ .EQ. 1) THEN
    #                 cylZ = 0
    #             ELSE
    #                 cylZ = -Environment%Depth*(1.-COS(PI/2.*(i-1.)/(cyldZ-1.)))
    #             END IF
    #             WRITE(11,'(3(X,E14.6))') cylR*COS(2.*PI*(j-1)/cyldTheta),cylR*SIN(2.*PI*(j-1)/cyldTheta),cylZ
    #         END DO
    #     END DO  
    mesh=[]
    for i in range(Nz):
        for j in range(Ntheta):
            if Nz==1:
                z=0
            else : 
                z=-depth*(1-np.cos(np.pi/2*(i/(Nz-1))))
            coord=[R*np.cos(2*np.pi*j/Ntheta), R*np.sin(2*np.pi*j/Ntheta), z]
            mesh.append(coord)

    return mesh