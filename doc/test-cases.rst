
##########
Test cases
##########

The following test cases are provided for verification with the original Aquaplus software (which is the ancestor of Nemoh) and/or HydroSTAR commercial software :cite:p:`HYDROSTAR`.
The input data for each case is available in the ``TestCases`` directory of the repository.

.. note:: Tecplot layout files (``.lay``) are provided in the relevant test case folders for easier plotting of the comparisons in Tecplot.

**********
1_Cylinder
**********

This example is a floating vertical cylinder with the usual 6 degrees of freedom (surge, sway, heave, roll, pitch, yaw).
It features a half-symmetric body mesh, deep water, and a wave direction :math:`\beta = 0^{\circ}`.
The results are shown in :numref:`fig:Cylinder`.

.. _`fig:Cylinder`:
.. figure:: figures/Ver_Cylinder.svg
   :align: center

   Comparison of the first order results between Nemoh and AQUAPLUS for the test case **1_Cylinder**

*********
2_2Bodies
*********

This example deals with the hydrodynamic interactions between two bodies.
The bodies are both cubes, and on of them is rotated by 45° degrees around the z axis.
Each body has the usual 6 degrees of freedom (surge, sway, heave, roll, pitch, yaw).
It features half-symmetric body meshes, two different bodies, a water depth :math:`D=20` m and a wave direction :math:`\beta=45^{\circ}`. The results are shown in :numref:`fig:2Bodies`.

.. _`fig:2Bodies`:
.. figure:: figures/Ver_2Bodies.svg
   :align: center

   Comparison of the first order results between Nemoh and AQUAPLUS for the test case **2_2Bodies**

****************
3_Nonsymmetrical
****************

This example deals with a body for which the (xOz) plane is not a symmetry plane.
The body has the usual 6 degrees of freedom (surge, sway, heave, roll, pitch, yaw).
It features a full non-symmetrical body mesh, deep water and a wave direction :math:`\beta=0^{\circ}`.
Comparisons of Nemoh results against Aquaplus are shown in :numref:`fig:NonSymmetrical_1` and :numref:`fig:NonSymmetrical_2`, a slight difference is observed in the results.
Added mass and damping coefficients comparison between Nemoh and HydroSTAR are shown in :numref:`fig:NonSymmetrical_mass` and :numref:`fig:NonSymmetrical_damp`, and for the excitation force is in :numref:`fig:NonSymmetrical_excforce`.
Good agreement between Nemoh and HydroSTAR is achieved.

.. _`fig:NonSymmetrical_1`:
.. figure:: figures/Ver_NonSymmetrical_1.svg
   :align: center

   Comparison of the first order results between Nemoh and AQUAPLUS for the test case **3_Nonsymmetrical**

.. _`fig:NonSymmetrical_2`:
.. figure:: figures/Ver_NonSymmetrical_2.svg
   :align: center

   Comparison of the first order results between Nemoh and AQUAPLUS for the test case **3_Nonsymmetrical**

.. _`fig:NonSymmetrical_mass`:
.. figure:: figures/Ver_NonSymmetrical_addedmass.svg
   :align: center

   Comparison of added mass coefficients between Nemoh, red dashed-line, and HydroSTAR, blue solid-line, for the test case **3_Nonsymmetrical**

.. _`fig:NonSymmetrical_damp`:
.. figure:: figures/Ver_NonSymmetrical_dampcoef.svg
   :align: center

   Comparison of damping coefficients between Nemoh, red dashed-line, and HydroSTAR, blue solid-line, for the test case **3_Nonsymmetrical**

.. _`fig:NonSymmetrical_excforce`:
.. figure:: figures/Ver_NonSymmetrical_excitationforce.svg
   :align: center

   Comparison of excitation force between Nemoh, red dashed-line, and HydroSTAR, blue solid-line, for the test case **3_Nonsymmetrical**

****************
4_Postprocessing
****************

This example displays some of the postprocessing capabilities of Nemoh, namely:

- Computation of free surface elevation
- Computation of far field coefficients (Kochin functions)
- Computation of the IRF (Impulse Response Function)

.. note:: Nemoh's postprocessing capabilities also include RAO (Response Amplitude Operator) computation and pressure outputting, not featured in this test.

The body is a cube with the usual 6 degrees of freedom (surge, sway, heave, roll, pitch, yaw).
It featutres a half-symmetric body mesh, a water depth :math:`D=10` m and a wave direction :math:`\beta=0^{\circ}`.
This test case shows a comparison of the free-surface elevation and the Kochin function. The results are shown in :numref:`fig:PostProcessing`.
The phase difference :math:`\pm \pi/2` of wave elevation between Nemoh and AQUAPLUS is due to different conventions of the incident potential.

.. _`fig:PostProcessing`:
.. figure:: figures/Ver_PostProcessing.svg
   :align: center

   Comparison of the diffracted wave elevation, the diffraction Kochin function between Nemoh and AQUAPLUS, test case **4_Postprocessing**

***********
5_QuickTest
***********

Shows a quantitative comparison of force and free-surface for the first-frequency diffraction problem. The comparison results are shown in the command window for all the test cases inside the directory ``5_QuickTest``.

****************
6_box_coarsemesh
****************

Shows the procedure for running the code starting with the executable ``mesh`` with a coarse description mesh file, ``meshbox``. No reference data is given in this test case.

*******************
7_Solvers_Check_OC3
*******************

Tests the performance of the three difference linear solvers, Gauss elimination, LU decomposition and GMRES. Reference logfiles reporting the computational time of the solvers are provided.

*********************
8a_Cylinder_irregfreq
*********************

Shows the results with and without irregular frequencies removal (IRR) method. The results are verified against HydroSTAR with IRR and shown in :numref:`fig:Cylinder_IRR_addedmass` and :numref:`fig:Cylinder_IRR_dampcoef` for the added mass and damping coefficients and in :numref:`fig:Cylinder_IRR_excforce` for the excitation forces. The mesh used was obtained using GMSH :cite:p:`GMSH` and is shown in :numref:`fig:meshesCylinder`.

.. _`fig:meshesCylinder`:
.. figure:: figures/Cylinder/mesh.svg
   :align: center

   Body boundary mesh for the Cylinder used for test case **8a_Cylinder_irregfreq** and **8b_QTF_Cylinder**.

.. _`fig:Cylinder_IRR_addedmass`:
.. figure:: figures/Cylinder/addedmass.svg
   :align: center

   Comparison of added masscoefficients between Nemoh without irregular frequencies removal (IRR), green dash-dotted line, Nemoh with IRR, red dashed-line and HydroSTAR with IRR, blue solid-line, for the test-case **8a_Cylinder_irregfreq**

.. _`fig:Cylinder_IRR_dampcoef`:
.. figure:: figures/Cylinder/dampcoef.svg
   :align: center

   Comparison of damping coefficients between Nemoh without irregular frequencies removal (IRR), green dash-dotted line, Nemoh with IRR, red dashed-line and HydroSTAR with IRR, blue solid-line, for the test-case **8a_Cylinder_irregfreq**

.. _`fig:Cylinder_IRR_excforce`:
.. figure:: figures/Cylinder/excForce.svg
   :align: center

   Comparison of excitation force between Nemoh without irregular frequencies removal (IRR), green dash-dotted line, Nemoh with IRR, red dashed-line and HydroSTAR with IRR, blue solid-line, for the test-case 8a_Cylinder_irregfreq

The following test cases are provided for the QTF verification with HydroSTAR software :cite:p:`HYDROSTAR`.

***************
8b_QTF_Cylinder
***************

Full body mesh with lid panels, CoG :math:`(0,0,0)`, deep water, wave direction :math:`0^{\circ}`, the difference-frequency QTF DUOK+HASBO. The results are shown in the density plot, :numref:`fig:QTFM_Cylinder_surge`, :numref:`fig:QTFM_Cylinder_heave` and :numref:`fig:QTFM_Cylinder_pitch`, and in the off-diagonal line plot, :numref:`fig:QTFM_diag_Cylinder_surge`, :numref:`fig:QTFM_diag_Cylinder_heave` and :numref:`fig:QTFM_diag_Cylinder_pitch`. The mesh used was obtained using GMSH :cite:p:`GMSH` and is shown in :numref:`fig:meshesCylinder`.

.. _`fig:QTFM_Cylinder_surge`:
.. figure:: figures/Cylinder/QTFsurge.svg
   :align: center

   Density plots of the normalized surge difference frequency QTF magnitude (without the free-surface integrals) for the floating Cylinder (test case **8b_QTF_Cylinder**. HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_Cylinder_heave`:
.. figure:: figures/Cylinder/QTFheave.svg
   :align: center

   Density plots of the normalized heave difference frequency QTF magnitude (without the free-surface integrals) for the floating Cylinder (test case **8b_QTF_Cylinder**. HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_Cylinder_pitch`:
.. figure:: figures/Cylinder/QTFpitch.svg
   :align: center

   Density plots of the normalized pitch difference frequency QTF magnitude (without the free-surface integrals) for the floating Cylinder (test case **8b_QTF_Cylinder**. HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_diag_Cylinder_surge`:
.. figure:: figures/Cylinder/QTFsurge_diag.svg
   :align: center

   Comparison of the surge off-diagonal difference frequency QTF for the Cylinder (test case **8b_QTF_Cylinder**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_Cylinder_heave`:
.. figure:: figures/Cylinder/QTFheave_diag.svg
   :align: center

   Comparison of the heave off-diagonal difference frequency QTF for the Cylinder (test case **8b_QTF_Cylinder**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_Cylinder_pitch`:
.. figure:: figures/Cylinder/QTFpitch_diag.svg
   :align: center

   Comparison of the pitch off-diagonal difference frequency QTF for the Cylinder (test case **8b_QTF_Cylinder**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

*************************
9_QTF_OC4_Semisubmersible
*************************

Full body mesh with lid panels, CoG :math:`(0,0,0)`, water depth 200 m, wave direction :math:`0^{\circ}` and :math:`30^{\circ}`, bi-directional QTF, the difference-frequency QTF DUOK+HASBO. The results are shown in the density plot, :numref:`fig:QTFM_OC4_surge`, :numref:`fig:QTFM_OC4_heave` and :numref:`fig:QTFM_OC4_pitch`, and in the off-diagonal line plot, :numref:`fig:QTFM_diag_OC4_surge`, :numref:`fig:QTFM_diag_OC4_heave` and :numref:`fig:QTFM_diag_OC4_pitch`, of the bi-directional QTF :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`. The mesh used was obtained using GMSH :cite:p:`GMSH` and is shown in :numref:`fig:meshesOC4`.

.. _`fig:meshesOC4`:
.. figure:: figures/OC4/bodymesh.svg
   :align: center

   Body boundary mesh for for the OC4-platform used for test case **9_QTF_OC4_Semisubmersible**.

.. _`fig:QTFM_OC4_surge`:
.. figure:: figures/OC4/QTFM_Surge_beta030.svg
   :align: center

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, surge difference frequency QTF magnitude (without the free-surface integrals) for the floating OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_OC4_heave`:
.. figure:: figures/OC4/QTFM_Heave_beta030.svg
   :align: center

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, heave difference frequency QTF magnitude (without the free-surface integrals) for the floating OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_OC4_pitch`:
.. figure:: figures/OC4/QTFM_Pitch_beta030.svg
   :align: center

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, pitch difference frequency QTF magnitude (without the free-surface integrals) for the floating OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_diag_OC4_surge`:
.. figure:: figures/OC4/QTFM_Surge_beta030_diag.svg
   :align: center

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, surge difference frequency QTF for the OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_OC4_heave`:
.. figure:: figures/OC4/QTFM_Heave_beta030_diag.svg
   :align: center

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, heave difference frequency QTF for the OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_OC4_pitch`:
.. figure:: figures/OC4/QTFM_Pitch_beta030_diag.svg
   :align: center

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, pitch difference frequency QTF for the OC4-semisubmersible platform (test case **9_QTF_OC4_Semisubmersible**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

****************
10a_QTF_SOFTWIND
****************

Half symmetric body mesh with lid panels, CoG :math:`(0,0,-71.56)`, water depth 200 m, wave direction :math:`0^{\circ}` and :math:`30^{\circ}`, bi-directional QTF, the difference-frequency QTF DUOK+HASBO. The results are shown in the density plot, :numref:`fig:QTFM_SOFTWIND_surge`, :numref:`fig:QTFM_SOFTWIND_heave` and :numref:`fig:QTFM_SOFTWIND_pitch`, and in the off-diagonal line plot, :numref:`fig:QTFM_diag_softwind_surge`, :numref:`fig:QTFM_diag_softwind_heave` and :numref:`fig:QTFM_diag_softwind_pitch`, of the bi-directional QTF :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`. The mesh used was obtained using GMSH :cite:p:`GMSH` and is shown in :numref:`fig:meshesSoftwind_body`.

.. _`fig:meshesSoftwind_body`:
.. figure:: figures/Softwind/bodymesh.svg
   :align: center

   Body boundary mesh for the SOFTWIND platform, used in test cases **10a_QTF_SOFTWIND** and  **10b_QTF_SOFTWIND_FS**

.. _`fig:QTFM_SOFTWIND_surge`:
.. figure:: figures/Softwind/QTFM_Surge_beta030.svg
   :align: center

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, surge difference frequency QTF magnitude (without the free-surface integrals) for the floating SOFTWIND platform (test case **10a_QTF_SOFTWIND**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_SOFTWIND_heave`:
.. figure:: figures/Softwind/QTFM_Heave_beta030.svg
   :align: center

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, heave difference frequency QTF magnitude (without the free-surface integrals) for the floating SOFTWIND platform (test case **10a_QTF_SOFTWIND**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_SOFTWIND_pitch`:
.. figure:: figures/Softwind/QTFM_Pitch_beta030.svg
   :align: center

   Density plots of the normalized bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, pitch difference frequency QTF magnitude (without the free-surface integrals) for the floating SOFTWIND platform (test case **10a_QTF_SOFTWIND**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_diag_softwind_surge`:
.. figure:: figures/Softwind/QTFM_Surge_beta030_diag.svg
   :align: center

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, surge difference frequency QTF for the SOFTWIND platform (test case **10a_QTF_SOFTWIND**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_softwind_heave`:
.. figure:: figures/Softwind/QTFM_Heave_beta030_diag.svg
   :align: center

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, heave difference frequency QTF for the SOFTWIND platform (test case **10a_QTF_SOFTWIND**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_softwind_pitch`:
.. figure:: figures/Softwind/QTFM_Pitch_beta030_diag.svg
   :align: center

   Comparison of the off-diagonal bi-directional, :math:`(\beta_1,\beta_2)=(0^{\circ},30^{\circ})`, pitch difference frequency QTF for the SOFTWIND platform (test case **10a_QTF_SOFTWIND**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

*******************
10b_QTF_SOFTWIND_FS
*******************

Half symmetric body mesh without lid panels, half symmetric free-surface mesh, CoG :math:`(0,0,-71.56)`, water depth 200 m, wave direction :math:`0^{\circ}`, the sum-frequency total QTF DUOK+HASBO+HASFS+ASYMP.
The results are shown in the density plot, :numref:`fig:QTFP_SOFTWIND_surge`, :numref:`fig:QTFP_SOFTWIND_heave` and :numref:`fig:QTFP_SOFTWIND_pitch` and in the off-diagonal line plot, :numref:`fig:QTFP_SOFTWIND_DIAG_surge`, :numref:`fig:QTFP_SOFTWIND_DIAG_heave` and :numref:`fig:QTFP_SOFTWIND_DIAG_pitch`.
The mesh used was obtained using GMSH :cite:p:`GMSH` and is shown in :numref:`fig:meshesSoftwind_body` (body mesh) and :numref:`fig:meshesSoftwind_FS` (free surface mesh).

.. _`fig:meshesSoftwind_FS`:
.. figure:: figures/Softwind/FSmesh.svg
   :align: center

   Free surface mesh for the SOFTWIND platform, used in test case **10b_QTF_SOFTWIND_FS**

.. _`fig:QTFP_SOFTWIND_surge`:
.. figure:: figures/Softwind/QTFP_Surge_beta00.svg
   :align: center

   Density plots of the normalized surge sum-frequency full QTF magnitude (including the free-surface integrals) for the floating SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference in the right column.

.. _`fig:QTFP_SOFTWIND_heave`:
.. figure:: figures/Softwind/QTFP_Heave_beta00.svg
   :align: center

   Density plots of the normalized heave sum-frequency full QTF magnitude (including the free-surface integrals) for the floating SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference in the right column.

.. _`fig:QTFP_SOFTWIND_pitch`:
.. figure:: figures/Softwind/QTFP_Pitch_beta00.svg
   :align: center

   Density plots of the normalized pitch sum-frequency full QTF magnitude (including the free-surface integrals) for the floating SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference in the right column.

.. _`fig:QTFP_SOFTWIND_DIAG_surge`:
.. figure:: figures/Softwind/QTFP_Surge_beta00_diag.svg
   :align: center

   Comparison of the off-diagonal surge sum-frequency full QTF for SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFP_SOFTWIND_DIAG_heave`:
.. figure:: figures/Softwind/QTFP_Heave_beta00_diag.svg
   :align: center

   Comparison of the off-diagonal heave sum-frequency full QTF for SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFP_SOFTWIND_DIAG_pitch`:
.. figure:: figures/Softwind/QTFP_Pitch_beta00_diag.svg
   :align: center

   Comparison of the off-diagonal pitch sum-frequency full QTF for SOFTWIND platform (test case **10b_QTF_SOFTWIND_FS**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

*****************
11_QTF_OC3_Hywind
*****************

Full body mesh with lid panels, CoG :math:`(0,0,0)`, water depth 320 m, wave direction :math:`0^{\circ}`, Nemoh1 uses GMRES solver, the difference-frequency QTF DUOK+HASBO. The results are shown in the density plot, :numref:`fig:QTFM_OC3_HYWIND_surge`, :numref:`fig:QTFM_OC3_HYWIND_heave` and :numref:`fig:QTFM_OC3_HYWIND_pitch`, and in the off-diagonal line plot, :numref:`fig:QTFM_diag_OC3_HYWIND_surge`, :numref:`fig:QTFM_diag_OC3_HYWIND_heave` and :numref:`fig:QTFM_diag_OC3_HYWIND_pitch`, of the difference-frequency QTF. The mesh used was obtained using GMSH :cite:p:`GMSH` and is shown in :numref:`fig:meshesHYWIND`.

.. _`fig:meshesHYWIND`:
.. figure:: figures/OC3_HYWIND/bodyMesh.svg
   :align: center

   Body boundary mesh for OC3-HYWIND platform, test case **11_QTF_OC3_Hywind**.

.. _`fig:QTFM_OC3_HYWIND_surge`:
.. figure:: figures/OC3_HYWIND/QTFM_Surge.svg
   :align: center

   Density plots of the normalized surge difference frequency QTF magnitude (without the free-surface integrals) for the floating OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_OC3_HYWIND_heave`:
.. figure:: figures/OC3_HYWIND/QTFM_Heave.svg
   :align: center

   Density plots of the normalized heave difference frequency QTF magnitude (without the free-surface integrals) for the floating OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_OC3_HYWIND_pitch`:
.. figure:: figures/OC3_HYWIND/QTFM_Pitch.svg
   :align: center

   Density plots of the normalized pitch difference frequency QTF magnitude (without the free-surface integrals) for the floating OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**). HydroSTAR results are on the left column, Nemoh results are on the middle column and the difference on the right column.

.. _`fig:QTFM_diag_OC3_HYWIND_surge`:
.. figure:: figures/OC3_HYWIND/QTFM_Surge_diag.svg
   :align: center

   Comparison of the off-diagonal surge difference frequency QTF for the OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_OC3_HYWIND_heave`:
.. figure:: figures/OC3_HYWIND/QTFM_Heave_diag.svg
   :align: center

   Comparison of the off-diagonal heave difference frequency QTF for the OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

.. _`fig:QTFM_diag_OC3_HYWIND_pitch`:
.. figure:: figures/OC3_HYWIND/QTFM_Pitch_diag.svg
   :align: center

   Comparison of the off-diagonal pitch difference frequency QTF for the OC3-HYWIND platform (test case **11_QTF_OC3_Hywind**) between HydroSTAR, real part (blue, solid-line), imaginary part (blue, dashed-dot line) and Nemoh, real part (red, dashed-line), imaginary part (red, dotted-line).

Full description of the QTF test-cases results is reported in :cite:t:`Kurnia22_JH,Kurnia22`. Note that the QTF comparisons between Nemoh and HydroSTAR for the bidirectional case are in good agreement only if the direction is switched, in Nemoh :math:`\beta=(\beta_1,\beta_2)` and in HydroSTAR :math:`\beta=(\beta_2,\beta_1)`; further investigation regarding this is needed. The imaginary part of QTFs have also a difference sign between Nemoh and HydroSTAR that may be due to different conventions of the incident potential.
