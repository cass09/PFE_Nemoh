################
#  Test cases  #
################
# TestCases directory
testdir=.

#NEMOH1 testcases
.PHONY: run_1_cylinder clean_1_cylinder
run_1_cylinder: 
	$(MAKE) -C $(testdir)/1_Cylinder/ run

clean_1_cylinder:
	$(MAKE) -C $(testdir)/1_Cylinder/ clean

.PHONY: run_2_2Bodies clean_2_2Bodies
run_2_2Bodies: 
	$(MAKE) -C $(testdir)/2_2Bodies/ run

clean_2_2Bodies:
	$(MAKE) -C $(testdir)/2_2Bodies/ clean

.PHONY: run_3_nonsymmetrical clean_3_nonsymmetrical
run_3_nonsymmetrical: 
	$(MAKE) -C $(testdir)/3_NonSymmetrical/ run

clean_3_nonsymmetrical:
	$(MAKE) -C $(testdir)/3_NonSymmetrical/ clean

.PHONY: run_4_Postprocessing clean_4_Postprocessing
run_4_Postprocessing: 
	$(MAKE) -C $(testdir)/4_Postprocessing/ run

clean_4_Postprocessing:
	$(MAKE) -C $(testdir)/4_Postprocessing/ clean


.PHONY: run_5_quicktest clean_5_quicktest
run_5_quicktest: 
	@echo ""
	@echo "Sphere"
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/1_Sphere/                     test
	@echo ""
	@echo "Sphere using y-symmetry"
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/2_SymmetricSphere/            test
	@echo ""
	@echo "Sphere in finite depth"
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/3_FiniteDepthSphere/          test
	@echo ""
	@echo "Sphere in finite depth using y-symmetry"
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/4_SymmetricFiniteDepthSphere/ test
	@echo ""
	@echo "Alien sphere"
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/5_AlienSphere/                test

clean_5_quicktest:
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/1_Sphere/                     clean
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/2_SymmetricSphere/            clean
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/3_FiniteDepthSphere/          clean
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/4_SymmetricFiniteDepthSphere/ clean
	@$(MAKE) --silent -C $(testdir)/5_QuickTests/5_AlienSphere/                clean

.PHONY: run_6_box_coarsemesh clean_6_box_coarsemesh
run_6_box_coarsemesh:
	$(MAKE) -C $(testdir)/6_box_coarsemesh/ run

clean_6_box_coarsemesh:
	$(MAKE) -C $(testdir)/6_box_coarsemesh/ clean

.PHONY: run_7_Solvers_Check_OC3 clean_Solvers_Check_OC3
run_7_Solvers_Check_OC3:
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh978_GaussElim/ 	run
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh978_LUdecomp/ 		run
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh978_GMRES/ 		run
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh4842_GaussElim/ 	run
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh4842_LUdecomp/ 	run
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh4842_GMRES/ 		run
	
clean_7_Solvers_Check_OC3:
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh978_GaussElim/ 	clean
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh978_LUdecomp/ 		clean
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh978_GMRES/ 		clean
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh4842_GaussElim/ 	clean
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh4842_LUdecomp/ 	clean
	$(MAKE) -C $(testdir)/7_Solvers_Check_OC3/mesh4842_GMRES/ 		clean

.PHONY: run_8a_Cylinder_irregfreq clean_8a_Cylinder_irregfreq
run_8a_Cylinder_irregfreq:
	$(MAKE) -C $(testdir)/8a_Cylinder_irregfreq/mesh634/ run
	$(MAKE) -C $(testdir)/8a_Cylinder_irregfreq/mesh634_nolid/ run

clean_8a_Cylinder_irregfreq:
	$(MAKE) -C $(testdir)/8a_Cylinder_irregfreq/mesh634/ clean
	$(MAKE) -C $(testdir)/8a_Cylinder_irregfreq/mesh634_nolid/ clean

.PHONY: clean_all_testsNEMOH1
clean_all_testsNEMOH1:	clean_1_cylinder clean_2_2Bodies clean_3_nonsymmetrical \
						clean_4_Postprocessing clean_5_quicktest clean_6_box_coarsemesh \
						clean_7_Solvers_Check_OC3 clean_8a_Cylinder_irregfreq
#NEMOH2 testcases	
.PHONY: run_8b_QTF_Cylinder clean_8b_QTF_Cylinder
run_8b_QTF_Cylinder:
	$(MAKE) -C $(testdir)/8b_QTF_Cylinder/mesh634_floating/ run

clean_8b_QTF_Cylinder:
	$(MAKE) -C $(testdir)/8b_QTF_Cylinder/mesh634_floating/ clean

.PHONY: run_9_oc4_semisub clean_9_oc4_semisub
run_9_oc4_semisub: 
	$(MAKE) -C $(testdir)/9_QTF_OC4_Semisubmersible/mesh2200_Floating/ run

clean_9_oc4_semisub:
	$(MAKE) -C $(testdir)/9_QTF_OC4_Semisubmersible/mesh2200_Floating/ clean
	$(MAKE) -C $(testdir)/9_QTF_OC4_Semisubmersible/mesh980_Floating/ clean

.PHONY: run_10a_softwind clean_10a_softwind
run_10a_softwind: 
	$(MAKE) -C $(testdir)/10_QTF_SOFTWIND/mesh_1872half_QTF_Floating/ run

clean_10a_softwind:
	$(MAKE) -C $(testdir)/10_QTF_SOFTWIND/mesh_1872half_QTF_Floating/ clean
	$(MAKE) -C $(testdir)/10_QTF_SOFTWIND/mesh_610half_QTF_Floating/ clean

.PHONY: run_10b_softwind_FS clean_10b_softwind_FS
run_10b_softwind_FS: 
	$(MAKE) -C $(testdir)/10_QTF_SOFTWIND/mesh_1872half_nolid_QTF_Floating_FS/ run

clean_10b_softwind_FS:
	$(MAKE) -C $(testdir)/10_QTF_SOFTWIND/mesh_1872half_nolid_QTF_Floating_FS/ clean
	$(MAKE) -C $(testdir)/10_QTF_SOFTWIND/mesh_610half_nolid_QTF_Floating_FS/ clean

.PHONY: run_11_QTF_OC3_Hywind clean_11_QTF_OC3_Hywind
run_11_QTF_OC3_Hywind: 
	$(MAKE) -C $(testdir)/11_QTF_OC3_Hywind/mesh4842_Floating/ run

clean_11_QTF_OC3_Hywind:
	$(MAKE) -C $(testdir)/11_QTF_OC3_Hywind/mesh4842_Floating/ clean
	$(MAKE) -C $(testdir)/11_QTF_OC3_Hywind/mesh978_Floating/ clean

.PHONY: clean_all_testsNEMOH2
clean_all_testsNEMOH2: clean_8b_QTF_Cylinder clean_9_oc4_semisub clean_10a_softwind \
						clean_10b_softwind_FS clean_11_QTF_OC3_Hywind

.PHONY: clean_all_tests
clean_all_tests: clean_all_testsNEMOH1 clean_all_testsNEMOH2


