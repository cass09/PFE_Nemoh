################
#  My Test cases  #
################
# TestCases directory
testdir=.

#NEMOH1 testcases
.PHONY: run_1_cylinder clean_1_cylinder
run_1_cylinder: 
	$(MAKE) -C $(testdir)/1_Cylinder/ run

clean_1_cylinder:
	$(MAKE) -C $(testdir)/1_Cylinder/ clean

.PHONY: run_my_cylinder clean_my_cylinder
run_my_cylinder: 
	$(MAKE) -C $(testdir)/my_Cylinder/ run

clean_my_cylinder:
	$(MAKE) -C $(testdir)/my_Cylinder/ clean



.PHONY: clean_all_tests
clean_all_tests: clean_1_cylinder clean_my_cylinder


