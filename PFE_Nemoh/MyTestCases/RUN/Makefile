bindir=../../bin/

.PHONY: run clean run_hydro

run_hydro:
			$(bindir)/preProc .
			$(bindir)/hydrosCal .
			$(bindir)/postProc .
# Build and run Sphere reference test case.
run:
			$(bindir)/preProc .
			$(bindir)/hydrosCal .
			$(bindir)/solver .
			$(bindir)/postProc .
			# cp -r reference_results/Aquaplus results/
			# cp reference_results/cylinder.lay results/
			rm -rf mesh
			rm -rf Normalvelocities.dat
			
# Delete results of the computations
clean:
			rm -rf mesh
			rm -rf results/
			rm -rf Normalvelocities.dat
			rm -rf logfile.txt
			rm -rf Motion/
			rm -rf Mechanics/
