bindir=../../bin/

.PHONY: run clean

# Build and run Sphere reference test case.
run:
			$(bindir)/preProc .
			$(bindir)/solver .
			$(bindir)/postProc .
			cp -r reference_results/Aquaplus results/
			cp reference_results/NonSymmetrical-1.lay results/
			cp reference_results/NonSymmetrical-2.lay results/
			
# Delete results of the computations
clean:
			rm -rf mesh
			rm -rf results/
			rm -rf Normalvelocities.dat
			rm -rf logfile.txt
