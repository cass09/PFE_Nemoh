bindir=../../bin/

.PHONY: run clean mesh

mesh:
			$(bindir)/mesh .

# Build and run Sphere reference test case.
run:
			$(bindir)/preProc .
			$(bindir)/hydrosCal .
			$(bindir)/solver .
			$(bindir)/postProc .
			# cp -r reference_results/Aquaplus results/
			# cp reference_results/cylinder.lay results/
			
# Delete results of the computations
clean:
			rm -rf mesh
			rm -rf results/
			rm -rf Normalvelocities.dat
			rm -rf logfile.txt
			rm -rf Motion/
			rm -rf Mechanics/
			rm -rf input_IT.dat
