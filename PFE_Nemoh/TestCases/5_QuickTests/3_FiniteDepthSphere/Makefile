SHELL=/bin/bash

bindir=../../../bin/
logfile=log

.PHONY: run silent_run test clean

# Build and run Sphere reference test case.
run:
	$(bindir)/preProc .
	$(bindir)/solver .
	$(bindir)/postProc .

silent_run:
	@test -f $(logfile) || echo '' > $(logfile)
	@$(bindir)/preProc  . >> $(logfile)
	@$(bindir)/solver   . >> $(logfile)
	@$(bindir)/postProc . >> $(logfile)

test: silent_run
	@echo "Forces.dat"
	@cat <(echo -n "    Present result") <(tr -d "\n\r" < results/Forces.dat)
	@cat <(echo -e "  \n  Reference ")  reference_results/Forces.dat
	@echo "Line 20 of freesurface.00001.dat"
	@cat <(echo -n "    Present result") <(head -n 20 results/freesurface.00001.dat | tail -1)
	@cat <(echo -n "    Reference     ") <(head -n 20 "reference_results/freesurface.    1.dat" | tail -1)

# Delete results of the computations
clean:
	rm -f $(logfile)
	rm -rf mesh
	rm -rf results/
	rm -rf Normalvelocities.dat
	rm -rf logfile.txt
