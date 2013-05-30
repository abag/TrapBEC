default: run
run:
	make -C src 
clean:
	make -C src clean
cleann:
	make -C src cleann
pristine: cleann
	rm -rf data/*
	rm -f *.eps *.png
	rm -f *.mat
	rm -f *.vtk
	rm -f *.vdf
	rm -f STOP ERROR
linkx:
	@for file in src/*.x; \
	do [ -e "`basename $$file`" ] || ln -s $$file .; \
	done
