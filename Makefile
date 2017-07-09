all: 1 2 3

1:
	@echo "Building LIBRCORRECT"
	cd librcorrect; make; cd ..

2:
	echo "Building LIBKMER"
	cd libkmer; make; cd ..


3:
	echo "Building GASSMv2"
	cd gassm2; make; cd ..


3:
	@echo "Building VCFFIX"
	cd vcffix; make; cd ..

4:
	@echo "Building VCFCMP"
	cd vcfcmp; make; cd ..

clean:
	cd gassm2; make clean; cd ..
	cd libkmer; make clean; cd ..
	cd librcorrect; make clean; cd ..
	cd vcffix; make clean; cd ..
	cd vcfcmp; make clean; cd ..