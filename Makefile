all: 1 2 3

1:
	@echo "Building LIBRCORRECT"
	cd librcorrect; make; cd ..

2:
	echo "Building LIBKMER"
	cd libkmer; make; cd ..


3: 1 2
	
	echo "Building GASSMv2"
	cd gassm2; make; cd ..



clean:
	cd gassm2; make clean; cd ..
	cd libkmer; make clean; cd ..
	cd librcorrect; make clean; cd ..
	cd lib; rm *; cd ..