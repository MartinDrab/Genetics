default:
	cd gassm2; make; cd ..
	cd rcorrect; make; cd ..
	cd vcffix; make; cd ..
	cd vcfcmp; make; cd ..

clean:
	cd gassm2; make clean; cd ..
	cd rcorrect; make clean; cd ..
	cd vcffix; make clean; cd ..
	cd vcfcmp; make clean; cd ..


