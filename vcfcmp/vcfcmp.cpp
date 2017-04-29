
#include <omp.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "vcf.h"



void print_usage(void)
{
	std::cout << "Usage: vcfcmp <RefSeqFile> <VCFFile1> <VCFFile2>" << std::endl;

	return;
}



bool load_reference(const std::string & aFileName, std::string & Reference)
{
	bool ret = true;

	Reference.clear();
	try {
		std::string line;

		std::ifstream fs;
		fs.open(aFileName);
		try {
			while (!fs.eof()) {
				std::getline(fs, line);
				if (line.size() > 0 && line.at(0) == '>')
					continue;

				Reference = Reference.append(line);
			}
		} catch (int err) {
			printf("ERROR: Failed to read a line from the file: %i\n", err);
		}
	} catch (int err) {
		printf("ERROR: failed to open the reference sequence file: %i\n", err);
		ret = false;
	}

	return ret;
}


int main(int argc, char *argv[])
{
	std::vector<std::string> args(argv + 1, argv + argc);

	if (args.size() == 3) {
		std::string refseq;

		std::cerr << "Loading reference..." << std::endl;
		if (load_reference(args[0], refseq)) {
			std::cerr << "Loading the truth set VCF..." << std::endl;
			CVCFFile orig(args[1]);
			std::cerr << "Loading our VCF..." << std::endl;
			CVCFFile res(args[2]);
			CVCFGraph graph;

			graph.AddReference(refseq);
			std::cerr << "Adding ref parts of our set..." << std::endl;
			for (const auto & rec : res.Records)
				graph.AddVCFRecordRef(rec, refseq);

			std::cerr << "Adding ref parts of the truth set..." << std::endl;
			for (const auto & rec : orig.Records)
				graph.AddVCFRecordRef(rec, refseq);

			std::cerr << "Adding alt parts of our set..." << std::endl;
			for (const auto & rec : res.Records)
				graph.AddVCFRecordAlt(rec, refseq);

			std::cerr << "Comparing the sets..." << std::endl;
			int numProcs = omp_get_num_procs();
			size_t *insertionArray = new size_t[numProcs*4];
			memset(insertionArray, 0, 4*sizeof(size_t)*numProcs);
			size_t *deletionArray = insertionArray + numProcs;
			size_t *replaceArray = deletionArray + numProcs;
			size_t *SNPArray = replaceArray + numProcs;
			int index = 0;

			omp_lock_t printLock;

			omp_init_lock(&printLock);
#pragma omp parallel for shared(graph, refseq, orig, insertionArray, deletionArray, replaceArray, SNPArray)
			for (index = 0; index < (int)orig.Records.size(); ++index) {
				const CVCFRecord & rec = orig.Records[index];
				int tid = omp_get_thread_num();

				if (graph.CheckVCFRecord(rec, refseq)) {
					switch (rec.Type) {
						case vcfrtInsertion: ++insertionArray[tid]; break;
						case vcfrtDeletion: ++deletionArray[tid]; break;
						case vcfrtReplace: ++replaceArray[tid]; break;
						case vcfrtSNP: ++SNPArray[tid]; break;
						default: throw std::exception(); break;
					}
				} else {
					omp_set_lock(&printLock);
					std::cout << rec.Pos << "\t" << rec.Ref << "\t"  << rec.Alt <<  "\t" << rec.SSW << std::endl;
					omp_unset_lock(&printLock);
				}
			}

			omp_destroy_lock(&printLock);

			size_t insertions = 0;
			size_t deletions = 0;
			size_t replaces = 0;
			size_t SNPs = 0;
			for (int i = 0; i < numProcs; ++i) {
				insertions += insertionArray[i];
				deletions += deletionArray[i];
				replaces += replaceArray[i];
				SNPs += SNPArray[i];
			}

			std::cout << "Insertions: " << insertions << " / " << orig.Insertions << " (" << ((insertions + 1)*100 / (orig.Insertions + 1)) << " %)" <<  std::endl;
			std::cout << "Deletions:  " << deletions << " / " << orig.Deletions << " (" << ((deletions + 1)*100 / (orig.Deletions + 1)) << " %)" << std::endl;
			std::cout << "Replaces:   " << replaces << " / " << orig.Replaces << " (" << ((replaces + 1)*100 / (orig.Replaces + 1)) << " %)" << std::endl;
			std::cout << "SNPs:       " << SNPs << " / " << orig.SNPs << " (" << ((SNPs + 1)*100 / (orig.SNPs + 1)) << " %)" << std::endl;
			
			delete[] insertionArray;
		}
	} else if (args.size() == 2) {
		std::cerr << "Loading the truth set VCF..." << std::endl;
		CVCFFile orig(args[0]);
		std::cerr << "Loading our VCF..." << std::endl;
		CVCFFile res(args[1]);
		std::cerr << "Comparing the sets..." << std::endl;
		auto truthIt = orig.Records.cbegin();
		auto testIt = res.Records.cbegin();
		std::vector<CVCFRecord> TPs;
		std::vector<CVCFRecord> FPs;
		std::vector<CVCFRecord> FNs;

		while (truthIt != orig.Records.cend() && testIt != res.Records.cend()) {
			if (truthIt->Pos < testIt->Pos) {
				FNs.push_back(*truthIt);
				++truthIt;
			} else if (truthIt->Pos > testIt->Pos) {
				FPs.push_back(*testIt);
				++testIt;
			} else {
				if (truthIt->Ref == testIt->Ref &&
					truthIt->Alt == testIt->Alt) {
					TPs.push_back(*truthIt);
				} else {
					FNs.push_back(*truthIt);
					FPs.push_back(*testIt);
				}

				++truthIt;
				++testIt;
			}
		}

		std::cerr << "Finishing false negatives..." << std::endl;
		while (truthIt != orig.Records.cend()) {
			FNs.push_back(*truthIt);
			++truthIt;
		}

		std::cerr << "Finishing false positives..." << std::endl;
		while (testIt != res.Records.cend()) {
			FPs.push_back(*testIt);
			++testIt;
		}

		std::cerr << "TPs: " << TPs.size() << std::endl;
		std::cerr << "FNs: " << FNs.size() << std::endl;
		std::cerr << "FPs: " << FPs.size() << std::endl;
		
		std::ofstream FNsFile("FNs.vcf");
		for (auto & r : FNs) {
			if (r.Type == vcfrtSNP)
				r.write(FNsFile);
		}

		FNsFile.close();

		std::ofstream FPsFile("FPs.vcf");
		for (auto & r : FPs) {
			if (r.Type == vcfrtSNP)
				r.write(FPsFile);
		}

		FPsFile.close();
	} else if (args.size() == 1) {
		std::cerr << "Loading a VCF file..." << std::endl;
		CVCFFile vcf(args[0]);
		std::map<std::size_t, std::size_t> insertionCounts_;
		std::map<std::size_t, std::size_t> deletionCounts_;

		for (auto & r : vcf.Records) {
			switch (r.Type) {
				case vcfrtDeletion:
					deletionCounts_.insert(std::make_pair(r.Ref.size(), 0)).first->second++;
					break;
				case vcfrtInsertion:
					insertionCounts_.insert(std::make_pair(r.Alt.size(), 0)).first->second++;
					break;
			}
		}

		std::cerr << "Total SNPs:        " << vcf.SNPs << std::endl;
		std::cerr << "Total Replaces:    " << vcf.Replaces << std::endl;
		std::cerr << "Total insertions:  " << vcf.Insertions << std::endl;
		for (auto & p : insertionCounts_)
			std::cerr << "  " << p.first << ", " << p.second << std::endl;

		std::cerr << "Total deletions:   " << vcf.Deletions << std::endl;
		for (auto & p : deletionCounts_)
			std::cerr << "  " << p.first << ", " << p.second << std::endl;
	} else print_usage();
	
	return 0;
}
