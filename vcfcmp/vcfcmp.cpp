
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

		if (load_reference(args[0], refseq)) {
			CVCFFile orig(args[1]);
			CVCFFile res(args[2]);
			CVCFGraph graph;

			graph.AddReference(refseq);
			for (const auto & rec : res.Records)
				graph.AddVCFRecord(rec, refseq);

			int numProcs = omp_get_num_procs();
			size_t *insertionArray = new size_t[numProcs*4];
			memset(insertionArray, 0, 4*sizeof(size_t)*numProcs);
			size_t *deletionArray = insertionArray + numProcs;
			size_t *replaceArray = deletionArray + numProcs;
			size_t *SNPArray = replaceArray + numProcs;
			int index = 0;
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
				}
			}

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
	} else print_usage();

	return 0;
}
