
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

		auto fs = std::fstream(aFileName, std::fstream::in);
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
			CUpdatedRefSeq orig(args[1], refseq);
			CUpdatedRefSeq res(args[2], refseq);
			CVCFCompareStatistics stats;
			orig.compare(res, stats);
			std::cout << "Done" << std::endl;
			for (const auto & diff : stats.Differences) {
				std::cout << "POS:      " << diff.Pos << std::endl;
				std::cout << "ORIGINAL: " << diff.AltOriginal << std::endl;
				std::cout << "OTHER:    " << diff.AltOther << std::endl;
				std::cout << "REF:      " << diff.Ref << std::endl;
				std::cout << std::endl;
			}

			std::cout << "Insertions:   " << stats.Insertions << " / " << orig.Insertions() << std::endl;
			std::cout << "Deletions:    " << stats.Deletions << " / " << orig.Deletions() << std::endl;
			std::cout << "Replaces:     " << stats.Replaces << " / " << orig.Replaces << std::endl;
			std::cout << "SNPs:         " << stats.SNPs << " / " << orig.SNPs() << std::endl;
			std::cout << "Undiscovered: " << stats.Undiscovered << std::endl;
		}
	} else print_usage();

	return 0;
}
