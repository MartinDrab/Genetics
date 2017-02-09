
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>
#include <sstream>
#include "tinydir.h"



int process_directory(const std::string & aName)
{
	int ret = 0;
	tinydir_dir entry;
	std::map<unsigned long long, std::string> vcfMap;
	std::vector<std::string> fileNames;

	ret = tinydir_open(&entry, aName.data());
	std::cerr << "tinydir_open(): " << ret << std::endl;
	if (ret == 0) {
		tinydir_file fileEntry;

		do {
			ret = tinydir_readfile(&entry, &fileEntry);
			std::cerr << "tinydir_readfile(): " << ret << std::endl;
			std::string fn = std::string(fileEntry.name);

			if (fn.size() > 4 && fn.substr(fn.size() - 4) == ".vcf")
				fileNames.push_back(fn);

			ret = tinydir_next(&entry);
		} while (ret == 0);

		for (auto fn : fileNames) {
			std::string n = aName;
			std::ifstream fs(n.append("/").append(fn));
			std::string strPos = fn.substr(0, fn.size() - 4);
			std::string line;

			std::cerr << fn << std::endl;
			unsigned long long startPos = std::strtoull(strPos.data(), NULL, 10);
			while (!fs.eof()) {
				std::getline(fs, line);
				if (line.size() > 1) {
					if (line[0] == '#')
						continue;

					auto startIt = line.find_first_of('\t');
					++startIt;
					auto endIt = line.find_first_of('\t', startIt);

					unsigned long long recPos = std::strtoull(line.substr(startIt, endIt - startIt).data(), NULL, 10);
					std::stringstream ss;
					ss << "1\t";
					ss << (startPos + recPos);
					ss << line.substr(endIt);
					vcfMap.insert(std::make_pair(startPos + recPos, ss.str()));
				} else break;
			}
		}

		tinydir_close(&entry);
	}

	std::cout << "##fileformat=VCFv4.1" << std::endl;
	std::cout << "##fileDate = 20160525" << std::endl;
	std::cout << "##source=GASSMV2" << std::endl;
	std::cout << "##reference=1" << std::endl;
	std::cout << "##phasing=partial" << std::endl;
	std::cout << "##contig=<ID=1,length=249250621>" << std::endl;
	std::cout << "##contig=<ID=2,length=243199373>" << std::endl;
	std::cout << "##contig=<ID=3,length=198022430>" << std::endl;
	std::cout << "##contig=<ID=4,length=191154276>" << std::endl;
	std::cout << "##contig=<ID=5,length=180915260>" << std::endl;
	std::cout << "##contig=<ID=6,length=171115067>" << std::endl;
	std::cout << "##contig=<ID=7,length=159138663>" << std::endl;
	std::cout << "##contig=<ID=8,length=146364022>" << std::endl;
	std::cout << "##contig=<ID=9,length=141213431>" << std::endl;
	std::cout << "##contig=<ID=10,length=135534747>" << std::endl;
	std::cout << "##contig=<ID=11,length=135006516>" << std::endl;
	std::cout << "##contig=<ID=12,length=133851895>" << std::endl;
	std::cout << "##contig=<ID=13,length=115169878>" << std::endl;
	std::cout << "##contig=<ID=14,length=107349540>" << std::endl;
	std::cout << "##contig=<ID=15,length=102531392>" << std::endl;
	std::cout << "##contig=<ID=16,length=90354753>" << std::endl;
	std::cout << "##contig=<ID=17,length=81195210>" << std::endl;
	std::cout << "##contig=<ID=18,length=78077248>" << std::endl;
	std::cout << "##contig=<ID=19,length=59128983>" << std::endl;
	std::cout << "##contig=<ID=20,length=63025520>" << std::endl;
	std::cout << "##contig=<ID=21,length=48129895>" << std::endl;
	std::cout << "##contig=<ID=22,length=51304566>" << std::endl;
	std::cout << "##contig=<ID=X,length=155270560>" << std::endl;
	std::cout << "##contig=<ID=Y,length=59373566>" << std::endl;
	std::cout << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" << std::endl;
	std::cout << "##FORMAT=<ID=PS,Number=1,Type=String,Description=\"Phase number\">" << std::endl;
	std::cout << "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">" << std::endl;
	std::cout << "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	DNAcko" << std::endl;
	for (auto & p : vcfMap)
		std::cout << p.second << std::endl;

	return ret;
}



int main(int argc, char **argv)
{
	std::vector<std::string> args(argv + 1, argv + argc);

	process_directory(args[0]);

	return 0;
}