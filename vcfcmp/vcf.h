
#ifndef __VCF_H__
#define __VCF_H__


#include <cstdint>
#include <string>
#include <vector>
#include <cassert>
#include <cstring>
#include <fstream>


enum EVCFRecordType {
	vcfrtUnknown = 0,
	vcfrtReplace = 1,
	vcfrtSNP = 2,
	vcfrtDeletion = 3,
	vcfrtInsertion = 4,
	vcfrtMax = 5,
};

class CVCFRecord {
public:
	CVCFRecord()
		: Type(Type_), RefName(RefName_), Pos(Pos_), Ref(Ref_), Alt(Alt_)
	{ }
	void operator =(const CVCFRecord & aRecord)
	{
		Ref_ = aRecord.Ref;
		Alt_ = aRecord.Alt;
		Pos_ = aRecord.Pos;
		Type_ = aRecord.Type;
		RefName_ = aRecord.RefName;

		return;
	}
	CVCFRecord(const CVCFRecord & aRecord)
		: Pos_(aRecord.Pos), Ref_(aRecord.Ref), Alt_(aRecord.Alt),
		Index_(aRecord.Index), RefName_(aRecord.RefName), Type_(aRecord.Type), Index(Index_),
		Type(Type_), Pos(Pos_), RefName(RefName_), Ref(Ref_), Alt(Alt_)
	{ }
	CVCFRecord(const std::string & aLine, const bool aSecond, const std::size_t aIndex = 0)
		: Type_(vcfrtUnknown), Type(Type_), Ref(Ref_), RefName(RefName_), Pos(Pos_), Alt(Alt_), Index(Index_),
		Index_(aIndex)
	{
		size_t index = aLine.find('\t');
		
		RefName_ = aLine.substr(0, index);
		size_t endIndex = aLine.find('\t', index + 1);
		std::string strPos = aLine.substr(index + 1, endIndex - index - 1);
		Pos_ = (size_t)strtoull(strPos.data(), 0, 10);
		if (Pos_ == 0)
			throw std::exception();

		index = endIndex;
		endIndex = aLine.find('\t', index + 1);
		Ref_ = aLine.substr(index + 1, endIndex - index - 1);
		if (Ref_ == ".") {
			index = endIndex;
			endIndex = aLine.find('\t', index + 1);
			Ref_ = aLine.substr(index + 1, endIndex - index - 1);
		}

		index = endIndex;
		endIndex = aLine.find('\t', index + 1);
		Alt_ = aLine.substr(index + 1, endIndex - index - 1);
		index = Alt_.find(',');
		if (index != std::string::npos) {
			if (aSecond)
				Alt_ = Alt_.substr(index + 1);
			else Alt_ = Alt_.substr(0, index);
		} else if (aSecond)
			Alt_.clear();

		if (Alt_.size() > 0) {
			if (Ref_.size() == 1 && Alt_.size() == 1)
				Type_ = vcfrtSNP;
			else if (Ref_.size() == 1)
				Type_ = vcfrtInsertion;
			else if (Alt_.size() == 1)
				Type_ = vcfrtDeletion;
			else Type_ = vcfrtReplace;
		}

		return;
	}
	const std::string & RefName = RefName_;
	const std::string & Ref = Ref_;
	const std::string & Alt = Alt_;
	const std::size_t & Pos = Pos_;
	const EVCFRecordType & Type = Type_;
	const std::size_t & Index = Index_;
private:
	std::string RefName_;
	std::string Ref_;
	std::string Alt_;
	std::size_t Pos_;
	EVCFRecordType Type_;
	std::size_t Index_;
};


class CRefSeqInterval {
public:
	const char *Sequence() const { return Sequence_; }
	size_t Pos() const { return Pos_; }
	size_t Length() const { return Length_; }
	CRefSeqInterval(const char *aSequence, const size_t aPos, const size_t aLength)
		: Sequence_(aSequence), Pos_(aPos), Length_(aLength) { }
private:
	const char *Sequence_;
	size_t Pos_;
	size_t Length_;
};

struct CVCFDifference {
public:
	size_t Pos;
	std::string Ref;
	std::string AltOriginal;
	std::string AltOther;
	CVCFDifference(const size_t aPos, const std::string & aRef, const std::string & aOrig, const std::string & aOther)
		: Pos(aPos), Ref(aRef), AltOriginal(aOrig), AltOther(aOther) { }
};

struct CVCFCompareStatistics {
public:
	size_t Insertions;
	size_t Deletions;
	size_t SNPs;
	size_t Replaces;
	size_t Undiscovered;
	std::vector<CVCFDifference> Differences;
	CVCFCompareStatistics()
		: Deletions(0), SNPs(0), Replaces(0), Insertions(0), Undiscovered(0) 
	{ }
};

class CUpdatedRefSeq {
public:
	CUpdatedRefSeq() { }
	CUpdatedRefSeq(const std::string & aFileName, const std::string & aReference)
		: Reference_(0), SNPs_(0), Insertions_(0), Deletions_(0), ReferenceLength_(aReference.size()),
		Replaces(Replaces_), Replaces_(0),
		Alternate_(0), AlternateSize_(0)
	{
		std::ifstream fs;
		fs.open(aFileName);
		std::string line;

		while (!fs.eof()) {
			std::getline(fs, line);
			if (line.empty() || (line.size() > 0 && line[0] == '#'))
				continue;

			CVCFRecord vr = CVCFRecord(line, false);
			vcfRecords_.push_back(vr);
			vr = CVCFRecord(line, true);
			if (vr.Alt.size() > 0)
				vcfRecords_.push_back(vr);
		}

		Reference_ = new char[ReferenceLength_ + 1];
		Reference_[ReferenceLength_] = '\0';
		memcpy(Reference_, aReference.data(), ReferenceLength_*sizeof(char));

		size_t lastPos = 0;
		for (auto & record : vcfRecords_) {
			switch (record.Type) {
				case vcfrtSNP: ++SNPs_; break;
				case vcfrtInsertion: ++Insertions_; break;
				case vcfrtDeletion: ++Deletions_; break;
				case vcfrtReplace: ++Replaces_; break;
			}

			if (lastPos >= record.Pos - 1)
				continue;

			AlternateSize_ += (record.Pos - 1 - lastPos);
			AlternateSize_ += record.Alt.size();
			lastPos = record.Pos - 1 + record.Ref.size();
		}

		Alternate_ = new char[AlternateSize_ + 1];
		Alternate_[AlternateSize_] = '\0';
		lastPos = 0;
		size_t altPos = 0;
		for (auto & record : vcfRecords_) {
			if (lastPos >= record.Pos - 1)
				continue;

			const size_t refSeqPartSize = (record.Pos - 1 - lastPos);

			Intervals_.push_back(CRefSeqInterval(Alternate_ + altPos, lastPos, refSeqPartSize));
			memcpy(Alternate_ + altPos, Reference_ + lastPos, refSeqPartSize*sizeof(char));
			altPos += refSeqPartSize;
			lastPos = record.Pos -1 + record.Ref.size();
			Intervals_.push_back(CRefSeqInterval(Alternate_ + altPos, record.Pos - 1, record.Ref.size()));
			memcpy(Alternate_ + altPos, record.Alt.data(), record.Alt.size()*sizeof(char));
			altPos += record.Alt.size();
		}

		return;
	}
	~CUpdatedRefSeq()
	{
		delete[] Alternate_;
		delete[] Reference_;

		return;
	}
	const char *Reference() const { return Reference_; }
	const char *Alternate() const { return Alternate_; }
	std::size_t Insertions() const { return Insertions_; }
	std::size_t Deletions() const { return Deletions_; }
	std::size_t SNPs() const { return SNPs_; }
	const std::size_t & Replaces = Replaces_;
	void compare(const CUpdatedRefSeq & aOther, CVCFCompareStatistics & aStatistics)
	{
		auto & ointer = aOther.Intervals();
		size_t matchCount = 0;

		for (auto & record : vcfRecords_) {			
			const char *altSeq = NULL;
			size_t recordPos = record.Pos - 1;

			for (size_t index = 0; index < ointer.size(); ++index) {
				if (ointer[index].Pos() <= recordPos && recordPos < ointer[index].Pos() + ointer[index].Length())
					altSeq = ointer[index].Sequence() + recordPos - ointer[index].Pos();
				else if (index < ointer.size() && recordPos >= ointer[index].Pos() + ointer[index].Length() && recordPos < ointer[index + 1].Pos())
					altSeq = ointer[index].Sequence() + recordPos - (ointer[index].Pos() + ointer[index].Length());
			
				if (altSeq != NULL)
					break;
			}

			if (altSeq != NULL) {
				if (memcmp(altSeq, record.Alt.data(), (record.Alt.size())*sizeof(char)) == 0) {
					switch (record.Type) {
						case vcfrtDeletion: {
							const char *refStart = Reference_ + record.Pos - 1 + 1;

							assert(record.Alt.size() == 1);
							if (memcmp(refStart, altSeq + 1, 5 * sizeof(char)) == 0)
								++aStatistics.Deletions;
						} break;
						case vcfrtInsertion:
							++aStatistics.Insertions;
							break;
						case vcfrtReplace:
							++aStatistics.Replaces;
							break;
						case vcfrtSNP:
							++aStatistics.SNPs;
							break;
					}
				} else {
					const char *refStart = Reference_ + record.Pos - 1;
					std::string refSeqpart(refStart, refStart + record.Alt.size());
					std::string alt(altSeq, altSeq + record.Alt.size());

					if (!alt.compare(refSeqpart))
						++aStatistics.Undiscovered;

					aStatistics.Differences.push_back(CVCFDifference(record.Pos, refSeqpart, record.Alt, alt));
				}
			}
		}

		return;
	}

	const std::vector<CVCFRecord> & VCFRecords() const { return vcfRecords_; }
	const std::vector<CRefSeqInterval> & Intervals() const { return Intervals_;  }
private:
	std::size_t Insertions_;
	std::size_t Deletions_;
	std::size_t SNPs_;
	std::size_t Replaces_;
	std::vector<CVCFRecord> vcfRecords_;
	char *Reference_;
	std::size_t ReferenceLength_;
	char *Alternate_;
	size_t AlternateSize_;
	std::vector<CRefSeqInterval> Intervals_;
};




#endif
