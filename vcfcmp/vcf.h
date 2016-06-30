
#ifndef __VCF_H__
#define __VCF_H__


#include <cstdint>
#include <string>
#include <vector>
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
		: Type(Type_)
	{ }
	void operator =(const CVCFRecord & aRecord)
	{
		Ref_ = aRecord.Ref();
		Alt_ = aRecord.Alt();
		Pos_ = aRecord.Pos();
		Type_ = aRecord.Type;
		RefName_ = aRecord.RefName();

		return;
	}
	CVCFRecord(const CVCFRecord & aRecord)
		: Pos_(aRecord.Pos()), Ref_(aRecord.Ref()), Alt_(aRecord.Alt()),
		RefName_(aRecord.RefName_), Type_(aRecord.Type), Type(Type_)
	{ }
	CVCFRecord(const std::string & aLine, const bool aSecond)
		: Type_(vcfrtUnknown), Type(Type_)
	{
		size_t index = aLine.find('\t');
		
		RefName_ = aLine.substr(0, index);
		size_t endIndex = aLine.find('\t', index + 1);
		std::string strPos = aLine.substr(index + 1, endIndex - index - 1);
		Pos_ = (size_t)strtoull(strPos.data(), 0, 10);
		if (Pos_ == 0)
			throw std::exception("POS is not a number");

		index = endIndex;
		endIndex = aLine.find('\t', index + 1);
		Ref_ = aLine.substr(index + 1, endIndex - index - 1);

		index = endIndex;
		endIndex = aLine.find('\t', index + 1);
		Alt_ = aLine.substr(index + 1, endIndex - index - 1);
		index = Alt_.find(',');
		if (index != std::string::npos) {
			if (aSecond)
				Alt_ = Alt_.substr(index + 1);
			else Alt_ = Alt_.substr(0, index);
		} else Alt_.clear();

		if (Ref_.size() == 1 && Alt_.size() == 1)
			Type_ = vcfrtSNP;
		else if (Ref_.size() == 1)
			Type_ = vcfrtInsertion;
		else if (Alt_.size() == 1)
			Type_ = vcfrtDeletion;
		else Type_ = vcfrtReplace;

		return;
	}
	const std::string RefName() const { return RefName_;  }
	const std::string Ref() const { return Ref_; }
	const std::string Alt() const { return Alt_; }
	std::size_t Pos() const { return Pos_; }
	const EVCFRecordType & Type = Type_;
private:
	std::string RefName_;
	std::string Ref_;
	std::string Alt_;
	std::size_t Pos_;
	EVCFRecordType Type_;
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
		auto fs = std::fstream(aFileName, std::fstream::in);
		std::string line;

		while (!fs.eof()) {
			std::getline(fs, line);
			if (line.empty())
				continue;

			CVCFRecord vr = CVCFRecord(line, false);
			vcfRecords_.push_back(vr);
			vr = CVCFRecord(line, true);
			if (vr.Alt().size() > 0)
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

			if (lastPos >= record.Pos() - 1)
				continue;

			if (memcmp(Reference_ + record.Pos() - 1, record.Ref().data(), record.Ref().size()*sizeof(char)) != 0)
				__debugbreak();

			AlternateSize_ += (record.Pos() - 1 - lastPos);
			AlternateSize_ += record.Alt().size();
			lastPos = record.Pos() + record.Ref().size() - 1;
		}

		Alternate_ = new char[AlternateSize_ + 1];
		Alternate_[AlternateSize_] = '\0';
		lastPos = 0;
		size_t altPos = 0;
		for (auto & record : vcfRecords_) {
			if (lastPos >= record.Pos() - 1)
				continue;

			const size_t refSeqPartSize = (record.Pos() - 1 - lastPos);

			Intervals_.push_back(CRefSeqInterval(Alternate_ + altPos, lastPos, refSeqPartSize));
			memcpy(Alternate_ + altPos, Reference_ + lastPos, refSeqPartSize*sizeof(char));
			altPos += refSeqPartSize;
			lastPos = record.Pos() + record.Ref().size() - 1;
			Intervals_.push_back(CRefSeqInterval(Alternate_ + altPos, record.Pos() - 1, record.Ref().size()));
			memcpy(Alternate_ + altPos, record.Alt().data(), record.Alt().size()*sizeof(char));
			altPos += record.Alt().size();
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
		auto ointer = aOther.Intervals();
		size_t matchCount = 0;

		for (auto & record : vcfRecords_) {			
			const char *altSeq = NULL;

			for (size_t index = 0; index < ointer.size(); ++index) {
				const size_t recPos = record.Pos() - 1;

				if (ointer[index].Pos() <= recPos && recPos < ointer[index].Pos() + ointer[index].Length())
					altSeq = ointer[index].Sequence() + recPos - ointer[index].Pos();
				else if (index < ointer.size() && recPos >= ointer[index].Pos() + ointer[index].Length() && recPos < ointer[index + 1].Pos())
					altSeq = ointer[index].Sequence() + recPos - (ointer[index].Pos() + ointer[index].Length());
			
				if (altSeq != NULL)
					break;
			}

			if (altSeq != NULL) {
				if (memcmp(altSeq, record.Alt().data(), (record.Alt().size())*sizeof(char)) == 0) {
					switch (record.Type) {
						case vcfrtDeletion:
							++aStatistics.Deletions;
							break;
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
					const char *refStart = Reference_ + record.Pos() - 1;
					std::string refSeqpart(refStart, refStart + record.Alt().size());
					std::string alt(altSeq, altSeq + record.Alt().size());

					if (alt.compare(refStart))
						++aStatistics.Undiscovered;

					aStatistics.Differences.push_back(CVCFDifference(record.Pos(), refSeqpart, record.Alt(), alt));
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
