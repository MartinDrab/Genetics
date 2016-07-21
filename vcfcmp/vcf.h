
#ifndef __VCF_H__
#define __VCF_H__


#include <cstdint>
#include <omp.h>
#include <string>
#include <vector>
#include <map>
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
	CVCFRecord(const std::string & aRefName, const std::string & aRef, const std::string & aAlt, const size_t aPosition, const std::size_t aIndex = 0)
		: Type_(vcfrtUnknown), Type(Type_), Ref(Ref_), RefName(RefName_), Pos(Pos_), Alt(Alt_), Index(Index_),
		Index_(aIndex), RefName_(aRefName), Ref_(aRef), Alt_(aAlt), Pos_(aPosition)
	{
		if (Ref_.size() == 1 && Alt_.size() == 1)
			Type_ = vcfrtSNP;
		else if (Ref_.size() == 1)
			Type_ = vcfrtInsertion;
		else if (Alt_.size() == 1)
			Type_ = vcfrtDeletion;
		else Type_ = vcfrtReplace;

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

enum EVCFVertexType {
	vcfvtReference,
	vcfvtVariantCall,
};

struct CVCFVertex {
public:
	EVCFVertexType Type() const { return Type_; }
	size_t Position() const { return Position_; }
	CVCFVertex *Son(const char aBase) const 
	{
		size_t index = (size_t)-1;

		switch (aBase) {
			case 'A':
				index = 0;
				break;
			case 'C':
				index = 1;
				break;
			case 'G':
				index = 2;
				break;
			case 'T':
				index = 3;
				break;
		}

		return Sons_[index];
	}
	bool AddSon(const char aBase, CVCFVertex *aVertex)
	{
		bool ret = false;
		size_t index = (size_t)-1;

		switch (aBase) {
		case 'A':
			index = 0;
			break;
		case 'C':
			index = 1;
			break;
		case 'G':
			index = 2;
			break;
		case 'T':
			index = 3;
			break;
		}

		ret = (Sons_[index] == 0 || Sons_[index] == aVertex);
		if (Sons_[index] == 0)
			Sons_[index] = aVertex;
		
		return ret;
	}
	CVCFVertex(const EVCFVertexType aType, const size_t aPosition = (size_t)-1)
		: Type_(aType), Position_(aPosition)
	{
		memset(Sons_, 0, sizeof(Sons_));
		
		return;
	}
private:
	CVCFVertex *Sons_[4];
	EVCFVertexType Type_;
	size_t Position_;
};

class CVCFGraph {
public:
	bool CheckVCFRecord(const CVCFRecord & aRecord, const std::string & aReference)
	{
		bool ret = true;
		size_t pos = aRecord.Pos - 1;
		CVCFVertex *tmp = 0;

		for (size_t i = 0; i < aRecord.Ref.size(); ++i) {
			auto v = Vertices_[pos + i];
			if (v == NULL)
				throw std::exception();

			if (tmp != 0 && tmp->Son(aRecord.Ref[i - 1]) != v)
				throw std::exception();

			if (v->Type() != vcfvtReference)
				throw std::exception();

			tmp = v;
		}

		tmp = Vertices_[pos];
		for (size_t i = 0; i < aRecord.Alt.size(); ++i) {
			tmp = tmp->Son(aRecord.Alt[i]);
			ret = (tmp != 0);
			if (!ret)
				break;
		}

		if (ret && aRecord.Type == vcfrtDeletion) {
			for (size_t i = 0; i < 10; ++i) {
				const char base = aReference[pos + aRecord.Ref.size() + i];

				if (base == 'N')
					break;

				tmp = tmp->Son(base);
				ret = (tmp != 0);
				if (!ret)
					break;
			}
		}

		return ret;
	}
	void AddVCFRecord(const CVCFRecord & aRecord, const std::string & aReference)
	{
		size_t pos = aRecord.Pos - 1;
		CVCFVertex *tmp = 0;

		for (size_t i = 0; i < aRecord.Ref.size(); ++i) {			
			auto v = Vertices_[pos + i];
			if (v == NULL)
				throw std::exception();
		
			if (tmp != 0 && tmp->Son(aRecord.Ref[i - 1]) != v)
				throw std::exception();

			if (v->Type() != vcfvtReference)
				throw std::exception();

			tmp = v;
		}

		CVCFVertex *rsVertex = Vertices_[pos];
		for (size_t i = 0; i < aRecord.Alt.size() - 1; ++i) {
			const char base = aRecord.Alt[i];

			tmp = rsVertex->Son(base);
			if (tmp == 0) {
				tmp = new CVCFVertex(vcfvtVariantCall);
				if (!rsVertex->AddSon(base, tmp))
					__debugbreak();
			}

			rsVertex = tmp;
		}

		{
			tmp = Vertices_[pos + aRecord.Ref.size()];
			if (tmp != 0) {
				size_t index = 0;
				char base = aRecord.Alt[aRecord.Alt.size() - 1];

				while (!rsVertex->AddSon(base, tmp)) {
					rsVertex = rsVertex->Son(base);
					++index;
					tmp = Vertices_[pos + aRecord.Ref.size() + index];
					base = aReference[pos + aRecord.Ref.size() + index - 1];
				}
			}
		}

		return;
	}
	void AddReference(const std::string & aReference)
	{
		int baseIndex = 0;
		CVCFVertex *prev = 0;
		Vertices_ = new CVCFVertex *[aReference.size() + 1];

		memset(Vertices_, 0, sizeof(CVCFVertex *)*aReference.size() + 1);
		Vertices_[0] = new CVCFVertex(vcfvtReference, 0);
#pragma omp parallel for shared(aReference)
		for (baseIndex = 1; baseIndex <= (int)aReference.size(); ++baseIndex) {
			if (aReference[baseIndex - 1] != 'N')
				Vertices_[baseIndex] = new CVCFVertex(vcfvtReference, baseIndex);
		}

		int pos = 0;
#pragma omp parallel for shared(aReference)
		for (pos = 1; pos <= (int)aReference.size(); ++pos) {
			const char base = aReference[pos - 1];
			
			if (base != 'N') {
				if (Vertices_[pos - 1] == NULL)
					Vertices_[pos - 1] = new CVCFVertex(vcfvtReference, pos - 1);

				Vertices_[pos - 1]->AddSon(base, Vertices_[pos]);
			}
		}

		return;
	}
private:
	CVCFVertex  **Vertices_;
};


class CVCFFile {
public:
	CVCFFile(const std::string & aFileName)
		: Records(vcfRecords_), Insertions(Insertions_), Deletions(Deletions_),
		Replaces(Replaces_), SNPs(SNPs_), Insertions_(0), Deletions_(0), Replaces_(0),
		SNPs_(0)
	{
		std::ifstream fs;
		fs.open(aFileName);
		std::string line;

		while (!fs.eof()) {
			std::getline(fs, line);
			if (line.empty() || (line.size() > 0 && line[0] == '#'))
				continue;

			ParseLine(line);
		}

		for (auto & rec : vcfRecords_) {
			switch (rec.Type) {
				case vcfrtInsertion: ++Insertions_; break;
				case vcfrtDeletion: ++Deletions_; break;
				case vcfrtReplace: ++Replaces_; break;
				case vcfrtSNP: ++SNPs_; break;
			}
		}

		return;
	}
	const std::vector<CVCFRecord> & Records = vcfRecords_;
	const std::size_t & Insertions = Insertions_;
	const std::size_t & Deletions = Deletions_;
	const std::size_t & Replaces = Replaces_;
	const std::size_t & SNPs = SNPs_;
private:
	CVCFFile() = delete;
	CVCFFile(const CVCFFile & aFile) = delete;
	CVCFFile & operator = (const CVCFFile &) = delete;
	
	std::size_t Insertions_;
	std::size_t Deletions_;
	std::size_t Replaces_;
	std::size_t SNPs_;
	std::vector<CVCFRecord> vcfRecords_;
	void ParseLine(const std::string & aLine)
	{
		size_t index = aLine.find('\t');

		std::string refName = aLine.substr(0, index);
		size_t endIndex = aLine.find('\t', index + 1);
		std::string strPos = aLine.substr(index + 1, endIndex - index - 1);
		size_t pos = (size_t)strtoull(strPos.data(), 0, 10);
		if (pos == 0)
			throw std::exception();

		index = endIndex;
		endIndex = aLine.find('\t', index + 1);
		std::string ref = aLine.substr(index + 1, endIndex - index - 1);
		if (ref == ".") {
			index = endIndex;
			endIndex = aLine.find('\t', index + 1);
			ref = aLine.substr(index + 1, endIndex - index - 1);
		}

		index = endIndex;
		endIndex = aLine.find('\t', index + 1);
		std::string altString = aLine.substr(index + 1, endIndex - index - 1);
		index = altString.find(',');
		if (index != std::string::npos) {
			while (index != std::string::npos) {
				vcfRecords_.push_back(CVCFRecord(refName, ref, altString.substr(0, index), pos));
				altString = altString.substr(index + 1);
				index = altString.find(',');
			}

			vcfRecords_.push_back(CVCFRecord(refName, ref, altString, pos));
		} else vcfRecords_.push_back(CVCFRecord(refName, ref, altString, pos));

		return;
	}
};




#endif
