//
// Created by hlf on 6/13/22.
//

#include "VCFWriter.h"
#include <utility>
#include "VCFConstants.h"

VCFWriter::VCFWriter(const std::string &outPath, SAMSequenceDictionary sequenceDictionary) : outPath(outPath),
                                                                                             refDict(sequenceDictionary) {
	assert(!outPath.empty());

	if (outPath.length() < 4 || outPath.rfind(".vcf") != (outPath.length() - 4))
		throw std::invalid_argument("something wrong with output file name.");

	if (sequenceDictionary.getSequences().empty())
		throw std::invalid_argument("A reference dictionary is required for creating Tribble indices on the fly");

	char outPath_cstr[outPath.length() + 1];
	strcpy(outPath_cstr, outPath.c_str());
	outFile = hts_open(outPath_cstr, "w");
	if (outFile == nullptr)
		throw std::invalid_argument("something wrong when opening output file.");

	hdr = bcf_hdr_init("w");
	if (hdr == nullptr)
		throw std::invalid_argument("something wrong when initializing vcf header.");
}

void VCFWriter::appendHeader_c(char str[]) {
	if (bcf_hdr_append(hdr, str) < 0)
		throw std::invalid_argument("something wrong when appending vcf header.");
}

void VCFWriter::appendHeader(const std::string &str) {
	char c_str[str.length() + 1];
	strcpy(c_str, str.c_str());
	appendHeader_c(c_str);
}

std::string VCFWriter::getCommandLine(int argc, char *argv[]) {
	std::string commandLine;
	for (int i = 0; i < argc; ++i) {
		if (i != 0)
			commandLine.append(" ");
		commandLine.append(argv[i]);
	}
	return commandLine;
}


void VCFWriter::add(std::shared_ptr<VariantContext>& vc) {
    bcf1_t * hts_vc = bcf_init();
    hts_vc->pos = vc->getStart();
    hts_vc->rid = vc->getContigInt();
    hts_vc->n_sample = bcf_hdr_nsamples(hdr);
    bcf_float_set_missing(hts_vc->qual);
    hts_vc->rlen = vc->getReference()->getLength();
    int len = 0;
    len += (vc->getReference()->getLength()+1);
    for(auto allele : vc->getAlternateAlleles()) {
        len += (allele->getLength()+1);
    }
    char * tmp = new char [len]{0};
    int sum = 0;
    memcpy(tmp, vc->getReference()->getBases().get(), vc->getReference()->getLength());
    sum += vc->getReference()->getLength();
    tmp[sum] = ',';
    sum++;
    for(auto allele : vc->getAlternateAlleles()) {
        memcpy(tmp + sum, allele->getBases().get(), allele->getLength());
        sum += allele->getLength();
        tmp[sum] = ',';
        sum++;
    }
    tmp[sum-1] = '\0';
    bcf_update_alleles_str(hdr, hts_vc, tmp);

    //info
    for(auto & kv : vc->getAttributes()) {
        if(kv.first == VCFConstants::REPEAT_UNIT_KEY) {
            string x = kv.second.getAttributeAsString();
            bcf_update_info_string(hdr, hts_vc, VCFConstants::REPEAT_UNIT_KEY.c_str(), x.c_str());
            continue;
        }
        if(kv.first == VCFConstants::STR_PRESENT_KEY) {
            bool x = kv.second.getAttributeAsBool();
            string s = "False";
            if(x) {
                s = "True";
            }
            bcf_update_info_string(hdr, hts_vc, VCFConstants::STR_PRESENT_KEY.c_str(), s.c_str());
            continue;
        }
        if(kv.first == VCFConstants::DEPTH_KEY) {
            int x = kv.second.getAttributeAsInt();
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::DEPTH_KEY.c_str(), &x, 1);
            continue;
        }
        if(kv.first == VCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY) {
            int x = kv.second.getAttributeAsInt();
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::EVENT_COUNT_IN_HAPLOTYPE_KEY.c_str(), &x, 1);
            continue;
        }
        if(kv.first == VCFConstants::MEDIAN_BASE_QUALITY_KEY) {
            std::vector<int> x = kv.second.getAttributeAsIntVector();
            int list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::MEDIAN_BASE_QUALITY_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY) {
            std::vector<int> x = kv.second.getAttributeAsIntVector();
            int list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::MEDIAN_FRAGMENT_LENGTH_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::MEDIAN_MAPPING_QUALITY_KEY) {
            std::vector<int> x = kv.second.getAttributeAsIntVector();
            int list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::MEDIAN_MAPPING_QUALITY_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::REPEATS_PER_ALLELE_KEY) {
            std::vector<int> x = kv.second.getAttributeAsIntVector();
            int list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::REPEATS_PER_ALLELE_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::MEDIAN_READ_POSITON_KEY) {
            std::vector<int> x = kv.second.getAttributeAsIntVector();
            int list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_int32(hdr, hts_vc, VCFConstants::MEDIAN_READ_POSITON_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::NORMAL_ARTIFACT_LOG_10_ODDS_KEY) {
            std::vector<double> x = kv.second.getAttributeAsDoubleVector();
            float list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_float(hdr, hts_vc, VCFConstants::NORMAL_ARTIFACT_LOG_10_ODDS_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::NORMAL_LOG_10_ODDS_KEY) {
            std::vector<double> x = kv.second.getAttributeAsDoubleVector();
            float list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_float(hdr, hts_vc, VCFConstants::NORMAL_LOG_10_ODDS_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::POPULATION_AF_KEY) {
            std::vector<double> x = kv.second.getAttributeAsDoubleVector();
            float list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_float(hdr, hts_vc, VCFConstants::POPULATION_AF_KEY.c_str(), list, x.size());
            continue;
        }
        if(kv.first == VCFConstants::TUMOR_LOG_10_ODDS_KEY) {
            std::vector<double> x = kv.second.getAttributeAsDoubleVector();
            float list[x.size()];
            for(int i = 0; i < x.size(); i++) {
                list[i] = x[i];
            }
            bcf_update_info_float(hdr, hts_vc, VCFConstants::TUMOR_LOG_10_ODDS_KEY.c_str(), list, x.size());
            continue;
        }
    }

    unordered_map<string, std::shared_ptr<Genotype>> sample2genotype;
    for(auto & geno : *vc->getGenotypes()->getGenotypes()) {
        sample2genotype.insert({geno->getSampleName(), geno});
    }

    //format:AD
    int maxSize = -1;
    for(string & str : sampleNamesInOrder) {
		std::vector<int> AD = sample2genotype.at(str)->getAD();
        maxSize = max(maxSize, (int) AD.size());
    }
    int ADs[maxSize * sample2genotype.size()];
    int j = 0;
    for(string &str : sampleNamesInOrder) {
	    std::vector<int> AD = sample2genotype.at(str)->getAD();
        for(int i = 0; i < AD.size(); i++) {
            ADs[i+j] = AD[i];
        }
        for(int i = (int) AD.size(); i < maxSize; i++) {
            ADs[i+j] = bcf_int32_missing;
        }
        j+=maxSize;
    }
    bcf_update_format_int32(hdr, hts_vc, VCFConstants::GENOTYPE_ALLELE_DEPTHS.c_str(), ADs, maxSize * sample2genotype.size());

    //format:DP
    int Dps[sample2genotype.size()];
    j = 0;
    for(string & str : sampleNamesInOrder) {
        Dps[j] = sample2genotype.at(str)->getDP();
        j++;
    }
    bcf_update_format_int32(hdr, hts_vc, VCFConstants::DEPTH_KEY.c_str(), Dps, sample2genotype.size());

    //format:AF
    maxSize = -1;
    for(string & str : sampleNamesInOrder) {
        int l;
        if(sample2genotype.at(str)->hasExtendedAttribute(VCFConstants::ALLELE_FRACTION_KEY)) {
            l = (int) sample2genotype.at(str)->getExtendedAttribute(VCFConstants::ALLELE_FRACTION_KEY).getAttributeAsDoubleVector().size();
        }
        maxSize = max(maxSize, l);
    }
    float fdata[maxSize * sample2genotype.size()];
    j = 0;
    for(string & str : sampleNamesInOrder) {
        std::vector<double> x;
        if(sample2genotype.at(str)->hasExtendedAttribute(VCFConstants::ALLELE_FRACTION_KEY)) {
            x = sample2genotype.at(str)->getExtendedAttribute(VCFConstants::ALLELE_FRACTION_KEY).getAttributeAsDoubleVector();
        }
        int l = (int) x.size();
        for(int i = 0; i < l; i++) {
            fdata[i+j] = x[i];
        }
        for(int i = l; i < maxSize; i++) {
            fdata[i+j] = bcf_float_missing;
        }
        j+=maxSize;
    }
    bcf_update_format_float(hdr, hts_vc, VCFConstants::ALLELE_FRACTION_KEY.c_str(), fdata, maxSize * sample2genotype.size());

    //format:F1R2
    maxSize = -1;
    for(string & str : sampleNamesInOrder) {
        int l;
        if(sample2genotype.at(str)->hasExtendedAttribute(VCFConstants::F1R2_KEY)) {
            l = (int) sample2genotype.at(str)->getExtendedAttribute(VCFConstants::F1R2_KEY).getAttributeAsIntVector().size();
        }
        maxSize = max(maxSize, l);
    }
    int* idata = new int[maxSize * sample2genotype.size()];
    j = 0;
    for(string & str : sampleNamesInOrder) {
        std::vector<int> x;
        if(sample2genotype.at(str)->hasExtendedAttribute(VCFConstants::F1R2_KEY)) {
            x = sample2genotype.at(str)->getExtendedAttribute(VCFConstants::F1R2_KEY).getAttributeAsIntVector();
        }
        int l = (int) x.size();
        for(int i = 0; i < l; i++) {
            idata[i+j] = x[i];
        }
        for(int i = l; i < maxSize; i++) {
            idata[i+j] = bcf_int32_missing;
        }
        j+=maxSize;
    }
    bcf_update_format_int32(hdr, hts_vc, VCFConstants::F1R2_KEY.c_str(), idata, maxSize * sample2genotype.size());
    delete[] idata;

    //format:F2R1
    maxSize = -1;
    for(string & str : sampleNamesInOrder) {
        size_t l;
        if(sample2genotype.at(str)->hasExtendedAttribute(VCFConstants::F2R1_KEY)) {
            l = sample2genotype.at(str)->getExtendedAttribute(VCFConstants::F2R1_KEY).getAttributeAsIntVector().size();
        }
        maxSize = max(maxSize, (int) l);
    }
    idata = new int[maxSize * sample2genotype.size()];
    j = 0;
    for(string & str : sampleNamesInOrder) {
        std::vector<int> x;
        if(sample2genotype.at(str)->hasExtendedAttribute(VCFConstants::F2R1_KEY)) {
            x = sample2genotype.at(str)->getExtendedAttribute(VCFConstants::F2R1_KEY).getAttributeAsIntVector();
        }
        int l = (int) x.size();
        for(int i = 0; i < l; i++) {
            idata[i+j] = x[i];
        }
        for(int i = l; i < maxSize; i++) {
            idata[i+j] = bcf_int32_missing;
        }
        j+=maxSize;
    }
    bcf_update_format_int32(hdr, hts_vc, VCFConstants::F2R1_KEY.c_str(), idata, maxSize * sample2genotype.size());
    delete[] idata;

    //format:SB
    maxSize = -1;
    for(string & str : sampleNamesInOrder) {
        int l;
        if(sample2genotype.at(str)->hasExtendedAttribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY)) {
            l = (int) sample2genotype.at(str)->getExtendedAttribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY).getAttributeAsIntVector().size();
        }
        maxSize = max(maxSize, l);
    }
    idata = new int[maxSize * sample2genotype.size()];
    j = 0;
    for(string & str : sampleNamesInOrder) {
        std::vector<int> x;
        if(sample2genotype.at(str)->hasExtendedAttribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY)) {
            x = sample2genotype.at(str)->getExtendedAttribute(VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY).getAttributeAsIntVector();
        }
        int l = (int) x.size();
        for(int i = 0; i < l; i++) {
            idata[i+j] = x[i];
        }
        for(int i = l; i < maxSize; i++) {
            idata[i+j] = bcf_int32_missing;
        }
        j+=maxSize;
    }
    bcf_update_format_int32(hdr, hts_vc, VCFConstants::STRAND_BIAS_BY_SAMPLE_KEY.c_str(), idata, maxSize * sample2genotype.size());
    delete[] idata;

    //format:GT
    int gt[2] = {bcf_gt_missing, bcf_gt_missing};
    bcf_update_format_int32(hdr, hts_vc, VCFConstants::GENOTYPE_KEY.c_str(), gt, 2);

    //filter
    if(vc->getFilters().empty()) {
        bcf_add_filter(hdr, hts_vc, 0);
    } else {
        for(int i : vc->getFilters()) {
            bcf_add_filter(hdr, hts_vc, i);
        }
    }

    bcf_write(outFile, hdr, hts_vc);
    bcf_destroy(hts_vc);
    delete[] tmp;
}

void VCFWriter::close() {
	bcf_hdr_destroy(hdr);
	hts_close(outFile);
}

void
VCFWriter::writeHeader(const string &cmdLine, const vector<SAMReadGroupRecord> &readGroup, const string &normalSample) {
    //version
    std::string version = "VCFv4.2";
    bcf_hdr_set_version(hdr, version.c_str());

    //filter
    for(auto & k : mFilterMetaData) {
        int len = k.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, k.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

    //format
    for(auto & val : mFormatMetaData) {
        string format = "##" + val.second;
        int len = format.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, format.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

    //info
    for(auto & val : mInfoMetaData) {
        string info = "##" + val.second;
        int len = info.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, info.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

    //contig
    for (auto &sequence: refDict.getSequences()) {
        contigMetaData.push_back("contig=<ID=" + sequence.getSequenceName() + ",length=" + std::to_string(sequence.getSequenceLength()) +">");
	}
    for(auto & val : contigMetaData) {
        string contig = "##" + val;
        int len = contig.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, contig.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

    std::set<std::string> samples;
    for (const auto &rg: readGroup) {
        std::string sampleName = rg.getAttributesNochange()[SAMReadGroupRecord::READ_GROUP_SAMPLE_TAG];
        if (sampleName == normalSample) {
            mOtherMetaData.insert(std::make_pair("normal_sample", sampleName));
        } else {
            mOtherMetaData.insert(std::make_pair("tumor_sample", sampleName));
        }
        samples.insert(sampleName);
    }
    for(auto & tmp : samples) {
        sampleNamesInOrder.emplace_back(tmp);
    }
    std::sort(sampleNamesInOrder.begin(), sampleNamesInOrder.end());
	for (const std::string &sampleName: sampleNamesInOrder) {
		char c_str[sampleName.length() + 1];
		strcpy(c_str, sampleName.c_str());
		bcf_hdr_add_sample(hdr, c_str);
	}
    for(auto & k : mOtherMetaData) {
        std::string sample = "##" + k.first + "=" + k.second;
        int len = sample.size();
        bcf_hrec_t * bht = bcf_hdr_parse_line(hdr, sample.c_str(), &len);
        bcf_hdr_add_hrec(hdr, bht);
    }

    bcf_hdr_sync(hdr);
    vcf_hdr_write(outFile, hdr);
}
