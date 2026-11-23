#include "alignment.h"
#include <string>

extern "C" {
struct AlignmentResult {
    int start;
    int end;
    double partial;
    double full;
};

AlignmentResult porechop_align(const char *read, size_t read_len, const char *adapter,
                               size_t adapter_len, const int *scores) {
    if (!read_len || !adapter_len) return {-1, 0, 0.0, 0.0};
    Dna5String sequence_h = std::string(read, read_len);
    Dna5String sequence_v = std::string(adapter, adapter_len);
    Align<Dna5String, ArrayGaps> alignment;
    resize(rows(alignment), 2);
    assignSource(row(alignment, 0), sequence_h);
    assignSource(row(alignment, 1), sequence_v);
    Score<int, Simple> scoring(scores[0], scores[1], scores[3], scores[2]);
    int score = globalAlignment(alignment, scoring, AlignConfig<true, true, true, true>());
    ScoredAlignment result(alignment, read_len, adapter_len, score);
    if (result.m_readStartPos == -1) return {-1, 0, 0.0, 0.0};
    return {result.m_readStartPos, result.m_readEndPos + 1,
            result.m_alignedRegionPercentIdentity, result.m_fullAdapterPercentIdentity};
}
}
