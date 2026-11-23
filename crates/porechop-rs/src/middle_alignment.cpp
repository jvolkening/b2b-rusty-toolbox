#include <parasail.h>

#include <algorithm>
#include <cctype>
#include <climits>
#include <cstdint>
#include <exception>
#include <memory>
#include <new>
#include <string>
#include <utility>
#include <vector>

extern "C" {
struct MiddleAlignment {
    int start;
    int end;
    double partial;
    double full;
};

int porechop_middle_align(const char *read, size_t read_len, const char *adapter,
                          size_t adapter_len, int backend, MiddleAlignment *out) {
    *out = {-1, 0, 0.0, 0.0};
    if (!read_len || !adapter_len) return 1;
    if (read_len > INT_MAX || adapter_len > INT_MAX / 4) return 0;
    try {
        auto normalize = [](char c) {
            switch (std::toupper(static_cast<unsigned char>(c))) {
                case 'A': return 'A';
                case 'C': return 'C';
                case 'G': return 'G';
                case 'T': case 'U': return 'T';
                default: return 'N';
            }
        };
        std::string r(read, read_len), a(adapter, adapter_len);
        for (char &c : r) c = normalize(c);
        for (char &c : a) c = normalize(c);

        const size_t k = std::min<size_t>(8, a.size() / 4);
        const size_t span = a.size() - k;
        const size_t offsets[] = {0, span / 3, span - span / 3 - (span % 3 != 0), span};
        std::vector<std::pair<size_t, size_t>> intervals;
        for (size_t offset : offsets) {
            const std::string seed = a.substr(offset, k);
            size_t pos = 0;
            while ((pos = r.find(seed, pos)) != std::string::npos) {
                const size_t margin = a.size() + 24;
                const size_t seed_end = pos + k;
                intervals.emplace_back(pos > margin ? pos - margin : 0,
                                       seed_end + std::min(margin, r.size() - seed_end));
                ++pos;
            }
        }
        if (intervals.empty()) return 1;
        std::sort(intervals.begin(), intervals.end());
        std::vector<std::pair<size_t, size_t>> merged;
        for (const auto &interval : intervals) {
            if (!merged.empty() && interval.first <= merged.back().second)
                merged.back().second = std::max(merged.back().second, interval.second);
            else
                merged.push_back(interval);
        }

        using Matrix = std::unique_ptr<parasail_matrix_t, decltype(&parasail_matrix_free)>;
        using Result = std::unique_ptr<parasail_result_t, decltype(&parasail_result_free)>;
        using Trace = std::unique_ptr<parasail_traceback_t, decltype(&parasail_traceback_free)>;
        thread_local Matrix matrix(parasail_matrix_create("ACGTN", 3, -6), parasail_matrix_free);
        if (!matrix) return 0;
        Result best(nullptr, parasail_result_free);
        size_t read_offset = 0, window_len = 0;
        for (const auto &interval : merged) {
            const size_t length = interval.second - interval.first;
            if (length > static_cast<size_t>(INT_MAX - 32) - a.size()
                    || length > (SIZE_MAX / 4) / (a.size() + 8)) return 0;
            parasail_result_t *raw = nullptr;
#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64) || defined(_M_IX86)
            if (backend == 2)
                raw = parasail_sg_trace_striped_avx2_256_32(
                    a.data(), static_cast<int>(a.size()), r.data() + interval.first,
                    static_cast<int>(interval.second - interval.first), 5, 2, matrix.get());
            else if (backend == 1)
                raw = parasail_sg_trace_striped_sse2_128_32(
                    a.data(), static_cast<int>(a.size()), r.data() + interval.first,
                    static_cast<int>(interval.second - interval.first), 5, 2, matrix.get());
            else
#endif
                raw = parasail_sg_trace(a.data(), static_cast<int>(a.size()),
                                        r.data() + interval.first,
                                        static_cast<int>(interval.second - interval.first),
                                        5, 2, matrix.get());
            Result candidate(raw, parasail_result_free);
            if (!candidate || parasail_result_is_saturated(candidate.get())) return 0;
            if (!best || candidate->score > best->score) {
                best = std::move(candidate);
                read_offset = interval.first;
                window_len = interval.second - interval.first;
            }
        }
        Trace trace(parasail_result_get_traceback(best.get(), a.data(), static_cast<int>(a.size()),
                                                 r.data() + read_offset, static_cast<int>(window_len),
                                                 matrix.get(), '|', ':', '.'),
                    parasail_traceback_free);
        if (!trace) return 0;
        const std::string ar(trace->query), rr(trace->ref);
        const size_t af = ar.find_first_not_of('-'), al = ar.find_last_not_of('-');
        const size_t rf = rr.find_first_not_of('-'), rl = rr.find_last_not_of('-');
        if (af == std::string::npos || rf == std::string::npos) return 1;
        const size_t begin = std::max(af, rf), end = std::min(al, rl);
        if (begin > end) return 1;
        int start = 0, pos = 0, partial_matches = 0, full_matches = 0;
        for (size_t i = 0; i < ar.size(); ++i) {
            if (i == begin) start = pos;
            if (i >= begin && i <= end && ar[i] == rr[i]) ++partial_matches;
            if (i >= af && i <= al && ar[i] == rr[i]) ++full_matches;
            if (i == end) {
                *out = {start + static_cast<int>(read_offset),
                        pos + 1 + static_cast<int>(read_offset),
                        100.0 * partial_matches / (end - begin + 1),
                        100.0 * full_matches / (al - af + 1)};
                return 1;
            }
            if (rr[i] != '-') ++pos;
        }
        return 0;
    } catch (const std::exception &) {
        return 0;
    }
}
}
