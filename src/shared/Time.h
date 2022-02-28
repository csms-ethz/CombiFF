// Copyright 2022 Salomé Rieder, CSMS ETH Zürich

#include <chrono>

namespace combi_ff {

typedef std::chrono::high_resolution_clock Clock;
typedef std::chrono::time_point<std::chrono::system_clock> Time;

inline std::chrono::milliseconds GetDuration(const Time& start,
                                             const Time& end) {
  return std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
}

std::ostream& operator<<(std::ostream& stream,
                         const std::chrono::milliseconds& time) {
  stream << std::chrono::duration_cast<std::chrono::hours>(time).count()
         << " h "
         << std::chrono::duration_cast<std::chrono::minutes>(time).count() % 60
         << " min "
         << std::chrono::duration_cast<std::chrono::seconds>(time).count() % 60
         << " s " << time.count() % 1000 << " ms";
  return stream;
}

}  // namespace combi_ff
