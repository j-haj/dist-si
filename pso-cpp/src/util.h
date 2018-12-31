#ifndef UTIL_H__
#define UITL_H__

#include <cctype>
#include <vector>

namespace util {
  template <typename T>
  inline std::vector<T> uniform_rand_vec(std::size_t n, T min_x, T max_x) {
    std::vector<T> v(n);
    static std::random_device rd;
    static std::mt19937_64 gen(rd);
    std::uniform_real_distribution<T> dist(min_x, max_x);
    for (std::size_t i = 0; i < v; ++i) {
      v[i] = dist(gen);
    }
    return v;
  }
}

#endif // UTIL_H__
