#ifndef UTIL_H__
#define UTIL_H__

#include <cctype>
#include <random>
#include <vector>

namespace util {
template <typename T>
inline std::vector<T> uniform_rand_vec(std::size_t n, T min_x, T max_x) {
  std::vector<T> v(n);
  std::random_device rd;
  static std::mt19937_64 gen(rd());
  std::uniform_real_distribution<T> dist(min_x, max_x);
  for (std::size_t i = 0; i < n; ++i) {
    v[i] = dist(gen);
  }
  return v;
}

template <typename T>
inline T uniform_unit() {
  std::random_device rd;
  static std::mt19937_64 gen(rd());
  std::uniform_real_distribution<T> dist(0, 1);
  return dist(gen);
}

  template <typename T>
  inline std::vector<T> uniform_unit_vec(std::size_t n) {
    std::random_device rd;
    static std::mt19937_64 gen(rd());
    std::uniform_real_distribution<T> dist(0, 1);
    std::vector<T> v(n);
    for (std::size_t i = 0; i < n; ++i) {
      v[i] = dist(gen);
    }
    return v;
  }

  template <typename T>
  inline T squared_diff(const std::vector<T>& v1, const std::vector<T>& v2) {
    T diff = 0.0;
    for (std::size_t i = 0; i < v1.size(); ++i) {
      const auto d = v1[i] - v2[i];
      diff += d * d;
    }
    return diff;
  }
}

#endif  // UTIL_H__
