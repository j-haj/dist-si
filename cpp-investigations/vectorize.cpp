#include <cctype>
#include <chrono>
#include <iostream>
#include <random>
#include <vector>

#include <omp.h>

using namespace std;

const size_t DIM = 5000;

template <typename T>
void  mul(vector<T>& x, const vector<T>& v1,
	  const vector<T>& v2, const vector<T>& v3) {
  for (std::size_t i = 0; i < DIM; ++i) {
    x[i] = x[i] + 2 * (v1[i] - v2[i]) + 2 * (v1[i] - v3[i]);
  }
}

template <typename T>
void nested_loop(vector<vector<T>>& xs, const vector<vector<T>>& v1s,
		 const vector<vector<T>>& v2s, const vector<vector<T>>& v3s) {
#pragma omp parallel for
  for (std::size_t i = 0; i < xs.size(); ++i) {
    for (std::size_t j = 0; j < DIM; ++j) {
      xs[i][j] = xs[i][j] + 2 * (v1s[i][j] - v2s[i][j]) +
	2 * (v1s[i][j] - v3s[i][j]);
    }
  }
}

template <typename T>
void single_loop(vector<vector<T>>& xs, const vector<vector<T>>& v1s,
		 const vector<vector<T>>& v2s, const vector<vector<T>>& v3s) {
#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < xs.size(); ++i) {
    mul(xs[i], v1s[i], v2s[i], v3s[i]);
  }
}

template <typename T>
vector<T> random_vec() {
  vector<T> v(DIM);
  std::random_device rd;
  static std::mt19937_64 gen(rd());
  std::uniform_real_distribution<T> dist(-10, 10);
  for (std::size_t i = 0; i < DIM; ++i) {
    v[i] = dist(gen);
  }
  return v;
}

template <typename T>
vector<vector<T>> random_vectors(std::size_t n) {
  vector<vector<T>> vs(n);
  for (std::size_t i = 0; i < n; ++i) {
    vs[i] = random_vec<T>();
  }
  return vs;
}

template <typename T>
void run_experiment(vector<vector<T>>& xs, const vector<vector<T>>& v1s,
		    const vector<vector<T>>& v2s, const vector<vector<T>>& v3s) {

  using sc = std::chrono::steady_clock;
  const int N = 10;
  vector<double> single_times;
  vector<double> nested_times;
  for (int i = 0; i < N; ++i) {
    auto start = sc::now();
    single_loop(xs, v1s, v2s, v3s);
    auto stop = sc::now();
    std::chrono::duration<double> elapsed = stop - start;
    single_times.push_back(elapsed.count());

    start = sc::now();
    nested_loop(xs, v1s, v2s, v3s);
    stop = sc::now();
    elapsed = stop - start;
    nested_times.push_back(elapsed.count());

    std::cout << "single: " << single_times[i] << "\tnested: "
	      << nested_times[i]
	      << '\n';
  }
  auto single_avg = std::accumulate(single_times.begin(),
				    single_times.end(), 0.0) / N;
  auto nested_avg = std::accumulate(nested_times.begin(),
				    nested_times.end(), 0.0) / N;
  std::cout << "single avg: " << single_avg << "\tnested avg: "
	    << nested_avg << '\n';
}

int main() {
  using T = double;
  std::size_t N = 100000;
  auto xs1 = random_vectors<T>(N);
  auto xs2(xs1);
  auto v1s = random_vectors<T>(N);
  auto v2s = random_vectors<T>(N);
  auto v3s = random_vectors<T>(N);


  omp_set_num_threads(12);
  auto n_threads = omp_get_num_threads();
  std::cout << "Running on " << n_threads << " threads\n";
  run_experiment<T>(xs1, v1s, v2s, v3s);
  
  return 0;
}
