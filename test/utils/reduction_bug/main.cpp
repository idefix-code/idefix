#include <Kokkos_Core.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

int main(int argc, char *argv[]) {

  Kokkos::initialize(argc, argv);

  {
    std::cout << "Execution space: "
              << Kokkos::DefaultExecutionSpace::name()
              << "\n\n";

    constexpr int nRepeat = 5;

    std::vector<int> sizes = {
        1,     2,     3,     7,
        15,    31,    32,    33,
        63,    64,    65,
        127,   128,   129,
        255,   256,   257,
        511,   512,   513,
        1023,  1024,  1025,
        2047,  2048,  2049,
        4095,  4096,  4097,
        8191,  8192,  8193,
        16384,
        32768,
        65535, 65536, 65537,
        100000};

    const std::vector<int> moduli = {
        32, 64, 128, 256, 512, 1024};

    for (int N : sizes) {

      std::vector<int> failed;

      std::cout << "=========================================\n";
      std::cout << "Testing N = " << N << "\n";

      for (int target = 0; target < N; target++) {

        bool success = true;

        for (int rep = 0; rep < nRepeat; rep++) {

          double result = std::numeric_limits<double>::max();

          Kokkos::parallel_reduce(
              "MinTest",
              Kokkos::RangePolicy<>(0, N),

              KOKKOS_LAMBDA(const int idx, double &minval) {

                const double val =
                    (idx == target) ? 5.0e-3 : 1.0e-2;

                if (val < minval)
                  minval = val;

              },

              Kokkos::Min<double>(result));

          if (std::abs(result - 5.0e-3) > 1e-12) {
            success = false;
            break;
          }
        }

        if (!success)
          failed.push_back(target);
      }

      if (failed.empty()) {
        std::cout << "PASS\n\n";
        continue;
      }

      std::cout << "FAIL : "
                << failed.size()
                << " failing indices\n\n";

      std::cout << "First failing indices:\n";
      for (size_t i = 0;
           i < std::min<size_t>(20, failed.size());
           i++)
        std::cout << failed[i] << " ";
      std::cout << "\n\n";

      std::cout << "Last failing indices:\n";
      size_t first =
          failed.size() > 20 ? failed.size() - 20 : 0;

      for (size_t i = first; i < failed.size(); i++)
        std::cout << failed[i] << " ";
      std::cout << "\n\n";

      //------------------------------------------------------------------
      // Analyse residues modulo common reduction granularities
      //------------------------------------------------------------------

      for (int m : moduli) {

        std::vector<int> hist(m, 0);

        for (int idx : failed)
          hist[idx % m]++;

        int maxCount =
            *std::max_element(hist.begin(), hist.end());

        if (maxCount == 0)
          continue;

        std::cout << "Modulo " << m
                  << " residue histogram:\n";

        bool interesting = false;

        for (int r = 0; r < m; r++) {

          if (hist[r] > 0) {

            interesting = true;

            std::cout << "  r = "
                      << std::setw(3) << r
                      << " : "
                      << hist[r]
                      << "\n";
          }
        }

        if (!interesting)
          std::cout << "  none\n";

        std::cout << "\n";
      }
    }
  }

  Kokkos::finalize();

  return 0;
}
