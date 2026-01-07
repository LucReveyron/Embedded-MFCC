#include <array>
#include <cmath>
#include <cstdint>
#include <utility>

constexpr int32_t A_Q15 = int32_t(0.54f * 32768.f);
constexpr int32_t B_Q15 = int32_t(0.46f * 32768.f);

template<int N>
const int16_t hamming_value(int n)
{
    float angle = 2.0f * M_PI * n / (N - 1);
    float w = 0.54f - 0.46f * std::cos(angle);  // Hamming formula
    int32_t q15 = int32_t(w * 32767.f + 0.5f);
    return int16_t(q15);
}

template<int N, int... IDX>
const auto make_hamming_lut_impl(std::integer_sequence<int, IDX...>)
{
    return std::array<int16_t, N>{ { hamming_value<N>(IDX)... } };
}

template<int N>
inline const auto make_hamming_lut()
{
    return make_hamming_lut_impl<N>(std::make_integer_sequence<int, N>());
}
