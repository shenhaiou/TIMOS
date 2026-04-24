#pragma once
#include <cstdint>
#include <type_traits>

#ifdef HAS_ACCELERATE
#include <Accelerate/Accelerate.h>
#endif

namespace timos {

// SplitMix64 — simple, fast, passes statistical tests; used only for seeding.
inline uint64_t splitmix64(uint64_t& state)
{
    state += 0x9E3779B97F4A7C15ULL;
    uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

class Xoshiro256ss {
public:
    using result_type = uint64_t;

    Xoshiro256ss() { seed(0x1234567890ABCDEFULL); }
    explicit Xoshiro256ss(uint64_t seed_value) { seed(seed_value); }

    void seed(uint64_t seed_value)
    {
        uint64_t sm = seed_value;
        s_[0] = splitmix64(sm);
        s_[1] = splitmix64(sm);
        s_[2] = splitmix64(sm);
        s_[3] = splitmix64(sm);
    }

    static constexpr result_type min() { return 0; }
    static constexpr result_type max() { return UINT64_MAX; }

    result_type operator()()
    {
        const uint64_t result = rotl(s_[1] * 5ULL, 7) * 9ULL;
        const uint64_t t      = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3]  = rotl(s_[3], 45);
        return result;
    }

private:
    static inline uint64_t rotl(uint64_t x, int k)
    {
        return (x << k) | (x >> (64 - k));
    }
    uint64_t s_[4];
};

// Phase D.9: Batched RNG Pool. 
class RngPool {
public:
    static constexpr int SIZE = 4096;

    explicit RngPool(uint64_t seed_value) {
#ifdef HAS_ACCELERATE
        // Use high-performance AES-CTR algorithm
        BNNSRandomGeneratorMethod method = BNNSRandomGeneratorMethodAES_CTR;
        gen_ = BNNSCreateRandomGeneratorWithSeed(method, seed_value, nullptr);
#else
        rng_.seed(seed_value);
#endif
    }

    ~RngPool() {
#ifdef HAS_ACCELERATE
        if (gen_) BNNSDestroyRandomGenerator(gen_);
#endif
    }

    // Get -log(uniform01). Used for free-path step generation.
    inline double get_neg_log() {
        if (next_log_ >= SIZE) refill_logs();
        return logs_[next_log_++];
    }

    // Get a precomputed sin/cos pair for scattering.
    inline void get_sincos(double& s, double& c) {
        if (next_sincos_ >= SIZE) refill_sincos();
        s = sin_phi_[next_sincos_];
        c = cos_phi_[next_sincos_];
        next_sincos_++;
    }

    // Get raw uniform [0, 1).
    inline double get_uniform() {
#ifdef HAS_ACCELERATE
        // If we need a one-off uniform, we pull from the logs buffer or add a new one.
        // For simplicity in the hot path, we just pull from a dedicated uniform buffer.
        if (next_uniform_ >= SIZE) refill_uniforms();
        return uniform_buffer_[next_uniform_++];
#else
        return (rng_() >> 11) * 0x1.0p-53;
#endif
    }

    uint64_t operator()() { 
#ifdef HAS_ACCELERATE
        return 0; // Not used in Accelerate path
#else
        return rng_(); 
#endif
    }

private:
    void refill_logs() {
        double raw[SIZE];
        fill_buffer_uniform(raw);
#ifdef HAS_ACCELERATE
        const int n = SIZE;
        vvlog(logs_, raw, &n);
        for (int i = 0; i < SIZE; ++i) logs_[i] = -logs_[i];
#else
        for (int i = 0; i < SIZE; ++i) logs_[i] = -std::log(raw[i]);
#endif
        next_log_ = 0;
    }

    void refill_sincos() {
        double phi[SIZE];
        const double two_pi = 2.0 * 3.14159265358979323846;
        fill_buffer_uniform(phi);
        for (int i = 0; i < SIZE; ++i) phi[i] *= two_pi;
#ifdef HAS_ACCELERATE
        const int n = SIZE;
        vvsincos(sin_phi_, cos_phi_, phi, &n);
#else
        for (int i = 0; i < SIZE; ++i) {
            sin_phi_[i] = std::sin(phi[i]);
            cos_phi_[i] = std::cos(phi[i]);
        }
#endif
        next_sincos_ = 0;
    }

#ifdef HAS_ACCELERATE
    void refill_uniforms() {
        fill_buffer_uniform(uniform_buffer_);
        next_uniform_ = 0;
    }
#endif

    // Vectorized fill of double buffer with [0, 1)
    void fill_buffer_uniform(double* buf) {
#ifdef HAS_ACCELERATE
        // Generating the full range of bits is typically faster than a constrained range.
        int64_t bits[SIZE];
        BNNSNDArrayDescriptor desc = {};
        desc.layout = BNNSDataLayoutVector;
        desc.size[0] = SIZE;
        desc.data_type = BNNSDataTypeInt64;
        desc.data = bits;
        
        // Fill with full range of 64-bit integers.
        BNNSRandomFillUniformInt(gen_, &desc, INT64_MIN, INT64_MAX);

        // Convert to [0, 1) double: (bits >> 11) / 2^53
        // We cast to uint64_t to ensure logical shift.
        for (int i = 0; i < SIZE; ++i) {
            buf[i] = ((uint64_t)bits[i] >> 11) * 0x1.0p-53;
            if (buf[i] <= 0.0) buf[i] = 1e-300;
        }
#else
        for (int i = 0; i < SIZE; ++i) {
            buf[i] = (rng_() >> 11) * 0x1.0p-53;
            if (buf[i] <= 0.0) buf[i] = 1e-300;
        }
#endif
    }

#ifdef HAS_ACCELERATE
    BNNSRandomGenerator gen_ = nullptr;
    double uniform_buffer_[SIZE];
    int next_uniform_ = SIZE;
#else
    Xoshiro256ss rng_;
#endif
    double logs_[SIZE];
    double sin_phi_[SIZE];
    double cos_phi_[SIZE];
    int next_log_ = SIZE;
    int next_sincos_ = SIZE;
};


} // namespace timos

template<class R>
inline double uniform01(R& rng)
{
    if constexpr (std::is_same_v<R, timos::RngPool>) {
        return rng.get_uniform();
    } else {
        return (rng() >> 11) * 0x1.0p-53;
    }
}
