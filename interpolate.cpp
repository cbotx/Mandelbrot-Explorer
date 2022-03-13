#include <vector>
#include <math.h>

#include "interpolate.h"

void mono_cubic_interpolate(const float* xs, const float* ys, const int length, float* mp, const int mp_sz) {
    // Get consecutive differences and slopes
    std::vector<float> dys, dxs, ms;
    for (int i = 0; i < length - 1; i++) {
        float dx = xs[i + 1] - xs[i];
        float dy = ys[i + 1] - ys[i];
        dxs.push_back(dx); dys.push_back(dy); ms.push_back(dy / dx);
    }

    // Get degree-1 coefficients
    std::vector<float> c1s = { ms[0] };
    for (int i = 0; i < dxs.size() - 1; i++) {
        float m = ms[i], mNext = ms[i + 1];
        if (m * mNext <= 0) {
            c1s.push_back(0);
        }
        else {
            float dx_ = dxs[i], dxNext = dxs[i + 1], common = dx_ + dxNext;
            c1s.push_back(3 * common / ((common + dxNext) / m + (common + dx_) / mNext));
        }
    }
    c1s.push_back(ms[ms.size() - 1]);

    // Get degree-2 and degree-3 coefficients
    std::vector<float> c2s, c3s;
    for (int i = 0; i < c1s.size() - 1; i++) {
        float c1 = c1s[i], m_ = ms[i], invDx = 1 / dxs[i], common_ = c1 + c1s[i + 1] - m_ - m_;
        c2s.push_back((m_ - c1 - common_) * invDx); c3s.push_back(common_ * invDx * invDx);
    }

    // Return interpolant function
    auto func = [&](float x) -> float {
        // The rightmost point in the dataset should give an exact result
        int i = length - 1;
        if (x == xs[i]) { return ys[i]; }

        // Search for the interval x is in, returning the corresponding y if x is one of the original xs
        int low = 0, mid, high = c3s.size() - 1;
        while (low <= high) {
            mid = (low + high) / 2;
            float xHere = xs[mid];
            if (xHere < x) { low = mid + 1; }
            else if (xHere > x) { high = mid - 1; }
            else { return ys[mid]; }
        }
        i = std::max(0, high);

        // Interpolate
        float diff = x - xs[i], diffSq = diff * diff;
        return ys[i] + c1s[i] * diff + c2s[i] * diffSq + c3s[i] * diff * diffSq;
    };
    for (int i = 0; i < mp_sz; ++i) {
        mp[i] = func(1.f * i / mp_sz);
    }
}

void cos_interpolate(const float* xs, const float* ys, const int length, float* mp, const int mp_sz) {
    int p = -1;
    for (int i = 0; i < mp_sz; ++i) {
        float x = 1.0 / mp_sz * i;
        while (p < length - 1 && x >= xs[(p + 1) % length]) ++p;
        float dx = (p == length - 1) ? 1.0 - xs[p] : xs[p + 1] - xs[p];
        float dy = ys[(p + length) % length] - ys[(p + 1) % length];
        mp[i] = (std::cos(3.1416 * (x - xs[(p + length) % length]) / dx) - 1) * dy / 2 + ys[(p + length) % length];
    }
}