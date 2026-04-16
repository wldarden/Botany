// src/format.h — human-readable display formatting for ImGui
// Engine units: distance=dm, time=hours, mass=grams glucose
// These convert to natural units (m/cm/mm, kg/g/mg/µg) for display.
#pragma once

#include <cstdio>
#include <cmath>

namespace botany {

// Rotating buffer — each call returns a unique slot.
// Safe for up to 8 format calls in one expression.
inline char* fmt_buf() {
    static char bufs[8][32];
    static int idx = 0;
    return bufs[idx++ & 7];
}

// Distance (dm) → m / cm / mm
inline const char* fmt_dist(float dm) {
    char* buf = fmt_buf();
    float a = std::abs(dm);
    if (a < 1e-8f)        std::snprintf(buf, 32, "0");
    else if (a < 0.01f)   std::snprintf(buf, 32, "%.3f mm", dm * 100.0f);   // < 1 mm
    else if (a < 0.1f)    std::snprintf(buf, 32, "%.2f mm", dm * 100.0f);   // 1-10 mm
    else if (a < 10.0f)   std::snprintf(buf, 32, "%.2f cm", dm * 10.0f);    // 1 cm - 1 m
    else                   std::snprintf(buf, 32, "%.3f m", dm * 0.1f);      // >= 1 m
    return buf;
}

// Mass (grams glucose) → kg / g / mg / µg
inline const char* fmt_mass(float g) {
    char* buf = fmt_buf();
    float a = std::abs(g);
    if (a < 1e-8f)        std::snprintf(buf, 32, "0");
    else if (a < 0.001f)  std::snprintf(buf, 32, "%.1f \xC2\xB5g", g * 1e6f);  // µg
    else if (a < 1.0f)    std::snprintf(buf, 32, "%.1f mg", g * 1e3f);
    else if (a < 1000.0f) std::snprintf(buf, 32, "%.2f g", g);
    else                   std::snprintf(buf, 32, "%.2f kg", g * 1e-3f);
    return buf;
}

// Mass rate (g/hr) → kg/g/mg/µg per hr
inline const char* fmt_mass_rate(float g) {
    char* buf = fmt_buf();
    float a = std::abs(g);
    if (a < 1e-8f)        std::snprintf(buf, 32, "0");
    else if (a < 0.001f)  std::snprintf(buf, 32, "%.1f \xC2\xB5g/hr", g * 1e6f);
    else if (a < 1.0f)    std::snprintf(buf, 32, "%.1f mg/hr", g * 1e3f);
    else if (a < 1000.0f) std::snprintf(buf, 32, "%.2f g/hr", g);
    else                   std::snprintf(buf, 32, "%.2f kg/hr", g * 1e-3f);
    return buf;
}

// Volume (ml) → L / mL / µL
inline const char* fmt_vol(float ml) {
    char* buf = fmt_buf();
    float a = std::abs(ml);
    if (a < 1e-8f)        std::snprintf(buf, 32, "0");
    else if (a < 0.001f)  std::snprintf(buf, 32, "%.1f \xC2\xB5L", ml * 1e3f);
    else if (a < 1000.0f) std::snprintf(buf, 32, "%.2f mL", ml);
    else                   std::snprintf(buf, 32, "%.2f L", ml * 1e-3f);
    return buf;
}

// Dimensionless signaling units (auxin, cytokinin, GA, ethylene, stress)
inline const char* fmt_au(float val) {
    char* buf = fmt_buf();
    float a = std::abs(val);
    if (a < 1e-6f)      std::snprintf(buf, 32, "0");
    else if (a < 0.01f) std::snprintf(buf, 32, "%.4f", val);
    else if (a < 1.0f)  std::snprintf(buf, 32, "%.3f", val);
    else                 std::snprintf(buf, 32, "%.2f", val);
    return buf;
}

} // namespace botany
