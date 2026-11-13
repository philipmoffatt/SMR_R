// Optional Windows/MinGW scalar exp for the frozen GRASS CPU math profile.
// Build with -fno-fast-math -ffp-contract=off. No GRASS runtime is required.
#include <cmath>
#include <cstdint>
#include "xla/ffi/api/ffi.h"

#if !defined(__MINGW64__) || !defined(__x86_64__)
#error "The frozen CPU profile is qualified only for Windows x64 MinGW"
#endif

namespace ffi = xla::ffi;

// CPython starts x87 at 53-bit precision; native MinGW GRASS uses 64-bit
// precision internally in exp. Control is thread-local and must be restored
// before returning to XLA, including when a worker is reused by another call.
class NativePrecision {
 public:
  NativePrecision() {
    __asm__ volatile("fnstcw %0" : "=m"(previous_) : : "memory");
    unsigned short native = (previous_ & ~0x0300) | 0x0300;
    __asm__ volatile("fldcw %0" : : "m"(native) : "memory");
  }
  ~NativePrecision() {
    __asm__ volatile("fldcw %0" : : "m"(previous_) : "memory");
  }
 private:
  unsigned short previous_;
};

ffi::Error ScalarExp(ffi::Buffer<ffi::F64> input,
                     ffi::ResultBuffer<ffi::F64> output) {
  if (input.element_count() != output->element_count())
    return ffi::Error::InvalidArgument("exp input/output sizes differ");
  NativePrecision precision;
  const double* x = input.typed_data();
  double* y = output->typed_data();
  for (int64_t i = 0; i < input.element_count(); ++i) y[i] = std::exp(x[i]);
  return ffi::Error::Success();
}

XLA_FFI_DEFINE_HANDLER_SYMBOL(
    SmrLegacyExp, ScalarExp,
    ffi::Ffi::Bind().Arg<ffi::Buffer<ffi::F64>>().Ret<ffi::Buffer<ffi::F64>>());
